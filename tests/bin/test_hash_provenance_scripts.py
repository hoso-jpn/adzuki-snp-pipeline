"""Unit tests for the Issue #42 provenance-hashing bin/ scripts.

Run with: python3 -m unittest discover -s tests/bin -v

These three scripts exist so the run manifest never has to stage a whole
cohort's raw reads, reference bundle and gVCFs into one terminal task
just to write a JSON file. What they publish is deliberately narrow --
basenames, checksums and samplesheet metadata, never a path and never
sequence content -- so that is what these tests pin.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import contextlib
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# See tests/bin/test_build_run_manifest.py for why bin/ goes on sys.path
# here but needs no PYTHONPATH in production.
sys.path.insert(0, str(REPO_ROOT / "bin"))


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


hash_input_fastqs = _load_module(
    "hash_input_fastqs", REPO_ROOT / "bin" / "hash_input_fastqs.py"
)
hash_reference_bundle = _load_module(
    "hash_reference_bundle", REPO_ROOT / "bin" / "hash_reference_bundle.py"
)
hash_run_artifacts = _load_module(
    "hash_run_artifacts", REPO_ROOT / "bin" / "hash_run_artifacts.py"
)

# A real fixture read pair this repository already ships, so the expected
# checksums below are the checksums of actual pipeline input rather than
# of something invented for the test.
FIXTURE_READ_1 = REPO_ROOT / "tests" / "data" / "reads" / "sample_a_L001_R1.fastq.gz"
FIXTURE_READ_2 = REPO_ROOT / "tests" / "data" / "reads" / "sample_a_L001_R2.fastq.gz"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return f"sha256:{digest.hexdigest()}"


class HashInputFastqsTests(unittest.TestCase):
    def _row(self, **overrides: object) -> list[str]:
        arguments: dict[str, object] = {
            "rank": 3,
            "sample_id": "sample_a",
            "read_group_id": "sample_a_L001",
            "library_id": "lib_a",
            "platform": "ILLUMINA",
            "platform_unit": "flowcell1.L001.ATCACG",
            "fastq_1": FIXTURE_READ_1,
            "fastq_2": FIXTURE_READ_2,
        }
        arguments.update(overrides)
        return hash_input_fastqs.build_row(**arguments).split("\t")

    def test_checksums_match_an_independent_hash_of_the_real_fixture(self) -> None:
        fields = self._row()
        self.assertEqual(fields[7], _sha256(FIXTURE_READ_1))
        self.assertEqual(fields[9], _sha256(FIXTURE_READ_2))

    def test_row_has_the_declared_column_count(self) -> None:
        self.assertEqual(
            len(self._row()), len(hash_input_fastqs.INPUT_PROVENANCE_COLUMNS)
        )

    def test_records_basenames_never_paths(self) -> None:
        fields = self._row()
        self.assertEqual(fields[6], "sample_a_L001_R1.fastq.gz")
        self.assertEqual(fields[8], "sample_a_L001_R2.fastq.gz")
        for field in fields:
            self.assertNotIn("/", field)

    def test_rank_is_zero_padded_so_lexical_order_is_numeric_order(self) -> None:
        self.assertEqual(self._row(rank=0)[0], "00000000")
        self.assertEqual(self._row(rank=3)[0], "00000003")
        self.assertLess(self._row(rank=9)[0], self._row(rank=10)[0])

    def test_empty_platform_unit_is_allowed(self) -> None:
        # The samplesheet contract permits an empty platform_unit.
        self.assertEqual(self._row(platform_unit="")[5], "")

    def test_negative_rank_is_rejected(self) -> None:
        with self.assertRaises(hash_input_fastqs.MalformedProvenanceFieldError):
            self._row(rank=-1)

    def test_metadata_containing_a_tab_is_rejected(self) -> None:
        # A tab would silently add a column and shift every field after
        # it into the wrong position on the reading side.
        with self.assertRaises(hash_input_fastqs.MalformedProvenanceFieldError):
            self._row(sample_id="sample\ta")

    def test_metadata_containing_a_newline_is_rejected(self) -> None:
        with self.assertRaises(hash_input_fastqs.MalformedProvenanceFieldError):
            self._row(read_group_id="rg\nother")

    def test_missing_fastq_fails_with_a_clear_message(self) -> None:
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "row.tsv"
            with contextlib.redirect_stderr(stderr):
                exit_code = hash_input_fastqs.main(
                    [
                        "--rank", "0",
                        "--sample-id", "sample_a",
                        "--read-group-id", "rg",
                        "--library-id", "lib",
                        "--platform", "ILLUMINA",
                        "--platform-unit", "",
                        "--fastq-1", str(Path(tmp) / "missing_R1.fastq.gz"),
                        "--fastq-2", str(FIXTURE_READ_2),
                        "--output", str(output),
                    ]
                )
            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())
            self.assertIn("missing_R1.fastq.gz", stderr.getvalue())


class HashReferenceBundleTests(unittest.TestCase):
    def _bundle(self, directory: Path) -> dict[str, object]:
        names = [
            "synthetic.fa",
            "synthetic.fa.fai",
            "synthetic.dict",
            "synthetic.fa.0123",
            "synthetic.fa.amb",
            "synthetic.fa.ann",
            "synthetic.fa.bwt.2bit.64",
            "synthetic.fa.pac",
        ]
        for index, name in enumerate(names):
            (directory / name).write_text(f"content-{index}\n", encoding="utf-8")
        return {
            "fasta": directory / "synthetic.fa",
            "fai": directory / "synthetic.fa.fai",
            "dict_file": directory / "synthetic.dict",
            "bwa_indexes": [directory / name for name in names[3:]],
        }

    def test_records_every_reference_file_with_its_role(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = hash_reference_bundle.build_rows(**self._bundle(Path(tmp)))

            roles = [row.split("\t")[0] for row in rows]
            self.assertEqual(roles[:3], ["fasta", "fai", "dict"])
            self.assertEqual(roles[3:], ["bwa_index"] * 5)

    def test_checksums_match_an_independent_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self._bundle(Path(tmp))
            rows = hash_reference_bundle.build_rows(**bundle)

            _role, filename, checksum = rows[0].split("\t")
            self.assertEqual(filename, "synthetic.fa")
            self.assertEqual(checksum, _sha256(bundle["fasta"]))

    def test_records_basenames_never_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = hash_reference_bundle.build_rows(**self._bundle(Path(tmp)))
            for row in rows:
                self.assertNotIn(tmp, row)
                self.assertNotIn("/", row.split("\t")[1])

    def test_missing_bwa_index_is_rejected(self) -> None:
        # An empty index set means the caller did not wire the real
        # mapping input through, which would silently under-record.
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self._bundle(Path(tmp))
            bundle["bwa_indexes"] = []
            with self.assertRaises(hash_reference_bundle.MalformedReferenceBundleError):
                hash_reference_bundle.build_rows(**bundle)

    def test_duplicate_reference_filename_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self._bundle(Path(tmp))
            bundle["bwa_indexes"] = [*bundle["bwa_indexes"], bundle["fasta"]]
            with self.assertRaises(hash_reference_bundle.MalformedReferenceBundleError):
                hash_reference_bundle.build_rows(**bundle)


class HashRunArtifactsTests(unittest.TestCase):
    def test_records_filename_and_checksum_per_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "sample_a.g.vcf.gz"
            second = Path(tmp) / "cohort.raw.vcf.gz"
            first.write_bytes(b"gvcf")
            second.write_bytes(b"cohort")

            rows = hash_run_artifacts.build_rows([first, second])

            self.assertEqual(
                rows,
                [
                    f"sample_a.g.vcf.gz\t{_sha256(first)}",
                    f"cohort.raw.vcf.gz\t{_sha256(second)}",
                ],
            )

    def test_empty_group_is_rejected(self) -> None:
        with self.assertRaises(hash_run_artifacts.MalformedArtifactGroupError):
            hash_run_artifacts.build_rows([])

    def test_duplicate_basename_within_a_group_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "a"
            second = Path(tmp) / "b"
            first.mkdir()
            second.mkdir()
            (first / "cohort.raw.vcf.gz").write_bytes(b"one")
            (second / "cohort.raw.vcf.gz").write_bytes(b"two")

            with self.assertRaises(hash_run_artifacts.MalformedArtifactGroupError):
                hash_run_artifacts.build_rows(
                    [first / "cohort.raw.vcf.gz", second / "cohort.raw.vcf.gz"]
                )


if __name__ == "__main__":
    unittest.main()
