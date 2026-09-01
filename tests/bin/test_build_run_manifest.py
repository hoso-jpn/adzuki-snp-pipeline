"""Unit tests for bin/build_run_manifest.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import re
import sys
import tempfile
import types
import unittest
from datetime import UTC, datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "build_run_manifest.py"

# Issue #42: bin/ scripts import their shared helpers as a plain sibling
# module (`from manifest_utils import ...`), which resolves on its own
# both from a source checkout and inside a Nextflow task container --
# in each case Python puts the *running script's* own directory first on
# sys.path. Loading a script by file path from a test does not go
# through that mechanism, so bin/ is put on sys.path here. This is a
# test-harness detail only: production invocations need no PYTHONPATH.
sys.path.insert(0, str(REPO_ROOT / "bin"))

RUN_ID_PATTERN = re.compile(r"^\d{8}T\d{6}Z-[0-9a-f]{8}$")

REFERENCE_CLI_ARGS = [
    "--reference-id", "GCF_000000000.1",
    "--reference-species", "Vigna angularis",
    "--reference-cultivar", "Synthetic",
    "--reference-accession", "GCF_000000000.1",
]

PARAMETER_CLI_ARGS = [
    "--sample-ploidy", "2",
    "--genomicsdb-batch-size", "50",
    "--optical-duplicate-pixel-distance", "100",
    "--enable-gs-panel",
    "--snp-filter-qd-min", "2.0",
    "--snp-filter-qual-min", "30.0",
    "--snp-filter-sor-max", "3.0",
    "--snp-filter-fs-max", "60.0",
    "--snp-filter-mq-min", "40.0",
    "--snp-filter-mq-rank-sum-min", "-12.5",
    "--snp-filter-read-pos-rank-sum-min", "-8.0",
    "--indel-filter-qd-min", "2.0",
    "--indel-filter-qual-min", "30.0",
    "--indel-filter-fs-max", "200.0",
    "--indel-filter-read-pos-rank-sum-min", "-20.0",
]

CONTAINER_CLI_ARGS = [
    "--bwa-mem2-container", "bwa-mem2:test",
    "--samtools-container", "samtools:test",
    "--gatk-container", "gatk:test",
    "--python-container", "python:test",
]


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


manifest_module = _load_module("build_run_manifest", SCRIPT_PATH)


def _write_variant_qc_tsv(path: Path) -> None:
    rows = [
        "cohort_id\tstage\ttype\tmetric\tvalue",
        "cohort\traw\tall\tnumber_of_samples\t2",
        "cohort\traw\tall\tsample_names\tsample_a,sample_b",
        "cohort\traw\tall\tcohort_total_genotypes\t4",
        "cohort\traw\tall\tcohort_missing_genotypes\t0",
        "cohort\traw\tall\tcohort_missingness_rate\t0.000000",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _write_variant_type_accounting_tsv(path: Path) -> None:
    rows = [
        "cohort_id\tmetric\tvalue",
        "cohort\traw_all_records\t2",
        "cohort\traw_snp_records\t2",
        "cohort\traw_indel_records\t0",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _write_gs_panel_manifest(path: Path) -> None:
    payload = {
        "schema_version": 1,
        "run_id": "20260101T000000Z-deadbeef",
        "cohort_id": "cohort",
        "panel_status": "empty",
        "manifest_hash": "sha256:" + "0" * 64,
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_samplesheet(path: Path, rows: list[tuple[str, str, Path, Path, str, str, str]]) -> None:
    lines = ["sample_id,read_group_id,fastq_1,fastq_2,library_id,platform,platform_unit"]
    for sample_id, read_group_id, fastq_1, fastq_2, library_id, platform, platform_unit in rows:
        lines.append(
            f"{sample_id},{read_group_id},{fastq_1},{fastq_2},{library_id},{platform},{platform_unit}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class NewRunIdTests(unittest.TestCase):
    def test_format_with_injected_time_and_suffix(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        run_id = manifest_module.new_run_id(now=moment, suffix="a1b2c3d4")
        self.assertEqual(run_id, "20260814T123456Z-a1b2c3d4")

    def test_real_call_matches_expected_pattern(self) -> None:
        self.assertRegex(manifest_module.new_run_id(), RUN_ID_PATTERN)


class UtcNowIsoTests(unittest.TestCase):
    def test_format_with_injected_time(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        self.assertEqual(manifest_module.utc_now_iso(now=moment), "2026-08-14T12:34:56Z")


class ChecksumFilesTests(unittest.TestCase):
    def test_checksums_keyed_by_filename_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.txt"
            path.write_bytes(b"hello")
            checksums = manifest_module.checksum_files([path])
        self.assertEqual(
            checksums["a.txt"],
            "sha256:2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
        )

    def test_duplicate_filename_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dir_a = Path(tmp) / "a"
            dir_b = Path(tmp) / "b"
            dir_a.mkdir()
            dir_b.mkdir()
            (dir_a / "same.txt").write_bytes(b"one")
            (dir_b / "same.txt").write_bytes(b"two")
            with self.assertRaises(ValueError):
                manifest_module.checksum_files([dir_a / "same.txt", dir_b / "same.txt"])


class ParseSamplesheetTests(unittest.TestCase):
    def test_parses_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samplesheet.csv"
            fq1, fq2 = Path(tmp) / "r1.fastq.gz", Path(tmp) / "r2.fastq.gz"
            fq1.write_bytes(b"x")
            fq2.write_bytes(b"y")
            _write_samplesheet(
                path, [("sample_a", "sample_a_L001", fq1, fq2, "lib_a", "ILLUMINA", "fc.1")]
            )
            rows = manifest_module.parse_samplesheet(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample_id"], "sample_a")

    def test_missing_column_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samplesheet.csv"
            path.write_text("sample_id,fastq_1\nsample_a,r1.fastq.gz\n", encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedSamplesheetError):
                manifest_module.parse_samplesheet(path)

    def test_empty_file_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samplesheet.csv"
            path.write_text("", encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedSamplesheetError):
                manifest_module.parse_samplesheet(path)

    def test_no_rows_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samplesheet.csv"
            path.write_text(
                "sample_id,read_group_id,fastq_1,fastq_2,library_id,platform,platform_unit\n",
                encoding="utf-8",
            )
            with self.assertRaises(manifest_module.MalformedSamplesheetError):
                manifest_module.parse_samplesheet(path)


class BuildSampleEntriesTests(unittest.TestCase):
    def test_shared_fastq_across_read_groups_is_hashed_once(self) -> None:
        # Regression guard: a sample with >1 read group sharing R1/R2 must not
        # be hashed twice just because it appears in two samplesheet rows --
        # this test uses distinct files per row (the pipeline's own contract
        # requires distinct fastq pairs per read group), but confirms the
        # per-row entries still each carry the correct, independent checksum.
        with tempfile.TemporaryDirectory() as tmp:
            fq1a, fq2a = Path(tmp) / "a_r1.fastq.gz", Path(tmp) / "a_r2.fastq.gz"
            fq1b, fq2b = Path(tmp) / "b_r1.fastq.gz", Path(tmp) / "b_r2.fastq.gz"
            for f, content in [(fq1a, b"1"), (fq2a, b"2"), (fq1b, b"3"), (fq2b, b"4")]:
                f.write_bytes(content)

            rows = [
                {
                    "sample_id": "sample_a",
                    "read_group_id": "sample_a_L001",
                    "fastq_1": str(fq1a),
                    "fastq_2": str(fq2a),
                    "library_id": "lib_a",
                    "platform": "ILLUMINA",
                    "platform_unit": "fc.1",
                },
                {
                    "sample_id": "sample_a",
                    "read_group_id": "sample_a_L002",
                    "fastq_1": str(fq1b),
                    "fastq_2": str(fq2b),
                    "library_id": "lib_a",
                    "platform": "ILLUMINA",
                    "platform_unit": "fc.2",
                },
            ]
            entries = manifest_module.build_sample_entries(rows)

        self.assertEqual(len(entries), 2)
        self.assertNotEqual(
            entries[0]["fastq_1"]["checksum"], entries[1]["fastq_1"]["checksum"]
        )

    def test_identical_fastq_path_reused_across_rows_is_hashed_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fq1, fq2 = Path(tmp) / "r1.fastq.gz", Path(tmp) / "r2.fastq.gz"
            fq1.write_bytes(b"x")
            fq2.write_bytes(b"y")
            rows = [
                {
                    "sample_id": "sample_a",
                    "read_group_id": "sample_a_L001",
                    "fastq_1": str(fq1),
                    "fastq_2": str(fq2),
                    "library_id": "lib_a",
                    "platform": "ILLUMINA",
                    "platform_unit": "fc.1",
                }
            ]
            entries = manifest_module.build_sample_entries(rows)
        self.assertEqual(
            entries[0]["fastq_1"]["checksum"],
            "sha256:2d711642b726b04401627ca9fbac32f5c8530fb1903cc4db02258717921a4881",
        )


class ReadMetricTsvTests(unittest.TestCase):
    def test_reads_all_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            _write_variant_type_accounting_tsv(path)
            metrics = manifest_module.read_metric_tsv(path, required_metrics=())
        self.assertEqual(metrics["raw_all_records"], "2")

    def test_missing_required_metric_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            path.write_text("cohort_id\tmetric\tvalue\ncohort\tsome_metric\t1\n", encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_metric_tsv(path, required_metrics=("missing_metric",))

    def test_empty_file_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            path.write_text("", encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_metric_tsv(path, required_metrics=())

    def test_row_with_too_few_fields_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            path.write_text("cohort_id\tmetric\tvalue\ncohort\tsome_metric\n", encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_metric_tsv(path, required_metrics=())


class ReadVariantQcTsvTests(unittest.TestCase):
    def test_reads_metric_from_fourth_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_qc.tsv"
            _write_variant_qc_tsv(path)
            metrics = manifest_module.read_variant_qc_tsv(path, required_metrics=())
        self.assertEqual(metrics["number_of_samples"], "2")
        self.assertEqual(metrics["sample_names"], "sample_a,sample_b")

    def test_missing_required_metric_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_qc.tsv"
            _write_variant_qc_tsv(path)
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_variant_qc_tsv(path, required_metrics=("nonexistent",))

    def test_row_with_too_few_fields_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_qc.tsv"
            path.write_text(
                "cohort_id\tstage\ttype\tmetric\tvalue\ncohort\traw\tall\n", encoding="utf-8"
            )
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_variant_qc_tsv(path, required_metrics=())


class ReadGsPanelManifestTests(unittest.TestCase):
    def test_embeds_summary_fields_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            _write_gs_panel_manifest(path)
            embedded = manifest_module.read_gs_panel_manifest(path)
        self.assertEqual(embedded["panel_status"], "empty")
        self.assertNotIn("checksums", embedded)

    def test_missing_key_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            path.write_text(json.dumps({"schema_version": 1}), encoding="utf-8")
            with self.assertRaises(manifest_module.MalformedAccountingError):
                manifest_module.read_gs_panel_manifest(path)


class CanonicalJsonHashTests(unittest.TestCase):
    def test_deterministic_and_order_sensitive_to_content(self) -> None:
        payload_a = {"b": 1, "a": 2}
        payload_b = {"a": 2, "b": 1}
        self.assertEqual(
            manifest_module.canonical_json_hash(payload_a),
            manifest_module.canonical_json_hash(payload_b),
        )

    def test_different_content_differs(self) -> None:
        self.assertNotEqual(
            manifest_module.canonical_json_hash({"a": 1}),
            manifest_module.canonical_json_hash({"a": 2}),
        )


class WriteJsonAtomicTests(unittest.TestCase):
    def test_writes_readable_json_with_trailing_newline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "out.json"
            manifest_module.write_json_atomic(path, {"a": 1})
            text = path.read_text(encoding="utf-8")
        self.assertTrue(text.endswith("\n"))
        self.assertEqual(json.loads(text), {"a": 1})

    def test_no_leftover_tmp_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "out.json"
            manifest_module.write_json_atomic(path, {"a": 1})
            leftovers = [p for p in Path(tmp).iterdir() if p.name.startswith(".")]
        self.assertEqual(leftovers, [])


class MainCliTests(unittest.TestCase):
    def _run_main(self, tmp: str, *, gs_panel: bool = True) -> tuple[int, Path]:
        tmp_path = Path(tmp)
        fq1, fq2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        fq1.write_bytes(b"x")
        fq2.write_bytes(b"y")
        samplesheet = tmp_path / "samplesheet.csv"
        _write_samplesheet(
            samplesheet, [("sample_a", "sample_a_L001", fq1, fq2, "lib_a", "ILLUMINA", "fc.1")]
        )

        fasta, fai, seq_dict = tmp_path / "ref.fa", tmp_path / "ref.fa.fai", tmp_path / "ref.dict"
        fasta.write_bytes(b"fasta")
        fai.write_bytes(b"fai")
        seq_dict.write_bytes(b"dict")

        variant_qc = tmp_path / "variant_qc.tsv"
        _write_variant_qc_tsv(variant_qc)
        variant_type_accounting = tmp_path / "variant_type_accounting.tsv"
        _write_variant_type_accounting_tsv(variant_type_accounting)

        output = tmp_path / "run_manifest.json"

        argv = [
            "--cohort-id", "cohort",
            "--pipeline-version", "0.2.0-dev",
            "--git-commit", "abc123",
            "--nextflow-version", "26.04.6",
            *CONTAINER_CLI_ARGS,
            *REFERENCE_CLI_ARGS,
            "--reference-fasta", str(fasta),
            "--reference-fai", str(fai),
            "--reference-dict", str(seq_dict),
            *PARAMETER_CLI_ARGS,
            "--samplesheet", str(samplesheet),
            "--variant-qc-tsv", str(variant_qc),
            "--variant-type-accounting-tsv", str(variant_type_accounting),
            "--output", str(output),
        ]

        if gs_panel:
            gs_panel_manifest = tmp_path / "gs_panel_manifest.json"
            _write_gs_panel_manifest(gs_panel_manifest)
            argv += ["--gs-panel-manifest", str(gs_panel_manifest)]

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            exit_code = manifest_module.main(argv)
        return exit_code, output

    def test_writes_valid_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            exit_code, output = self._run_main(tmp)
            self.assertEqual(exit_code, 0)
            payload = json.loads(output.read_text(encoding="utf-8"))

        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["cohort_id"], "cohort")
        self.assertEqual(payload["git_commit"], "abc123")
        self.assertEqual(payload["nextflow_version"], "26.04.6")
        self.assertEqual(len(payload["samples"]), 1)
        self.assertEqual(payload["samples"][0]["sample_id"], "sample_a")
        self.assertEqual(payload["cohort_accounting"]["number_of_samples"], "2")
        self.assertEqual(payload["variant_type_accounting"]["raw_snp_records"], "2")
        self.assertEqual(payload["gs_panel"]["panel_status"], "empty")
        self.assertEqual(payload["parameters"]["sample_ploidy"], 2)
        self.assertIs(payload["parameters"]["enable_gs_panel"], True)
        self.assertIn("manifest_hash", payload)

    def test_empty_git_commit_recorded_as_null(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fq1, fq2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
            fq1.write_bytes(b"x")
            fq2.write_bytes(b"y")
            samplesheet = tmp_path / "samplesheet.csv"
            _write_samplesheet(
                samplesheet, [("sample_a", "sample_a_L001", fq1, fq2, "lib_a", "ILLUMINA", "fc.1")]
            )
            fasta, fai, seq_dict = tmp_path / "ref.fa", tmp_path / "ref.fa.fai", tmp_path / "ref.dict"
            fasta.write_bytes(b"fasta")
            fai.write_bytes(b"fai")
            seq_dict.write_bytes(b"dict")
            variant_qc = tmp_path / "variant_qc.tsv"
            _write_variant_qc_tsv(variant_qc)
            variant_type_accounting = tmp_path / "variant_type_accounting.tsv"
            _write_variant_type_accounting_tsv(variant_type_accounting)
            output = tmp_path / "run_manifest.json"

            argv = [
                "--cohort-id", "cohort",
                "--pipeline-version", "0.2.0-dev",
                "--nextflow-version", "26.04.6",
                *CONTAINER_CLI_ARGS,
                *REFERENCE_CLI_ARGS,
                "--reference-fasta", str(fasta),
                "--reference-fai", str(fai),
                "--reference-dict", str(seq_dict),
                *PARAMETER_CLI_ARGS,
                "--samplesheet", str(samplesheet),
                "--variant-qc-tsv", str(variant_qc),
                "--variant-type-accounting-tsv", str(variant_type_accounting),
                "--output", str(output),
            ]
            exit_code = manifest_module.main(argv)
            payload = json.loads(output.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 0)
        self.assertIsNone(payload["git_commit"])
        self.assertIsNone(payload["gs_panel"])

    def test_missing_samplesheet_column_fails_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            samplesheet = tmp_path / "samplesheet.csv"
            samplesheet.write_text("sample_id,fastq_1\nsample_a,r1.fastq.gz\n", encoding="utf-8")
            fasta, fai, seq_dict = tmp_path / "ref.fa", tmp_path / "ref.fa.fai", tmp_path / "ref.dict"
            fasta.write_bytes(b"fasta")
            fai.write_bytes(b"fai")
            seq_dict.write_bytes(b"dict")
            variant_qc = tmp_path / "variant_qc.tsv"
            _write_variant_qc_tsv(variant_qc)
            variant_type_accounting = tmp_path / "variant_type_accounting.tsv"
            _write_variant_type_accounting_tsv(variant_type_accounting)
            output = tmp_path / "run_manifest.json"

            argv = [
                "--cohort-id", "cohort",
                "--pipeline-version", "0.2.0-dev",
                "--nextflow-version", "26.04.6",
                *CONTAINER_CLI_ARGS,
                *REFERENCE_CLI_ARGS,
                "--reference-fasta", str(fasta),
                "--reference-fai", str(fai),
                "--reference-dict", str(seq_dict),
                *PARAMETER_CLI_ARGS,
                "--samplesheet", str(samplesheet),
                "--variant-qc-tsv", str(variant_qc),
                "--variant-type-accounting-tsv", str(variant_type_accounting),
                "--output", str(output),
            ]
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(argv)

        self.assertEqual(exit_code, 1)
        self.assertIn("missing required column", stderr.getvalue())
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
