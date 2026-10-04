"""Synthetic byte/lineage contracts; fake indexes do not test BWA accuracy."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from validate_reference_bundle import (  # noqa: E402
    BWA_INDEX_SUFFIXES,
    BundleError,
    fingerprint,
    inspect_fasta,
    main,
    validate_bundle,
)


class ReferenceBundleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "tiny.fa"
        self.fasta.write_bytes(b">chr1 description\nACGT\nNN\n>chr2\nacgt\n")
        self.fai = self.root / "tiny.fa.fai"
        self.dictionary = self.root / "tiny.dict"
        self.reset_sidecars()
        self.indexes = []
        for suffix in BWA_INDEX_SUFFIXES:
            path = self.root / ("tiny.fa" + suffix)
            path.write_text("synthetic index bytes " + suffix)
            self.indexes.append(path)
        self.receipt = self.root / "build.sha256"
        self.receipt.write_text(
            "".join(
                hashlib.sha256(path.read_bytes()).hexdigest() + "  " + path.name + "\n"
                for path in [self.fasta, *self.indexes]
            )
        )
        self.version = self.root / "version.txt"
        self.version.write_text("synthetic-test-tool\n")
        self.kwargs = dict(
            fasta=self.fasta, fai=self.fai, dictionary=self.dictionary, indexes=self.indexes
        )

    def reset_sidecars(self):
        contigs, rows = inspect_fasta(self.fasta)
        self.fai.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
        self.dictionary.write_text(
            "@HD\tVN:1.6\n"
            + "".join(f"@SQ\tSN:{c['name']}\tLN:{c['length']}\tM5:{c['md5']}\n" for c in contigs)
        )

    def build(self):
        return validate_bundle(
            **self.kwargs,
            receipt=self.receipt,
            tool_version=self.version,
            container="synthetic-test-container",
        )

    def test_known_fasta_offsets_and_sequence_md5(self):
        contigs, fai = inspect_fasta(self.fasta)
        self.assertEqual(fai, [["chr1", 6, 18, 4, 5], ["chr2", 4, 32, 4, 5]])
        self.assertEqual(contigs[0]["md5"], hashlib.md5(b"ACGTNN").hexdigest())
        self.assertEqual(contigs[1]["md5"], hashlib.md5(b"ACGT").hexdigest())

    def test_controlled_build_roundtrips_as_prebuilt(self):
        manifest = self.build()
        path = self.root / "bundle.json"
        path.write_text(json.dumps(manifest))
        self.assertEqual(validate_bundle(**self.kwargs, prebuilt_manifest=path), manifest)
        self.assertEqual(manifest["validation"]["result"], "passed")
        self.assertNotIn(str(self.root), json.dumps(manifest))

    def test_same_name_length_base_change_rejects_dict(self):
        self.fasta.write_bytes(self.fasta.read_bytes().replace(b"ACGT", b"TCGT", 1))
        with self.assertRaisesRegex(BundleError, "M5"):
            self.build()

    def test_changed_fasta_with_regenerated_sidecars_still_rejects_old_index(self):
        self.fasta.write_bytes(self.fasta.read_bytes().replace(b"ACGT", b"TCGT", 1))
        self.reset_sidecars()
        with self.assertRaisesRegex(BundleError, "receipt"):
            self.build()

    def test_fai_offset_mismatch(self):
        self.fai.write_text(self.fai.read_text().replace("\t18\t", "\t19\t"))
        with self.assertRaisesRegex(BundleError, "offset"):
            self.build()

    def test_missing_dictionary_m5_rejected(self):
        self.dictionary.write_text("@SQ\tSN:chr1\tLN:6\n@SQ\tSN:chr2\tLN:4\n")
        with self.assertRaisesRegex(BundleError, "missing M5"):
            self.build()

    def test_changed_index_rejected(self):
        self.indexes[0].write_text("different build")
        with self.assertRaisesRegex(BundleError, "receipt"):
            self.build()

    def test_missing_manifest_is_not_pass(self):
        with self.assertRaisesRegex(BundleError, "exactly one"):
            validate_bundle(**self.kwargs)

    def test_partial_manifest_rejected(self):
        path = self.root / "bundle.json"
        path.write_text(json.dumps({"contract": "reference_bundle_v1"}))
        with self.assertRaisesRegex(BundleError, "complete"):
            validate_bundle(**self.kwargs, prebuilt_manifest=path)

    def test_manifest_tampering_rejected(self):
        manifest = self.build()
        manifest["fai"]["sha256"] = "0" * 64
        path = self.root / "bundle.json"
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(BundleError, "fingerprint/content"):
            validate_bundle(**self.kwargs, prebuilt_manifest=path)

    def test_receipt_requires_all_five_indexes(self):
        self.kwargs["indexes"] = self.indexes[:-1]
        with self.assertRaisesRegex(BundleError, "five"):
            self.build()

    def test_crlf_and_final_line_without_newline(self):
        self.fasta.write_bytes(b">chr1\r\nACGT\r\nNN")
        contigs, rows = inspect_fasta(self.fasta)
        self.assertEqual(rows, [["chr1", 6, 7, 4, 6]])
        self.assertEqual(contigs[0]["md5"], hashlib.md5(b"ACGTNN").hexdigest())

    def test_samtools_fai_line_width_regressions(self):
        # Independently generated with samtools 1.24 (pysam 0.24.1), not
        # inspect_fasta: an unterminated single line still has width bases+1.
        cases = [
            (b">chr1\nACGT", [["chr1", 4, 6, 4, 5]]),
            (b">chr1\r\nACGT", [["chr1", 4, 7, 4, 5]]),
            (b">chr1\nACGT\n", [["chr1", 4, 6, 4, 5]]),
            (b">chr1\r\nACGT\r\n", [["chr1", 4, 7, 4, 6]]),
            (
                b">chr1\nACGT\n>chr2\nTGCA",
                [["chr1", 4, 6, 4, 5], ["chr2", 4, 17, 4, 5]],
            ),
            (
                b">chr1\r\nACGT\r\n>chr2\r\nTGCA",
                [["chr1", 4, 7, 4, 6], ["chr2", 4, 20, 4, 5]],
            ),
            (b">chr1\nACGT\nACGT", [["chr1", 8, 6, 4, 5]]),
            (b">chr1\nACGT\nNN", [["chr1", 6, 6, 4, 5]]),
            (b">chr1\r\nACGT\r\nACGT", [["chr1", 8, 7, 4, 6]]),
            (b">chr1\r\nACGT\r\nNN", [["chr1", 6, 7, 4, 6]]),
        ]
        for data, expected in cases:
            with self.subTest(data=data):
                self.fasta.write_bytes(data)
                self.assertEqual(inspect_fasta(self.fasta)[1], expected)
                self.reset_sidecars()
                self.fai.write_text("".join("\t".join(map(str, row)) + "\n" for row in expected))
                self.receipt.write_text(
                    "".join(
                        hashlib.sha256(path.read_bytes()).hexdigest() + "  " + path.name + "\n"
                        for path in [self.fasta, *self.indexes]
                    )
                )
                self.test_controlled_build_roundtrips_as_prebuilt()

    def cli_arguments(self):
        args = [
            "--fasta",
            str(self.fasta),
            "--fai",
            str(self.fai),
            "--dictionary",
            str(self.dictionary),
            "--output",
            str(self.root / "output.json"),
        ]
        for index in self.indexes:
            args += ["--bwa-index", str(index)]
        return args

    def assert_cli_rejects_without_disclosure(self, args, value):
        stderr = io.StringIO()
        stdout = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(stdout):
            self.assertEqual(main(self.cli_arguments() + args), 1)
        self.assertIn("value redacted", stderr.getvalue())
        self.assertNotIn(value, stderr.getvalue() + stdout.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())
        self.assertFalse((self.root / "output.json").exists())

    def test_generated_and_prebuilt_reject_unpublishable_provenance(self):
        valid = self.build()
        unsafe_values = [
            "/opt/private-images/bwa.sif",
            "./private-images/bwa.sif",
            "~/private-images/bwa.sif",
            "file:///opt/private-images/bwa.sif",
            "https://dummy-user:dummy-password@registry.example/bwa:1",
            "http://127.0.0.1:5000/bwa:1",
        ]
        for field in ("container", "tool_version"):
            for value in unsafe_values:
                with self.subTest(field=field, value=value):
                    container = value if field == "container" else "synthetic-test-container"
                    self.version.write_text(value if field == "tool_version" else "2.2.1")
                    self.assert_cli_rejects_without_disclosure(
                        [
                            "--receipt",
                            str(self.receipt),
                            "--tool-version",
                            str(self.version),
                            "--container",
                            container,
                        ],
                        value,
                    )
                    previous = json.loads(json.dumps(valid))
                    previous["index_build"][field] = value
                    previous["fingerprint"] = fingerprint(previous)
                    manifest = self.root / "prebuilt.json"
                    manifest.write_text(json.dumps(previous))
                    self.assert_cli_rejects_without_disclosure(
                        ["--prebuilt-manifest", str(manifest)], value
                    )

    def test_container_whitespace_is_rejected_in_both_modes(self):
        value = "registry.example/bwa:1 invalid"
        self.assert_cli_rejects_without_disclosure(
            [
                "--receipt",
                str(self.receipt),
                "--tool-version",
                str(self.version),
                "--container",
                value,
            ],
            value,
        )
        previous = self.build()
        previous["index_build"]["container"] = value
        previous["fingerprint"] = fingerprint(previous)
        manifest = self.root / "prebuilt.json"
        manifest.write_text(json.dumps(previous))
        self.assert_cli_rejects_without_disclosure(["--prebuilt-manifest", str(manifest)], value)

    def test_publishable_registry_container_roundtrips(self):
        for container in (
            "quay.io/biocontainers/bwa-mem2:2.2.1",
            "docker://registry.example/bwa@sha256:" + "a" * 64,
        ):
            with self.subTest(container=container):
                manifest = validate_bundle(
                    **self.kwargs,
                    receipt=self.receipt,
                    tool_version=self.version,
                    container=container,
                )
                path = self.root / "bundle.json"
                path.write_text(json.dumps(manifest))
                self.assertEqual(validate_bundle(**self.kwargs, prebuilt_manifest=path), manifest)

    def test_short_interior_sequence_line_rejected(self):
        self.fasta.write_bytes(b">chr1\nACGT\nNN\nACGT\n")
        with self.assertRaisesRegex(BundleError, "final sequence"):
            inspect_fasta(self.fasta)

    def test_duplicate_contig_rejected(self):
        self.fasta.write_bytes(b">chr1\nACGT\n>chr1\nACGT\n")
        with self.assertRaisesRegex(BundleError, "duplicate"):
            inspect_fasta(self.fasta)


if __name__ == "__main__":
    unittest.main()
