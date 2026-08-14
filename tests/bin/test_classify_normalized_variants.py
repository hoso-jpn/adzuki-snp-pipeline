"""Unit tests for bin/classify_normalized_variants.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "classify_normalized_variants.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


classify_module = _load_module("classify_normalized_variants", SCRIPT_PATH)


class ClassifyVariantTests(unittest.TestCase):
    """Tests for classify_variant: pure REF/ALT shape classification."""

    def test_single_base_ref_and_alt_is_snp(self) -> None:
        self.assertEqual(classify_module.classify_variant("A", "G"), "snp")

    def test_equal_length_multibase_is_mnp(self) -> None:
        self.assertEqual(classify_module.classify_variant("AT", "GC"), "mnp")

    def test_shorter_ref_is_indel(self) -> None:
        self.assertEqual(classify_module.classify_variant("A", "ATT"), "indel")

    def test_shorter_alt_is_indel(self) -> None:
        self.assertEqual(classify_module.classify_variant("ATT", "A"), "indel")

    def test_spanning_deletion_star_is_symbolic(self) -> None:
        self.assertEqual(classify_module.classify_variant("A", "*"), "symbolic_or_star")

    def test_angle_bracket_symbolic_allele(self) -> None:
        self.assertEqual(classify_module.classify_variant("A", "<DEL>"), "symbolic_or_star")

    def test_breakend_allele(self) -> None:
        self.assertEqual(classify_module.classify_variant("A", "A[chr2:100["), "symbolic_or_star")

    def test_still_multiallelic_alt_raises(self) -> None:
        with self.assertRaises(classify_module.MalformedVcfError):
            classify_module.classify_variant("A", "G,C")

    def test_no_alt_marker_is_not_miscounted_as_snp(self) -> None:
        # "." and a single-base REF are the same string length, so this
        # must be checked before the snp/mnp length comparison.
        self.assertEqual(classify_module.classify_variant("A", "."), "no_alt")


class ParseNormalizedVcfTests(unittest.TestCase):
    """Tests for parse_normalized_vcf: classification counts and eligibility."""

    def test_normal_mixed_fixture_counts_each_class(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_normal_mixed.vcf.gz"
        )

        self.assertEqual(result.total_input_records, 5)
        self.assertEqual(
            result.class_counts,
            {"snp": 2, "mnp": 1, "indel": 1, "symbolic_or_star": 1, "no_alt": 0},
        )
        self.assertEqual(result.duplicate_key_records, 0)
        self.assertEqual(len(result.output_records), 2)

    def test_normal_mixed_fixture_header_strips_stale_filter_line(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_normal_mixed.vcf.gz"
        )

        self.assertTrue(
            all(not line.startswith("##FILTER=") for line in result.header.meta_lines)
        )
        self.assertEqual(result.header.sample_names, ("sample_a", "sample_b"))

    def test_duplicate_key_fixture_excludes_all_colliding_occurrences(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_duplicate_key.vcf.gz"
        )

        self.assertEqual(result.total_input_records, 3)
        self.assertEqual(result.class_counts["snp"], 3)
        self.assertEqual(result.distinct_duplicate_keys, (("chrTest", "100", "A", "G"),))
        self.assertEqual(result.duplicate_key_records, 2)
        # only the unique chrTest:200 C>T record survives
        self.assertEqual(len(result.output_records), 1)
        self.assertEqual(result.output_records[0].key, ("chrTest", "200", "C", "T"))

    def test_empty_vcf_reports_all_zero(self) -> None:
        result = classify_module.parse_normalized_vcf(FIXTURES_DIR / "classify_empty.vcf.gz")

        self.assertEqual(result.total_input_records, 0)
        self.assertEqual(
            result.class_counts,
            {"snp": 0, "mnp": 0, "indel": 0, "symbolic_or_star": 0, "no_alt": 0},
        )
        self.assertEqual(result.duplicate_key_records, 0)
        self.assertEqual(result.output_records, ())
        # the sample header must survive even with zero data rows
        self.assertEqual(result.header.sample_names, ("sample_a", "sample_b"))

    def test_no_alt_fixture_excludes_no_alt_record_from_output(self) -> None:
        result = classify_module.parse_normalized_vcf(FIXTURES_DIR / "classify_no_alt.vcf.gz")

        self.assertEqual(result.total_input_records, 2)
        self.assertEqual(result.class_counts["snp"], 1)
        self.assertEqual(result.class_counts["no_alt"], 1)
        self.assertEqual(len(result.output_records), 1)
        self.assertEqual(result.output_records[0].key, ("chrTest", "100", "A", "G"))

    def test_still_multiallelic_row_raises_with_cause(self) -> None:
        with self.assertRaises(classify_module.MalformedVcfError) as raised:
            classify_module.parse_normalized_vcf(
                FIXTURES_DIR / "classify_still_multiallelic.vcf.gz"
            )

        self.assertIn("multiple ALT alleles", str(raised.exception))

    def test_malformed_row_raises_with_file_and_cause(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            broken_path = Path(tmp) / "broken.vcf.gz"
            with gzip.open(broken_path, "wt", encoding="utf-8") as handle:
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\n")
                handle.write("chrTest\t100\n")

            with self.assertRaises(classify_module.MalformedVcfError) as raised:
                classify_module.parse_normalized_vcf(broken_path)

            self.assertIn(str(broken_path), str(raised.exception))


class WriteVcfTests(unittest.TestCase):
    """Tests for write_vcf: FILTER reset and output shape."""

    def test_output_records_have_filter_reset_to_dot(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_normal_mixed.vcf.gz"
        )

        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "output.vcf"
            classify_module.write_vcf(output_path, result)

            data_lines = [
                line
                for line in output_path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")
            ]

        self.assertEqual(len(data_lines), 2)
        self.assertTrue(all(line.split("\t")[6] == "." for line in data_lines))

    def test_output_is_plain_text_not_gzip(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_normal_mixed.vcf.gz"
        )

        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "output.vcf"
            classify_module.write_vcf(output_path, result)

            with output_path.open("rb") as handle:
                magic = handle.read(2)

        # gzip/BGZF magic bytes are 0x1f 0x8b; plain text must not start with them
        self.assertNotEqual(magic, b"\x1f\x8b")


class OutputContractTests(unittest.TestCase):
    """Tests for output header and row structure."""

    def test_accounting_header(self) -> None:
        self.assertEqual(
            classify_module.ACCOUNTING_HEADER, ("cohort_id", "metric", "value")
        )

    def test_accounting_rows_include_every_class_and_duplicate_metrics(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_duplicate_key.vcf.gz"
        )
        rows = classify_module.build_accounting_rows("cohort", result)
        metrics = {row[1]: row[2] for row in rows}

        self.assertEqual(metrics["total_input_records"], "3")
        self.assertEqual(metrics["snp_records"], "3")
        self.assertEqual(metrics["mnp_records"], "0")
        self.assertEqual(metrics["indel_records"], "0")
        self.assertEqual(metrics["symbolic_or_star_records"], "0")
        self.assertEqual(metrics["no_alt_records"], "0")
        self.assertEqual(metrics["duplicate_key_records"], "2")
        self.assertEqual(metrics["distinct_duplicate_keys"], "1")
        self.assertEqual(metrics["output_records"], "1")

    def test_summary_text_lists_colliding_keys(self) -> None:
        result = classify_module.parse_normalized_vcf(
            FIXTURES_DIR / "classify_duplicate_key.vcf.gz"
        )
        summary = classify_module.build_summary_text("cohort", result)

        self.assertIn("chrTest:100 A>G", summary)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_all_three_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_path = tmp_path / "output.vcf"
            accounting_path = tmp_path / "accounting.tsv"
            summary_path = tmp_path / "summary.txt"

            exit_code = classify_module.main(
                [
                    "--normalized-vcf",
                    str(FIXTURES_DIR / "classify_normal_mixed.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--output",
                    str(output_path),
                    "--accounting-output",
                    str(accounting_path),
                    "--summary-output",
                    str(summary_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_path.exists())
            self.assertTrue(accounting_path.exists())
            self.assertTrue(summary_path.exists())

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_vcf = tmp_path / "does-not-exist.vcf.gz"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = classify_module.main(
                    [
                        "--normalized-vcf",
                        str(missing_vcf),
                        "--cohort-id",
                        "cohort",
                        "--output",
                        str(tmp_path / "o.vcf"),
                        "--accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_vcf), stderr.getvalue())

    def test_main_fails_clearly_for_still_multiallelic_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = classify_module.main(
                    [
                        "--normalized-vcf",
                        str(FIXTURES_DIR / "classify_still_multiallelic.vcf.gz"),
                        "--cohort-id",
                        "cohort",
                        "--output",
                        str(tmp_path / "o.vcf"),
                        "--accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("multiple ALT alleles", stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_path = tmp_path / "output.vcf"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--normalized-vcf",
                    str(FIXTURES_DIR / "classify_empty.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--output",
                    str(output_path),
                    "--accounting-output",
                    str(tmp_path / "a.tsv"),
                    "--summary-output",
                    str(tmp_path / "s.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("sample_a\tsample_b", output_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
