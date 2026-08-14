"""Unit tests for bin/reconcile_variant_type_counts.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "reconcile_variant_type_counts.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


reconcile_module = _load_module("reconcile_variant_type_counts", SCRIPT_PATH)


def _write_variant_qc_tsv(
    path: Path,
    *,
    number_of_mnps: int = 0,
    number_of_others: int = 0,
    number_of_multiallelic_sites: int = 0,
) -> None:
    """Write a minimal raw/all variant_qc.tsv fixture (see Issue #16's schema)."""
    rows = [
        "cohort_id\tstage\tvariant_type\tmetric\tvalue",
        f"cohort\traw\tall\tnumber_of_mnps\t{number_of_mnps}",
        f"cohort\traw\tall\tnumber_of_others\t{number_of_others}",
        f"cohort\traw\tall\tnumber_of_multiallelic_sites\t{number_of_multiallelic_sites}",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


class ParseVcfSitesTests(unittest.TestCase):
    """Tests for parse_vcf_sites: record counts and (contig, position) sets."""

    def test_counts_records_and_collects_positions(self) -> None:
        sites = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_clean_raw_all.vcf.gz"
        )

        self.assertEqual(sites.record_count, 2)
        self.assertEqual(sites.positions, frozenset({("chrTest", "100"), ("chrTest", "200")}))

    def test_empty_vcf_has_zero_records(self) -> None:
        sites = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_clean_raw_indel.vcf.gz"
        )

        self.assertEqual(sites.record_count, 0)
        self.assertEqual(sites.positions, frozenset())

    def test_malformed_row_raises_with_file_and_cause(self) -> None:
        import gzip

        with tempfile.TemporaryDirectory() as tmp:
            broken_path = Path(tmp) / "broken.vcf.gz"
            with gzip.open(broken_path, "wt", encoding="utf-8") as handle:
                handle.write("#CHROM\tPOS\nno_position_column\n")

            with self.assertRaises(reconcile_module.MalformedVcfError) as raised:
                reconcile_module.parse_vcf_sites(broken_path)

            self.assertIn(str(broken_path), str(raised.exception))


class ReadCrossReferencedMetricsTests(unittest.TestCase):
    """Tests for reading MNP/other/multiallelic metrics from variant_qc.tsv."""

    def test_reads_all_three_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cohort.raw.all.variant_qc.tsv"
            _write_variant_qc_tsv(
                path, number_of_mnps=1, number_of_others=2, number_of_multiallelic_sites=3
            )

            metrics = reconcile_module.read_cross_referenced_metrics(path)

            self.assertEqual(
                metrics,
                {"number_of_mnps": "1", "number_of_others": "2", "number_of_multiallelic_sites": "3"},
            )

    def test_missing_metric_raises_with_file_and_cause(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cohort.raw.all.variant_qc.tsv"
            path.write_text(
                "cohort_id\tstage\tvariant_type\tmetric\tvalue\n"
                "cohort\traw\tall\tnumber_of_mnps\t0\n",
                encoding="utf-8",
            )

            with self.assertRaises(reconcile_module.MalformedVariantQcError) as raised:
                reconcile_module.read_cross_referenced_metrics(path)

            message = str(raised.exception)
            self.assertIn(str(path), message)
            self.assertIn("number_of_others", message)


class ReconcileTests(unittest.TestCase):
    """Tests for the reconciliation arithmetic across the three scenarios."""

    def test_clean_split_has_zero_not_selected_and_no_overlap(self) -> None:
        raw_all = reconcile_module.parse_vcf_sites(FIXTURES_DIR / "reconcile_clean_raw_all.vcf.gz")
        raw_snp = reconcile_module.parse_vcf_sites(FIXTURES_DIR / "reconcile_clean_raw_snp.vcf.gz")
        raw_indel = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_clean_raw_indel.vcf.gz"
        )
        cross_referenced = {"number_of_mnps": "0", "number_of_others": "0", "number_of_multiallelic_sites": "0"}

        result = reconcile_module.reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

        self.assertEqual(result.raw_all_records, 2)
        self.assertEqual(result.raw_snp_records, 2)
        self.assertEqual(result.raw_indel_records, 0)
        self.assertEqual(result.records_not_selected, 0)
        self.assertFalse(result.records_not_selected_is_negative)
        self.assertEqual(result.snp_indel_overlap_records, 0)

    def test_excluded_record_yields_positive_not_selected_without_overlap(self) -> None:
        raw_all = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_excluded_raw_all.vcf.gz"
        )
        raw_snp = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_excluded_raw_snp.vcf.gz"
        )
        raw_indel = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_excluded_raw_indel.vcf.gz"
        )
        cross_referenced = {"number_of_mnps": "1", "number_of_others": "0", "number_of_multiallelic_sites": "0"}

        result = reconcile_module.reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

        # 3 raw/all - 1 snp - 1 indel = 1 record (the MNP-like site at
        # position 200) selected into neither type-specific VCF.
        self.assertEqual(result.records_not_selected, 1)
        self.assertFalse(result.records_not_selected_is_negative)
        self.assertEqual(result.snp_indel_overlap_records, 0)

    def test_multiallelic_overlap_yields_negative_not_selected_and_is_flagged(self) -> None:
        raw_all = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_all.vcf.gz"
        )
        raw_snp = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_snp.vcf.gz"
        )
        raw_indel = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_indel.vcf.gz"
        )
        cross_referenced = {"number_of_mnps": "0", "number_of_others": "0", "number_of_multiallelic_sites": "1"}

        result = reconcile_module.reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

        # 3 raw/all - 2 snp - 2 indel = -1: the site at position 300 is
        # counted in both type-specific selections.
        self.assertEqual(result.raw_all_records, 3)
        self.assertEqual(result.raw_snp_records, 2)
        self.assertEqual(result.raw_indel_records, 2)
        self.assertEqual(result.records_not_selected, -1)
        self.assertTrue(result.records_not_selected_is_negative)
        self.assertEqual(result.snp_indel_overlap_records, 1)

    def test_negative_value_is_not_silently_omitted_from_output(self) -> None:
        raw_all = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_all.vcf.gz"
        )
        raw_snp = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_snp.vcf.gz"
        )
        raw_indel = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_overlap_raw_indel.vcf.gz"
        )
        cross_referenced = {"number_of_mnps": "0", "number_of_others": "0", "number_of_multiallelic_sites": "1"}
        result = reconcile_module.reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

        rows = reconcile_module.build_output_rows("cohort", result)
        values_by_metric = {row[1]: row[2] for row in rows}

        self.assertEqual(values_by_metric["records_not_selected"], "-1")
        self.assertEqual(values_by_metric["records_not_selected_is_negative"], "true")
        self.assertEqual(values_by_metric["snp_indel_overlap_records"], "1")

        summary_text = reconcile_module.build_summary_text("cohort", result)
        self.assertIn("WARNING: records_not_selected is negative", summary_text)
        self.assertIn("snp_indel_overlap_records = 1", summary_text)


class OutputContractTests(unittest.TestCase):
    """Tests for output header and row structure."""

    def test_output_header(self) -> None:
        self.assertEqual(reconcile_module.OUTPUT_HEADER, ("cohort_id", "metric", "value"))

    def test_output_rows_include_cross_referenced_metrics(self) -> None:
        raw_all = reconcile_module.parse_vcf_sites(FIXTURES_DIR / "reconcile_clean_raw_all.vcf.gz")
        raw_snp = reconcile_module.parse_vcf_sites(FIXTURES_DIR / "reconcile_clean_raw_snp.vcf.gz")
        raw_indel = reconcile_module.parse_vcf_sites(
            FIXTURES_DIR / "reconcile_clean_raw_indel.vcf.gz"
        )
        cross_referenced = {"number_of_mnps": "0", "number_of_others": "0", "number_of_multiallelic_sites": "0"}
        result = reconcile_module.reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

        rows = reconcile_module.build_output_rows("cohort", result)
        metrics = [row[1] for row in rows]

        self.assertIn("raw_all_number_of_mnps", metrics)
        self.assertIn("raw_all_number_of_others", metrics)
        self.assertIn("raw_all_number_of_multiallelic_sites", metrics)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_both_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            variant_qc_path = tmp_path / "cohort.raw.all.variant_qc.tsv"
            _write_variant_qc_tsv(variant_qc_path)
            output_path = tmp_path / "variant_type_accounting.tsv"
            summary_path = tmp_path / "variant_type_accounting.summary.txt"

            exit_code = reconcile_module.main(
                [
                    "--cohort-id",
                    "cohort",
                    "--raw-all-vcf",
                    str(FIXTURES_DIR / "reconcile_clean_raw_all.vcf.gz"),
                    "--raw-snp-vcf",
                    str(FIXTURES_DIR / "reconcile_clean_raw_snp.vcf.gz"),
                    "--raw-indel-vcf",
                    str(FIXTURES_DIR / "reconcile_clean_raw_indel.vcf.gz"),
                    "--raw-all-variant-qc",
                    str(variant_qc_path),
                    "--output",
                    str(output_path),
                    "--summary-output",
                    str(summary_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_path.exists())
            self.assertTrue(summary_path.exists())

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            variant_qc_path = tmp_path / "cohort.raw.all.variant_qc.tsv"
            _write_variant_qc_tsv(variant_qc_path)
            missing_vcf = tmp_path / "does-not-exist.vcf.gz"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = reconcile_module.main(
                    [
                        "--cohort-id",
                        "cohort",
                        "--raw-all-vcf",
                        str(missing_vcf),
                        "--raw-snp-vcf",
                        str(FIXTURES_DIR / "reconcile_clean_raw_snp.vcf.gz"),
                        "--raw-indel-vcf",
                        str(FIXTURES_DIR / "reconcile_clean_raw_indel.vcf.gz"),
                        "--raw-all-variant-qc",
                        str(variant_qc_path),
                        "--output",
                        str(tmp_path / "o.tsv"),
                        "--summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_vcf), stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            variant_qc_path = tmp_path / "cohort.raw.all.variant_qc.tsv"
            _write_variant_qc_tsv(variant_qc_path)
            output_path = tmp_path / "o.tsv"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--cohort-id",
                    "cohort",
                    "--raw-all-vcf",
                    str(FIXTURES_DIR / "reconcile_overlap_raw_all.vcf.gz"),
                    "--raw-snp-vcf",
                    str(FIXTURES_DIR / "reconcile_overlap_raw_snp.vcf.gz"),
                    "--raw-indel-vcf",
                    str(FIXTURES_DIR / "reconcile_overlap_raw_indel.vcf.gz"),
                    "--raw-all-variant-qc",
                    str(variant_qc_path),
                    "--output",
                    str(output_path),
                    "--summary-output",
                    str(tmp_path / "s.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("records_not_selected\t-1", output_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
