"""Unit tests for bin/reconcile_gs_panel_accounting.py.

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
SCRIPT_PATH = REPO_ROOT / "bin" / "reconcile_gs_panel_accounting.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


reconcile_module = _load_module("reconcile_gs_panel_accounting", SCRIPT_PATH)


class CountVcfRecordsTests(unittest.TestCase):
    """Tests for count_vcf_records."""

    def test_counts_data_rows(self) -> None:
        count = reconcile_module.count_vcf_records(
            FIXTURES_DIR / "reconcile_gs_clean_raw_all.vcf.gz"
        )
        self.assertEqual(count, 3)

    def test_empty_vcf_has_zero_records(self) -> None:
        count = reconcile_module.count_vcf_records(
            FIXTURES_DIR / "reconcile_gs_empty_pass.vcf.gz"
        )
        self.assertEqual(count, 0)


class ReadAccountingMetricTests(unittest.TestCase):
    """Tests for read_accounting_metric."""

    def test_reads_output_records_metric(self) -> None:
        value = reconcile_module.read_accounting_metric(
            FIXTURES_DIR / "reconcile_gs_clean_normalization_accounting.tsv",
            "output_records",
        )
        self.assertEqual(value, "2")

    def test_missing_metric_raises_with_cause(self) -> None:
        with self.assertRaises(reconcile_module.MalformedAccountingError) as raised:
            reconcile_module.read_accounting_metric(
                FIXTURES_DIR / "reconcile_gs_clean_normalization_accounting.tsv",
                "no_such_metric",
            )
        self.assertIn("no_such_metric", str(raised.exception))


class CountMetadataRowsTests(unittest.TestCase):
    """Tests for count_metadata_rows."""

    def test_counts_data_rows_excluding_header(self) -> None:
        count = reconcile_module.count_metadata_rows(
            FIXTURES_DIR / "reconcile_gs_clean_sample_metadata.tsv"
        )
        self.assertEqual(count, 2)

    def test_zero_data_rows_is_valid(self) -> None:
        count = reconcile_module.count_metadata_rows(
            FIXTURES_DIR / "reconcile_gs_empty_variant_metadata.tsv"
        )
        self.assertEqual(count, 0)


class ReconcileTests(unittest.TestCase):
    """Tests for reconcile: the full-lineage arithmetic."""

    def test_clean_populated_scenario(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=3,
            normalized_records=3,
            classified_biallelic_snp_records=2,
            gs_pass_records=1,
            final_matrix_variant_records=1,
            final_matrix_sample_count=2,
        )

        self.assertEqual(result.gs_hard_filter_excluded_records, 1)
        self.assertFalse(result.gs_hard_filter_excluded_is_negative)
        self.assertTrue(result.final_matrix_variant_count_matches_pass)
        self.assertEqual(result.panel_status, "populated")

    def test_empty_panel_is_not_flagged_as_an_error(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=2,
            normalized_records=2,
            classified_biallelic_snp_records=2,
            gs_pass_records=0,
            final_matrix_variant_records=0,
            final_matrix_sample_count=2,
        )

        self.assertEqual(result.panel_status, "empty")
        self.assertTrue(result.final_matrix_variant_count_matches_pass)
        self.assertEqual(result.final_matrix_sample_count, 2)

    def test_mismatched_matrix_variant_count_is_flagged(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=3,
            normalized_records=3,
            classified_biallelic_snp_records=2,
            gs_pass_records=1,
            final_matrix_variant_records=2,
            final_matrix_sample_count=2,
        )

        self.assertFalse(result.final_matrix_variant_count_matches_pass)

    def test_negative_excluded_records_is_flagged_not_hidden(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=3,
            normalized_records=3,
            classified_biallelic_snp_records=1,
            gs_pass_records=2,
            final_matrix_variant_records=2,
            final_matrix_sample_count=2,
        )

        self.assertEqual(result.gs_hard_filter_excluded_records, -1)
        self.assertTrue(result.gs_hard_filter_excluded_is_negative)


class OutputContractTests(unittest.TestCase):
    """Tests for output header, rows, and summary text."""

    def test_output_header(self) -> None:
        self.assertEqual(reconcile_module.OUTPUT_HEADER, ("cohort_id", "metric", "value"))

    def test_summary_warns_on_mismatch(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=3,
            normalized_records=3,
            classified_biallelic_snp_records=2,
            gs_pass_records=1,
            final_matrix_variant_records=2,
            final_matrix_sample_count=2,
        )
        summary = reconcile_module.build_summary_text("cohort", result)
        self.assertIn("WARNING: final_matrix_variant_records does not match", summary)

    def test_summary_explains_empty_panel_is_normal(self) -> None:
        result = reconcile_module.reconcile(
            raw_all_records=2,
            normalized_records=2,
            classified_biallelic_snp_records=2,
            gs_pass_records=0,
            final_matrix_variant_records=0,
            final_matrix_sample_count=2,
        )
        summary = reconcile_module.build_summary_text("cohort", result)
        self.assertIn("normal outcome, not an error", summary)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def _run_main(self, prefix: str, variant_metadata_override: Path | None = None) -> tuple[int, Path, Path]:
        tmp_path = Path(tempfile.mkdtemp())
        output_path = tmp_path / "accounting.tsv"
        summary_path = tmp_path / "summary.txt"

        variant_metadata = variant_metadata_override or (
            FIXTURES_DIR / f"reconcile_gs_{prefix}_variant_metadata.tsv"
        )

        exit_code = reconcile_module.main(
            [
                "--cohort-id",
                "cohort",
                "--raw-all-vcf",
                str(FIXTURES_DIR / f"reconcile_gs_{prefix}_raw_all.vcf.gz"),
                "--normalized-vcf",
                str(FIXTURES_DIR / f"reconcile_gs_{prefix}_normalized.vcf.gz"),
                "--normalization-accounting",
                str(FIXTURES_DIR / f"reconcile_gs_{prefix}_normalization_accounting.tsv"),
                "--gs-pass-vcf",
                str(FIXTURES_DIR / f"reconcile_gs_{prefix}_pass.vcf.gz"),
                "--variant-metadata",
                str(variant_metadata),
                "--sample-metadata",
                str(FIXTURES_DIR / f"reconcile_gs_{prefix}_sample_metadata.tsv"),
                "--output",
                str(output_path),
                "--summary-output",
                str(summary_path),
            ]
        )

        return exit_code, output_path, summary_path

    def test_main_succeeds_for_clean_scenario(self) -> None:
        exit_code, output_path, summary_path = self._run_main("clean")

        self.assertEqual(exit_code, 0)
        values = {
            line.split("\t")[1]: line.split("\t")[2]
            for line in output_path.read_text(encoding="utf-8").splitlines()[1:]
        }
        self.assertEqual(values["raw_all_records"], "3")
        self.assertEqual(values["classified_biallelic_snp_records"], "2")
        self.assertEqual(values["gs_pass_records"], "1")
        self.assertEqual(values["panel_status"], "populated")
        self.assertTrue(summary_path.exists())

    def test_main_reports_empty_panel_status(self) -> None:
        exit_code, output_path, _summary_path = self._run_main("empty")

        self.assertEqual(exit_code, 0)
        values = {
            line.split("\t")[1]: line.split("\t")[2]
            for line in output_path.read_text(encoding="utf-8").splitlines()[1:]
        }
        self.assertEqual(values["gs_pass_records"], "0")
        self.assertEqual(values["panel_status"], "empty")
        self.assertEqual(values["final_matrix_sample_count"], "2")

    def test_main_never_silently_hides_a_mismatch(self) -> None:
        exit_code, output_path, summary_path = self._run_main(
            "clean",
            variant_metadata_override=FIXTURES_DIR / "reconcile_gs_mismatch_variant_metadata.tsv",
        )

        self.assertEqual(exit_code, 0)
        values = {
            line.split("\t")[1]: line.split("\t")[2]
            for line in output_path.read_text(encoding="utf-8").splitlines()[1:]
        }
        self.assertEqual(values["final_matrix_variant_count_matches_pass"], "false")
        self.assertIn("WARNING", summary_path.read_text(encoding="utf-8"))

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        tmp_path = Path(tempfile.mkdtemp())
        missing_vcf = tmp_path / "does-not-exist.vcf.gz"
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr):
            exit_code = reconcile_module.main(
                [
                    "--cohort-id",
                    "cohort",
                    "--raw-all-vcf",
                    str(missing_vcf),
                    "--normalized-vcf",
                    str(FIXTURES_DIR / "reconcile_gs_clean_normalized.vcf.gz"),
                    "--normalization-accounting",
                    str(FIXTURES_DIR / "reconcile_gs_clean_normalization_accounting.tsv"),
                    "--gs-pass-vcf",
                    str(FIXTURES_DIR / "reconcile_gs_clean_pass.vcf.gz"),
                    "--variant-metadata",
                    str(FIXTURES_DIR / "reconcile_gs_clean_variant_metadata.tsv"),
                    "--sample-metadata",
                    str(FIXTURES_DIR / "reconcile_gs_clean_sample_metadata.tsv"),
                    "--output",
                    str(tmp_path / "o.tsv"),
                    "--summary-output",
                    str(tmp_path / "s.txt"),
                ]
            )

        self.assertEqual(exit_code, 1)
        self.assertIn(str(missing_vcf), stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        tmp_path = Path(tempfile.mkdtemp())
        output_path = tmp_path / "o.tsv"

        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--cohort-id",
                "cohort",
                "--raw-all-vcf",
                str(FIXTURES_DIR / "reconcile_gs_clean_raw_all.vcf.gz"),
                "--normalized-vcf",
                str(FIXTURES_DIR / "reconcile_gs_clean_normalized.vcf.gz"),
                "--normalization-accounting",
                str(FIXTURES_DIR / "reconcile_gs_clean_normalization_accounting.tsv"),
                "--gs-pass-vcf",
                str(FIXTURES_DIR / "reconcile_gs_clean_pass.vcf.gz"),
                "--variant-metadata",
                str(FIXTURES_DIR / "reconcile_gs_clean_variant_metadata.tsv"),
                "--sample-metadata",
                str(FIXTURES_DIR / "reconcile_gs_clean_sample_metadata.tsv"),
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
        self.assertIn("panel_status\tpopulated", output_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
