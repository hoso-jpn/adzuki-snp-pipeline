"""Unit tests for bin/build_gs_panel.py.

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
SCRIPT_PATH = REPO_ROOT / "bin" / "build_gs_panel.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build_module = _load_module("build_gs_panel", SCRIPT_PATH)


class ClassifyGenotypeTests(unittest.TestCase):
    """Tests for classify_genotype: the dosage table and every non-standard case."""

    def test_hom_ref_is_dosage_negative_one(self) -> None:
        cell = build_module.classify_genotype("0/0")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "-1")

    def test_het_is_dosage_zero(self) -> None:
        cell = build_module.classify_genotype("0/1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "0")

    def test_het_reversed_order_is_also_dosage_zero(self) -> None:
        cell = build_module.classify_genotype("1/0")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "0")

    def test_hom_alt_is_dosage_one(self) -> None:
        cell = build_module.classify_genotype("1/1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "1")

    def test_fully_missing_dot(self) -> None:
        cell = build_module.classify_genotype(".")
        self.assertEqual(cell.category, "missing")
        self.assertEqual(cell.dosage, "nan")

    def test_fully_missing_slash(self) -> None:
        cell = build_module.classify_genotype("./.")
        self.assertEqual(cell.category, "missing")

    def test_partially_missing_is_missing_not_standard(self) -> None:
        cell = build_module.classify_genotype("0/.")
        self.assertEqual(cell.category, "missing")
        self.assertEqual(cell.dosage, "nan")

    def test_phased_is_its_own_category(self) -> None:
        cell = build_module.classify_genotype("0|1")
        self.assertEqual(cell.category, "phased")
        self.assertEqual(cell.dosage, "nan")

    def test_haploid_is_non_diploid(self) -> None:
        cell = build_module.classify_genotype("0")
        self.assertEqual(cell.category, "non_diploid")

    def test_triploid_is_non_diploid(self) -> None:
        cell = build_module.classify_genotype("0/0/1")
        self.assertEqual(cell.category, "non_diploid")

    def test_non_biallelic_index_is_its_own_category(self) -> None:
        cell = build_module.classify_genotype("0/2")
        self.assertEqual(cell.category, "non_biallelic_index")


class ParseGsPassVcfTests(unittest.TestCase):
    """Tests for parse_gs_pass_vcf: sample names, records, and GT extraction."""

    def test_standard_fixture_reads_samples_and_records_in_order(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")

        self.assertEqual(vcf.sample_names, ("sample_a", "sample_b"))
        self.assertEqual(len(vcf.records), 3)
        self.assertEqual(vcf.records[0].variant_key, "chrPanel:100:A:G")
        self.assertEqual(vcf.records[0].sample_genotypes, ("0/0", "0/1"))

    def test_empty_vcf_has_zero_records_but_keeps_sample_names(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")

        self.assertEqual(vcf.sample_names, ("sample_a", "sample_b"))
        self.assertEqual(vcf.records, ())

    def test_missing_gt_format_raises_with_cause(self) -> None:
        with self.assertRaises(build_module.MalformedVcfError) as raised:
            build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_no_gt_format.vcf.gz")

        self.assertIn("FORMAT field has no GT subfield", str(raised.exception))

    def test_zero_sample_columns_raises_with_cause(self) -> None:
        with self.assertRaises(build_module.MalformedVcfError) as raised:
            build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_no_samples.vcf.gz")

        self.assertIn("expected at least 10", str(raised.exception))


class BuildMatrixTests(unittest.TestCase):
    """Tests for build_matrix_rows and write_matrix: on-disk shape and header."""

    def test_standard_fixture_matrix_rows(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_matrix_rows(vcf)

        self.assertEqual(
            rows,
            [
                ["chrPanel:100:A:G", "-1", "0"],
                ["chrPanel:200:C:T", "1", "-1"],
                ["chrPanel:300:G:A", "0", "0"],
            ],
        )

    def test_missing_fixture_matrix_rows(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_matrix_rows(vcf)

        self.assertEqual(
            rows,
            [
                ["chrPanel:100:A:G", "nan", "0"],
                ["chrPanel:200:C:T", "nan", "1"],
            ],
        )

    def test_matrix_header_is_written_even_with_zero_variants(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")

        with tempfile.TemporaryDirectory() as tmp:
            matrix_path = Path(tmp) / "matrix.tsv.gz"
            build_module.write_matrix(matrix_path, vcf)

            with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
                lines = handle.read().splitlines()

        self.assertEqual(lines, ["variant_key\tsample_a\tsample_b"])

    def test_dosage_tokens_round_trip_through_float(self) -> None:
        for token in ("-1", "0", "1", "nan"):
            # every cell token must be float()-parseable without special-casing
            value = float(token)
            if token == "nan":
                self.assertNotEqual(value, value)  # NaN != NaN
            else:
                self.assertEqual(str(int(value)), token)


class SampleMetadataTests(unittest.TestCase):
    """Tests for build_sample_metadata_rows: per-sample missing/non-standard counts."""

    def test_standard_fixture_has_zero_missing(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "0", "0.000000", "0"])
        self.assertEqual(rows[1], ["cohort", "1", "sample_b", "0", "0.000000", "0"])

    def test_missing_fixture_counts_missing_per_sample(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        # sample_a: both variants missing (./. and 0/.) -> 2/2
        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "2", "1.000000", "0"])
        # sample_b: 0/2 variants missing -> 0/2
        self.assertEqual(rows[1], ["cohort", "1", "sample_b", "0", "0.000000", "0"])

    def test_non_standard_fixture_counts_non_standard_as_missing_too(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        # sample_a: phased, haploid, triploid, non-biallelic-index -> all 4 non-standard
        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "4", "1.000000", "4"])

    def test_empty_vcf_reports_na_rate_not_zero(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0][4], "NA")


class VariantMetadataTests(unittest.TestCase):
    """Tests for build_variant_metadata_rows: per-variant missing counts and order."""

    def test_standard_fixture_variant_order_and_keys(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_variant_metadata_rows("cohort", vcf)

        self.assertEqual(
            [row[2] for row in rows],
            ["chrPanel:100:A:G", "chrPanel:200:C:T", "chrPanel:300:G:A"],
        )
        self.assertEqual([row[1] for row in rows], ["0", "1", "2"])

    def test_missing_fixture_missing_counts_per_variant(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_variant_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0][8], "1")  # chrPanel:100 -- sample_a missing
        self.assertEqual(rows[0][9], "0.500000")


class GenotypeAccountingTests(unittest.TestCase):
    """Tests for build_genotype_accounting_rows: reason-specific totals."""

    def test_standard_fixture_accounting(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        metrics = {
            row[1]: row[2] for row in build_module.build_genotype_accounting_rows("cohort", vcf)
        }

        self.assertEqual(metrics["total_genotype_cells"], "6")
        self.assertEqual(metrics["standard_hom_ref_calls"], "2")
        self.assertEqual(metrics["standard_het_calls"], "3")
        self.assertEqual(metrics["standard_hom_alt_calls"], "1")
        self.assertEqual(metrics["total_treated_as_missing"], "0")

    def test_non_standard_fixture_accounting_breaks_down_every_reason(self) -> None:
        # var1: sample_a=0|1 (phased), sample_b=0/0 (standard)
        # var2: sample_a=0, sample_b=1 (both haploid -> non_diploid)
        # var3: sample_a=0/0/1 (triploid -> non_diploid), sample_b=0/1 (standard)
        # var4: sample_a=0/2 (non_biallelic_index), sample_b=0/0 (standard)
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        metrics = {
            row[1]: row[2] for row in build_module.build_genotype_accounting_rows("cohort", vcf)
        }

        self.assertEqual(metrics["phased_calls_treated_as_missing"], "1")
        self.assertEqual(metrics["non_diploid_calls_treated_as_missing"], "3")
        self.assertEqual(metrics["non_biallelic_index_calls_treated_as_missing"], "1")
        self.assertEqual(metrics["total_treated_as_missing"], "5")
        self.assertEqual(metrics["standard_hom_ref_calls"], "2")
        self.assertEqual(metrics["standard_het_calls"], "1")
        self.assertEqual(metrics["standard_hom_alt_calls"], "0")

    def test_summary_text_mentions_every_reason(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        summary = build_module.build_genotype_accounting_summary_text("cohort", vcf)

        self.assertIn("phased, treated as missing", summary)
        self.assertIn("non-diploid, treated as missing", summary)
        self.assertIn("non-biallelic-index, treated as missing", summary)


class OutputContractTests(unittest.TestCase):
    """Tests for output headers."""

    def test_headers(self) -> None:
        self.assertEqual(
            build_module.GENOTYPE_ACCOUNTING_HEADER, ("cohort_id", "metric", "value")
        )
        self.assertEqual(
            build_module.SAMPLE_METADATA_HEADER,
            (
                "cohort_id",
                "sample_index",
                "sample_id",
                "missing_genotype_count",
                "missing_genotype_rate",
                "non_standard_genotype_count",
            ),
        )
        self.assertEqual(
            build_module.VARIANT_METADATA_HEADER,
            (
                "cohort_id",
                "variant_index",
                "variant_key",
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "missing_genotype_count",
                "missing_genotype_rate",
            ),
        )


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_all_five_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "matrix.tsv.gz"
            sample_metadata_path = tmp_path / "sample_metadata.tsv"
            variant_metadata_path = tmp_path / "variant_metadata.tsv"
            accounting_path = tmp_path / "accounting.tsv"
            summary_path = tmp_path / "summary.txt"

            exit_code = build_module.main(
                [
                    "--gs-pass-vcf",
                    str(FIXTURES_DIR / "build_panel_standard.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--matrix-output",
                    str(matrix_path),
                    "--sample-metadata-output",
                    str(sample_metadata_path),
                    "--variant-metadata-output",
                    str(variant_metadata_path),
                    "--genotype-accounting-output",
                    str(accounting_path),
                    "--genotype-accounting-summary-output",
                    str(summary_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            for path in (
                matrix_path,
                sample_metadata_path,
                variant_metadata_path,
                accounting_path,
                summary_path,
            ):
                self.assertTrue(path.exists(), f"{path} was not written")

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_vcf = tmp_path / "does-not-exist.vcf.gz"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = build_module.main(
                    [
                        "--gs-pass-vcf",
                        str(missing_vcf),
                        "--cohort-id",
                        "cohort",
                        "--matrix-output",
                        str(tmp_path / "m.tsv.gz"),
                        "--sample-metadata-output",
                        str(tmp_path / "sm.tsv"),
                        "--variant-metadata-output",
                        str(tmp_path / "vm.tsv"),
                        "--genotype-accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--genotype-accounting-summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_vcf), stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "matrix.tsv.gz"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--gs-pass-vcf",
                    str(FIXTURES_DIR / "build_panel_empty.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--matrix-output",
                    str(matrix_path),
                    "--sample-metadata-output",
                    str(tmp_path / "sm.tsv"),
                    "--variant-metadata-output",
                    str(tmp_path / "vm.tsv"),
                    "--genotype-accounting-output",
                    str(tmp_path / "a.tsv"),
                    "--genotype-accounting-summary-output",
                    str(tmp_path / "s.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "variant_key\tsample_a\tsample_b\n")


if __name__ == "__main__":
    unittest.main()
