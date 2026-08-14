"""Unit tests for bin/reconcile_gs_panel_accounting.py.

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


def _vcf(
    record_count: int,
    sample_count: int,
    *,
    sample_ids: tuple[str, ...] | None = None,
    variant_keys: tuple[str, ...] | None = None,
    all_rows_ok: bool = True,
):
    return reconcile_module.VcfSummary(
        record_count=record_count,
        sample_count=sample_count,
        sample_ids=sample_ids if sample_ids is not None else tuple(
            f"sample_{i}" for i in range(sample_count)
        ),
        variant_keys=variant_keys if variant_keys is not None else tuple(
            f"chrTest:{100 * (i + 1)}:A:G" for i in range(record_count)
        ),
        all_rows_have_expected_sample_count=all_rows_ok,
    )


def _matrix(
    variant_count: int,
    sample_count: int,
    *,
    sample_ids: tuple[str, ...] | None = None,
    variant_keys: tuple[str, ...] | None = None,
    all_rows_ok: bool = True,
):
    return reconcile_module.MatrixSummary(
        variant_count=variant_count,
        sample_count=sample_count,
        sample_ids=sample_ids if sample_ids is not None else tuple(
            f"sample_{i}" for i in range(sample_count)
        ),
        variant_keys=variant_keys if variant_keys is not None else tuple(
            f"chrTest:{100 * (i + 1)}:A:G" for i in range(variant_count)
        ),
        all_rows_have_expected_width=all_rows_ok,
    )


class SummarizeVcfTests(unittest.TestCase):
    """Tests for summarize_vcf: record/sample counts, identity, and row-shape checks."""

    def test_counts_records_and_samples(self) -> None:
        summary = reconcile_module.summarize_vcf(FIXTURES_DIR / "reconcile_gs_clean_raw_all.vcf.gz")
        self.assertEqual(summary.record_count, 3)
        self.assertEqual(summary.sample_count, 2)
        self.assertTrue(summary.all_rows_have_expected_sample_count)

    def test_reads_sample_ids_in_header_order(self) -> None:
        summary = reconcile_module.summarize_vcf(FIXTURES_DIR / "reconcile_gs_clean_raw_all.vcf.gz")
        self.assertEqual(summary.sample_ids, ("sample_a", "sample_b"))

    def test_reads_variant_keys_in_file_order(self) -> None:
        summary = reconcile_module.summarize_vcf(FIXTURES_DIR / "reconcile_gs_clean_raw_all.vcf.gz")
        self.assertEqual(
            summary.variant_keys,
            ("chrTest:100:A:G", "chrTest:200:C:T", "chrTest:300:A:ATT"),
        )

    def test_empty_vcf_has_zero_records_but_keeps_sample_count(self) -> None:
        summary = reconcile_module.summarize_vcf(FIXTURES_DIR / "reconcile_gs_empty_pass.vcf.gz")
        self.assertEqual(summary.record_count, 0)
        self.assertEqual(summary.sample_count, 2)
        self.assertEqual(summary.variant_keys, ())
        self.assertTrue(summary.all_rows_have_expected_sample_count)

    def test_row_with_wrong_sample_field_count_is_flagged_not_rejected(self) -> None:
        summary = reconcile_module.summarize_vcf(
            FIXTURES_DIR / "reconcile_gs_malformed_row_raw_all.vcf.gz"
        )
        self.assertEqual(summary.record_count, 2)
        self.assertFalse(summary.all_rows_have_expected_sample_count)
        # the variant key is still captured even though the sample-field
        # count for this row is wrong -- CHROM/POS/REF/ALT are unaffected.
        self.assertEqual(summary.variant_keys, ("chrTest:100:A:G", "chrTest:200:C:T"))

    def test_missing_chrom_header_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "no_header.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")

            with self.assertRaises(reconcile_module.MalformedVcfError) as raised:
                reconcile_module.summarize_vcf(path)
            self.assertIn("no #CHROM header", str(raised.exception))

    def test_data_row_before_header_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "data_before_header.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("chrTest\t100\t.\tA\tG\t100\t.\t.\tGT\t0/1\n")

            with self.assertRaises(reconcile_module.MalformedVcfError) as raised:
                reconcile_module.summarize_vcf(path)
            self.assertIn("before #CHROM header", str(raised.exception))

    def test_chrom_header_with_no_samples_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "no_samples.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            with self.assertRaises(reconcile_module.MalformedVcfError) as raised:
                reconcile_module.summarize_vcf(path)
            self.assertIn("expected at least 10", str(raised.exception))

    def test_data_row_with_too_few_fields_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "short_row.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\n")
                handle.write("chrTest\t100\t.\tA\tG\n")

            with self.assertRaises(reconcile_module.MalformedVcfError) as raised:
                reconcile_module.summarize_vcf(path)
            self.assertIn("expected at least 10", str(raised.exception))


class SummarizeMatrixTests(unittest.TestCase):
    """Tests for summarize_matrix: header validation, identity, and row-width checks."""

    def test_counts_variant_and_sample_columns(self) -> None:
        summary = reconcile_module.summarize_matrix(FIXTURES_DIR / "reconcile_gs_clean_matrix.tsv.gz")
        self.assertEqual(summary.variant_count, 1)
        self.assertEqual(summary.sample_count, 2)
        self.assertTrue(summary.all_rows_have_expected_width)

    def test_reads_sample_ids_and_variant_keys(self) -> None:
        summary = reconcile_module.summarize_matrix(FIXTURES_DIR / "reconcile_gs_clean_matrix.tsv.gz")
        self.assertEqual(summary.sample_ids, ("sample_a", "sample_b"))
        self.assertEqual(summary.variant_keys, ("chrTest:100:A:G",))

    def test_header_only_matrix_has_zero_variants_but_keeps_samples(self) -> None:
        summary = reconcile_module.summarize_matrix(FIXTURES_DIR / "reconcile_gs_empty_matrix.tsv.gz")
        self.assertEqual(summary.variant_count, 0)
        self.assertEqual(summary.sample_count, 2)
        self.assertEqual(summary.variant_keys, ())

    def test_row_with_dropped_column_is_flagged_not_rejected(self) -> None:
        summary = reconcile_module.summarize_matrix(
            FIXTURES_DIR / "reconcile_gs_dropped_column_matrix.tsv.gz"
        )
        # the row is still counted and its variant_key still captured --
        # only the width flag records that something is wrong.
        self.assertEqual(summary.variant_count, 1)
        self.assertEqual(summary.variant_keys, ("chrTest:100:A:G",))
        self.assertFalse(summary.all_rows_have_expected_width)

    def test_empty_file_raises(self) -> None:
        with self.assertRaises(reconcile_module.MalformedMatrixError) as raised:
            reconcile_module.summarize_matrix(FIXTURES_DIR / "reconcile_gs_empty_file_matrix.tsv.gz")
        self.assertIn("file is empty", str(raised.exception))

    def test_wrong_first_header_column_raises(self) -> None:
        with self.assertRaises(reconcile_module.MalformedMatrixError) as raised:
            reconcile_module.summarize_matrix(
                FIXTURES_DIR / "reconcile_gs_malformed_header_matrix.tsv.gz"
            )
        self.assertIn("variant_key", str(raised.exception))


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


class ReadMetadataColumnTests(unittest.TestCase):
    """Tests for read_metadata_column: named-column extraction, in row order."""

    def test_reads_sample_id_column_in_order(self) -> None:
        values = reconcile_module.read_metadata_column(
            FIXTURES_DIR / "reconcile_gs_clean_sample_metadata.tsv", "sample_id"
        )
        self.assertEqual(values, ("sample_a", "sample_b"))

    def test_reads_variant_key_column_in_order(self) -> None:
        values = reconcile_module.read_metadata_column(
            FIXTURES_DIR / "reconcile_gs_clean_variant_metadata.tsv", "variant_key"
        )
        self.assertEqual(values, ("chrTest:100:A:G",))

    def test_zero_data_rows_is_valid(self) -> None:
        values = reconcile_module.read_metadata_column(
            FIXTURES_DIR / "reconcile_gs_empty_variant_metadata.tsv", "variant_key"
        )
        self.assertEqual(values, ())

    def test_missing_column_raises(self) -> None:
        with self.assertRaises(reconcile_module.MalformedAccountingError) as raised:
            reconcile_module.read_metadata_column(
                FIXTURES_DIR / "reconcile_gs_clean_sample_metadata.tsv", "no_such_column"
            )
        self.assertIn("no_such_column", str(raised.exception))

    def test_row_too_short_for_column_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "short_row.tsv"
            path.write_text("a\tb\tc\nx\ty\n", encoding="utf-8")

            with self.assertRaises(reconcile_module.MalformedAccountingError) as raised:
                reconcile_module.read_metadata_column(path, "c")
            self.assertIn(str(path), str(raised.exception))


class DescribeSequenceMismatchTests(unittest.TestCase):
    """Tests for _describe_sequence_mismatch: bounded, useful diagnostics."""

    def test_reports_length_difference(self) -> None:
        message = reconcile_module._describe_sequence_mismatch(
            [("a", ("x", "y")), ("b", ("x",))]
        )
        self.assertIn("a=2", message)
        self.assertIn("b=1", message)

    def test_reports_first_differing_index(self) -> None:
        message = reconcile_module._describe_sequence_mismatch(
            [("a", ("x", "y", "z")), ("b", ("x", "Y", "z"))]
        )
        self.assertIn("index 1", message)
        self.assertIn("a[1]='y'", message)
        self.assertIn("b[1]='Y'", message)


class ReconcileTests(unittest.TestCase):
    """Tests for reconcile: cross-file consistency checks and lineage arithmetic."""

    def test_clean_populated_scenario(self) -> None:
        result = reconcile_module.reconcile(
            raw_all=_vcf(3, 2),
            normalized=_vcf(3, 2),
            classified_biallelic_snp_records=2,
            gs_pass=_vcf(1, 2, variant_keys=("chrTest:100:A:G",)),
            matrix=_matrix(1, 2, variant_keys=("chrTest:100:A:G",)),
            variant_metadata_keys=("chrTest:100:A:G",),
            sample_metadata_ids=("sample_0", "sample_1"),
        )

        self.assertEqual(result.raw_all_records, 3)
        self.assertEqual(result.gs_hard_filter_excluded_records, 1)
        self.assertEqual(result.gs_pass_records, 1)
        self.assertEqual(result.matrix_variant_records, 1)
        self.assertEqual(result.matrix_sample_count, 2)
        self.assertEqual(result.panel_status, "populated")

    def test_multi_variant_matching_order_is_not_flagged(self) -> None:
        """A larger, genuinely-agreeing scenario: order preserved everywhere."""
        keys = ("chrTest:100:A:G", "chrTest:200:C:T", "chrTest:300:A:G")
        samples = ("sample_a", "sample_b", "sample_c")

        result = reconcile_module.reconcile(
            raw_all=_vcf(3, 3, sample_ids=samples),
            normalized=_vcf(3, 3, sample_ids=samples),
            classified_biallelic_snp_records=3,
            gs_pass=_vcf(3, 3, sample_ids=samples, variant_keys=keys),
            matrix=_matrix(3, 3, sample_ids=samples, variant_keys=keys),
            variant_metadata_keys=keys,
            sample_metadata_ids=samples,
        )

        self.assertEqual(result.matrix_variant_records, 3)
        self.assertEqual(result.panel_status, "populated")

    def test_empty_panel_is_not_flagged_as_an_error(self) -> None:
        result = reconcile_module.reconcile(
            raw_all=_vcf(2, 2),
            normalized=_vcf(2, 2),
            classified_biallelic_snp_records=2,
            gs_pass=_vcf(0, 2, variant_keys=()),
            matrix=_matrix(0, 2, variant_keys=()),
            variant_metadata_keys=(),
            sample_metadata_ids=("sample_0", "sample_1"),
        )

        self.assertEqual(result.panel_status, "empty")
        self.assertEqual(result.matrix_sample_count, 2)

    def test_negative_excluded_records_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=1,
                gs_pass=_vcf(2, 2),
                matrix=_matrix(2, 2),
                variant_metadata_keys=_vcf(2, 2).variant_keys,
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("negative", str(raised.exception))

    def test_matrix_variant_count_mismatch_raises(self) -> None:
        """The core matrix-corruption regression: a matrix that silently
        lost its one data row must be caught even though nothing else
        looks wrong."""
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2, variant_keys=("chrTest:100:A:G",)),
                matrix=_matrix(0, 2, variant_keys=()),
                variant_metadata_keys=("chrTest:100:A:G",),
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("variant identity and/or order disagree", str(raised.exception))

    def test_matrix_row_width_mismatch_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2, variant_keys=("chrTest:100:A:G",)),
                matrix=_matrix(1, 2, variant_keys=("chrTest:100:A:G",), all_rows_ok=False),
                variant_metadata_keys=("chrTest:100:A:G",),
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("column count", str(raised.exception))

    def test_variant_key_value_mismatch_raises_even_with_equal_counts(self) -> None:
        """Same variant *count* everywhere, but the matrix's own
        variant_key names a different variant than the PASS VCF."""
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2, variant_keys=("chrTest:100:A:G",)),
                matrix=_matrix(1, 2, variant_keys=("chrTest:999:A:G",)),
                variant_metadata_keys=("chrTest:100:A:G",),
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("variant identity and/or order disagree", str(raised.exception))

    def test_variant_order_mismatch_raises_even_with_equal_sets(self) -> None:
        """Same two variant keys everywhere, but variant_metadata lists
        them in the opposite order from the matrix and the PASS VCF."""
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(2, 2),
                normalized=_vcf(2, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(2, 2, variant_keys=("chrTest:100:A:G", "chrTest:200:C:T")),
                matrix=_matrix(2, 2, variant_keys=("chrTest:100:A:G", "chrTest:200:C:T")),
                variant_metadata_keys=("chrTest:200:C:T", "chrTest:100:A:G"),
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("same length (2), first disagreement at index 0", str(raised.exception))

    def test_sample_order_mismatch_raises_even_with_equal_sets(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(1, 2),
                normalized=_vcf(1, 2),
                classified_biallelic_snp_records=1,
                gs_pass=_vcf(1, 2, sample_ids=("sample_a", "sample_b")),
                matrix=_matrix(1, 2, sample_ids=("sample_b", "sample_a")),
                variant_metadata_keys=_vcf(1, 2).variant_keys,
                sample_metadata_ids=("sample_a", "sample_b"),
            )
        self.assertIn("sample identity and/or order disagree", str(raised.exception))

    def test_sample_metadata_mismatch_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(1, 2),
                normalized=_vcf(1, 2),
                classified_biallelic_snp_records=1,
                gs_pass=_vcf(1, 2),
                matrix=_matrix(1, 2),
                variant_metadata_keys=_vcf(1, 2).variant_keys,
                sample_metadata_ids=("sample_0",),
            )
        self.assertIn("sample identity and/or order disagree", str(raised.exception))

    def test_raw_all_row_shape_mismatch_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2, all_rows_ok=False),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2),
                matrix=_matrix(1, 2),
                variant_metadata_keys=_vcf(1, 2).variant_keys,
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("raw/all", str(raised.exception))

    def test_normalized_row_shape_mismatch_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2, all_rows_ok=False),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2),
                matrix=_matrix(1, 2),
                variant_metadata_keys=_vcf(1, 2).variant_keys,
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("normalized", str(raised.exception))

    def test_gs_pass_row_shape_mismatch_raises(self) -> None:
        with self.assertRaises(reconcile_module.InconsistentGsPanelError) as raised:
            reconcile_module.reconcile(
                raw_all=_vcf(3, 2),
                normalized=_vcf(3, 2),
                classified_biallelic_snp_records=2,
                gs_pass=_vcf(1, 2, all_rows_ok=False),
                matrix=_matrix(1, 2),
                variant_metadata_keys=_vcf(1, 2).variant_keys,
                sample_metadata_ids=("sample_0", "sample_1"),
            )
        self.assertIn("PASS VCF", str(raised.exception))


class OutputContractTests(unittest.TestCase):
    """Tests for output header, rows, and summary text."""

    def _clean_result(self):
        return reconcile_module.reconcile(
            raw_all=_vcf(3, 2),
            normalized=_vcf(3, 2),
            classified_biallelic_snp_records=2,
            gs_pass=_vcf(1, 2, variant_keys=("chrTest:100:A:G",)),
            matrix=_matrix(1, 2, variant_keys=("chrTest:100:A:G",)),
            variant_metadata_keys=("chrTest:100:A:G",),
            sample_metadata_ids=("sample_0", "sample_1"),
        )

    def test_output_header(self) -> None:
        self.assertEqual(reconcile_module.OUTPUT_HEADER, ("cohort_id", "metric", "value"))

    def test_build_output_rows_includes_every_metric(self) -> None:
        rows = reconcile_module.build_output_rows("cohort", self._clean_result())
        values = {row[1]: row[2] for row in rows}
        self.assertEqual(values["raw_all_records"], "3")
        self.assertEqual(values["matrix_variant_records"], "1")
        self.assertEqual(values["matrix_sample_count"], "2")
        self.assertEqual(values["panel_status"], "populated")

    def test_summary_states_cross_check_passed(self) -> None:
        summary = reconcile_module.build_summary_text("cohort", self._clean_result())
        self.assertIn("Cross-checked directly against the panel's own artifacts", summary)

    def test_summary_explains_empty_panel_is_normal(self) -> None:
        result = reconcile_module.reconcile(
            raw_all=_vcf(2, 2),
            normalized=_vcf(2, 2),
            classified_biallelic_snp_records=2,
            gs_pass=_vcf(0, 2, variant_keys=()),
            matrix=_matrix(0, 2, variant_keys=()),
            variant_metadata_keys=(),
            sample_metadata_ids=("sample_0", "sample_1"),
        )
        summary = reconcile_module.build_summary_text("cohort", result)
        self.assertIn("normal outcome, not an error", summary)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and hard-failure exit codes."""

    def _run_main(
        self,
        prefix: str,
        *,
        matrix_override: Path | None = None,
        variant_metadata_override: Path | None = None,
    ) -> tuple[int, Path, Path, str]:
        tmp_path = Path(tempfile.mkdtemp())
        output_path = tmp_path / "accounting.tsv"
        summary_path = tmp_path / "summary.txt"
        stderr = io.StringIO()

        matrix = matrix_override or (FIXTURES_DIR / f"reconcile_gs_{prefix}_matrix.tsv.gz")
        variant_metadata = variant_metadata_override or (
            FIXTURES_DIR / f"reconcile_gs_{prefix}_variant_metadata.tsv"
        )

        with contextlib.redirect_stderr(stderr):
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
                    "--matrix",
                    str(matrix),
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

        return exit_code, output_path, summary_path, stderr.getvalue()

    def test_main_succeeds_for_clean_scenario(self) -> None:
        exit_code, output_path, summary_path, _stderr = self._run_main("clean")

        self.assertEqual(exit_code, 0)
        values = {
            line.split("\t")[1]: line.split("\t")[2]
            for line in output_path.read_text(encoding="utf-8").splitlines()[1:]
        }
        self.assertEqual(values["raw_all_records"], "3")
        self.assertEqual(values["classified_biallelic_snp_records"], "2")
        self.assertEqual(values["gs_pass_records"], "1")
        self.assertEqual(values["matrix_variant_records"], "1")
        self.assertEqual(values["matrix_sample_count"], "2")
        self.assertEqual(values["panel_status"], "populated")
        self.assertTrue(summary_path.exists())

    def test_main_reports_empty_panel_status(self) -> None:
        exit_code, output_path, _summary_path, _stderr = self._run_main("empty")

        self.assertEqual(exit_code, 0)
        values = {
            line.split("\t")[1]: line.split("\t")[2]
            for line in output_path.read_text(encoding="utf-8").splitlines()[1:]
        }
        self.assertEqual(values["gs_pass_records"], "0")
        self.assertEqual(values["matrix_variant_records"], "0")
        self.assertEqual(values["matrix_sample_count"], "2")
        self.assertEqual(values["panel_status"], "empty")

    def test_main_fails_when_matrix_disagrees_with_metadata_that_still_looks_correct(self) -> None:
        """The core matrix-corruption regression: a matrix that silently
        lost its one data row must be caught even though variant/sample
        metadata still report the original, correct-looking counts."""
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean", matrix_override=FIXTURES_DIR / "reconcile_gs_matrix_corrupt_matrix.tsv.gz"
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("variant identity and/or order disagree", stderr)

    def test_main_fails_when_matrix_drops_a_sample_column(self) -> None:
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean", matrix_override=FIXTURES_DIR / "reconcile_gs_sample_mismatch_matrix.tsv.gz"
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("sample identity and/or order disagree", stderr)

    def test_main_fails_when_matrix_row_has_a_dropped_dosage_column(self) -> None:
        """Same row/sample counts as the clean fixture; only the row's
        own column width is wrong."""
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean", matrix_override=FIXTURES_DIR / "reconcile_gs_dropped_column_matrix.tsv.gz"
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("column count", stderr)

    def test_main_fails_when_matrix_sample_order_is_swapped(self) -> None:
        """Same sample set and count as the clean fixture; only the
        matrix header's column order is wrong."""
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean", matrix_override=FIXTURES_DIR / "reconcile_gs_sample_reordered_matrix.tsv.gz"
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("sample identity and/or order disagree", stderr)

    def test_main_fails_when_matrix_variant_key_is_swapped(self) -> None:
        """Same row/sample counts as the clean fixture; the matrix's
        variant_key names a different variant than the PASS VCF."""
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean", matrix_override=FIXTURES_DIR / "reconcile_gs_variant_key_swapped_matrix.tsv.gz"
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("variant identity and/or order disagree", stderr)

    def test_main_fails_when_variant_metadata_is_reordered(self) -> None:
        """Same two variant keys everywhere, but variant_metadata lists
        them in the opposite order from the matrix and the PASS VCF --
        counts alone cannot catch this."""
        exit_code, output_path, _summary_path, stderr = self._run_main("reordered")

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("variant identity and/or order disagree", stderr)

    def test_main_fails_when_variant_metadata_disagrees(self) -> None:
        exit_code, output_path, _summary_path, stderr = self._run_main(
            "clean",
            variant_metadata_override=FIXTURES_DIR / "reconcile_gs_mismatch_variant_metadata.tsv",
        )

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("variant identity and/or order disagree", stderr)

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
                    "--matrix",
                    str(FIXTURES_DIR / "reconcile_gs_clean_matrix.tsv.gz"),
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
                "--matrix",
                str(FIXTURES_DIR / "reconcile_gs_clean_matrix.tsv.gz"),
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
