"""Unit tests for bin/summarize_filter_qc.py.

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
SCRIPT_PATH = REPO_ROOT / "bin" / "summarize_filter_qc.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses inspects sys.modules[cls.__module__] while the class
    # body executes, so the module must be registered before exec_module
    # runs, not just returned afterward.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


qc = _load_module("summarize_filter_qc", SCRIPT_PATH)


class ParseFilteredVcfTests(unittest.TestCase):
    """Tests for parse_filtered_vcf against real and hand-crafted VCFs."""

    def test_real_filtered_snp_vcf(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp.vcf.gz")

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].chrom, "chrSynthetic1")
        self.assertEqual(records[0].pos, "1501")
        self.assertEqual(records[0].filter_value, "SNP_SOR_HIGH")
        self.assertEqual(records[0].info["SOR"], "4.407")
        self.assertNotIn("MQRankSum", records[0].info)

    def test_real_empty_indel_vcf(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_indel_empty.vcf.gz")

        self.assertEqual(records, [])

    def test_multi_tag_fixture_parses_pass_and_multi_tag_rows(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")

        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].filter_value, "PASS")
        self.assertEqual(records[1].filter_value, "SNP_QD_LOW")
        self.assertEqual(records[2].filter_value, "SNP_QD_LOW;SNP_SOR_HIGH")
        self.assertEqual(records[1].info["MQRankSum"], "0.300")

    def test_malformed_row_raises_with_file_and_cause(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            broken_path = Path(tmp) / "broken.vcf.gz"
            _write_gzip_vcf(broken_path, ["#CHROM\tPOS", "chr1\t100\t.\tA\tT"])

            with self.assertRaises(qc.MalformedVcfError) as raised:
                qc.parse_filtered_vcf(broken_path)

            self.assertIn(str(broken_path), str(raised.exception))


class FilterBreakdownTests(unittest.TestCase):
    """Tests for compute_filter_breakdown: totals, PASS, tags, combinations."""

    def test_real_filtered_snp_all_non_pass_single_tag(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp.vcf.gz")
        breakdown = qc.compute_filter_breakdown(records)

        self.assertEqual(breakdown.total_records, 2)
        self.assertEqual(breakdown.pass_records, 0)
        self.assertEqual(breakdown.non_pass_records, 2)
        self.assertEqual(breakdown.multi_tag_records, 0)
        self.assertEqual(breakdown.combination_counts, {"SNP_SOR_HIGH": 2})
        self.assertEqual(breakdown.tag_counts, {"SNP_SOR_HIGH": 2})

    def test_empty_vcf_reports_all_zero(self) -> None:
        breakdown = qc.compute_filter_breakdown([])

        self.assertEqual(breakdown.total_records, 0)
        self.assertEqual(breakdown.pass_records, 0)
        self.assertEqual(breakdown.non_pass_records, 0)
        self.assertEqual(breakdown.multi_tag_records, 0)
        self.assertEqual(breakdown.combination_counts, {})
        self.assertEqual(breakdown.tag_counts, {})

    def test_multi_tag_fixture_reconciles_total_pass_nonpass(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")
        breakdown = qc.compute_filter_breakdown(records)

        self.assertEqual(breakdown.total_records, 3)
        self.assertEqual(breakdown.pass_records, 1)
        self.assertEqual(breakdown.non_pass_records, 2)
        self.assertEqual(breakdown.total_records, breakdown.pass_records + breakdown.non_pass_records)
        self.assertEqual(breakdown.multi_tag_records, 1)

    def test_multi_tag_fixture_tag_counts_can_exceed_non_pass_records(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")
        breakdown = qc.compute_filter_breakdown(records)

        # SNP_QD_LOW appears on two records, SNP_SOR_HIGH on one: 3 tag
        # occurrences total, spread across only 2 non-PASS records.
        self.assertEqual(breakdown.tag_counts, {"SNP_QD_LOW": 2, "SNP_SOR_HIGH": 1})
        self.assertGreater(sum(breakdown.tag_counts.values()), breakdown.non_pass_records)

    def test_multi_tag_fixture_combination_counts(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")
        breakdown = qc.compute_filter_breakdown(records)

        self.assertEqual(
            breakdown.combination_counts,
            {"PASS": 1, "SNP_QD_LOW": 1, "SNP_QD_LOW;SNP_SOR_HIGH": 1},
        )


class AnnotationCoverageTests(unittest.TestCase):
    """Tests distinguishing missing-annotation-unevaluated from filter-passed."""

    def test_real_filtered_snp_ranksum_missing_is_not_conflated_with_passed(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp.vcf.gz")
        coverages = {c.annotation: c for c in qc.compute_annotation_coverage(records, "snp")}

        mq_rank_sum = coverages["MQRankSum"]
        self.assertEqual(mq_rank_sum.present_records, 0)
        self.assertEqual(mq_rank_sum.missing_records, 2)
        self.assertEqual(mq_rank_sum.evaluable_rate, "0.000000")
        self.assertEqual(mq_rank_sum.filter_tag, "SNP_MQRANKSUM_LOW")
        self.assertEqual(mq_rank_sum.filter_tagged_records, "0")
        # NA, not "0.000000": 0 tagged out of 0 present is undefined, not
        # "evaluated and none failed".
        self.assertEqual(mq_rank_sum.filter_hit_rate, "NA")

        read_pos_rank_sum = coverages["ReadPosRankSum"]
        self.assertEqual(read_pos_rank_sum.evaluable_rate, "0.000000")
        self.assertEqual(read_pos_rank_sum.filter_hit_rate, "NA")

    def test_real_filtered_snp_sor_present_and_hit_for_both_records(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp.vcf.gz")
        coverages = {c.annotation: c for c in qc.compute_annotation_coverage(records, "snp")}

        sor = coverages["SOR"]
        self.assertEqual(sor.present_records, 2)
        self.assertEqual(sor.evaluable_rate, "1.000000")
        self.assertEqual(sor.filter_tagged_records, "2")
        self.assertEqual(sor.filter_hit_rate, "1.000000")

    def test_annotation_present_but_did_not_trigger_filter_is_zero_not_na(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")
        coverages = {c.annotation: c for c in qc.compute_annotation_coverage(records, "snp")}

        mq_rank_sum = coverages["MQRankSum"]
        # Present on exactly one record (pos 200) and evaluable there,
        # but that record's FILTER never got SNP_MQRANKSUM_LOW: a real
        # 0-out-of-1 hit rate, correctly distinct from "not evaluable".
        self.assertEqual(mq_rank_sum.present_records, 1)
        self.assertEqual(mq_rank_sum.filter_tagged_records, "0")
        self.assertEqual(mq_rank_sum.filter_hit_rate, "0.000000")

    def test_indel_has_no_filter_for_sor_mq_mqranksum(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_indel_empty.vcf.gz")
        coverages = {c.annotation: c for c in qc.compute_annotation_coverage(records, "indel")}

        for annotation in ("SOR", "MQ", "MQRankSum"):
            with self.subTest(annotation=annotation):
                self.assertEqual(coverages[annotation].filter_tag, "NA")
                self.assertEqual(coverages[annotation].filter_tagged_records, "NA")
                self.assertEqual(coverages[annotation].filter_hit_rate, "NA")

        # QD/QUAL/FS/ReadPosRankSum do apply to indels.
        self.assertEqual(coverages["QD"].filter_tag, "INDEL_QD_LOW")
        self.assertEqual(coverages["ReadPosRankSum"].filter_tag, "INDEL_READPOSRANKSUM_LOW")

    def test_empty_vcf_evaluable_rate_is_na_not_zero(self) -> None:
        coverages = qc.compute_annotation_coverage([], "snp")

        for coverage in coverages:
            with self.subTest(annotation=coverage.annotation):
                self.assertEqual(coverage.total_records, 0)
                self.assertEqual(coverage.evaluable_rate, "NA")


class OutputContractTests(unittest.TestCase):
    """Tests for output headers and row structure."""

    def test_filter_breakdown_header(self) -> None:
        self.assertEqual(
            qc.FILTER_BREAKDOWN_HEADER,
            ("cohort_id", "stage", "variant_type", "category", "key", "record_count"),
        )

    def test_annotation_qc_header(self) -> None:
        self.assertEqual(
            qc.ANNOTATION_QC_HEADER,
            (
                "cohort_id",
                "stage",
                "variant_type",
                "annotation",
                "total_records",
                "present_records",
                "missing_records",
                "evaluable_rate",
                "filter_tag",
                "filter_tagged_records",
                "filter_hit_rate",
            ),
        )

    def test_annotation_qc_rows_cover_all_seven_annotations_in_order(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp.vcf.gz")
        coverages = qc.compute_annotation_coverage(records, "snp")
        rows = qc.build_annotation_qc_rows("cohort", "filtered", "snp", coverages)

        self.assertEqual(
            [row[3] for row in rows],
            ["QD", "QUAL", "SOR", "FS", "MQ", "MQRankSum", "ReadPosRankSum"],
        )

    def test_summary_text_states_total_equals_pass_plus_nonpass(self) -> None:
        records = qc.parse_filtered_vcf(FIXTURES_DIR / "filtered_snp_multi_tag.vcf.gz")
        breakdown = qc.compute_filter_breakdown(records)
        coverages = qc.compute_annotation_coverage(records, "snp")

        summary_text = qc.build_filter_qc_summary_text("cohort", "filtered", "snp", breakdown, coverages)

        self.assertIn("Reconciliation: total (3) = PASS (1) + non-PASS (2)", summary_text)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_all_three_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            filter_breakdown_output = tmp_path / "filter_breakdown.tsv"
            annotation_qc_output = tmp_path / "annotation_qc.tsv"
            summary_output = tmp_path / "summary.txt"

            exit_code = qc.main(
                [
                    "--filtered-vcf",
                    str(FIXTURES_DIR / "filtered_snp.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--stage",
                    "filtered",
                    "--variant-type",
                    "snp",
                    "--filter-breakdown-output",
                    str(filter_breakdown_output),
                    "--annotation-qc-output",
                    str(annotation_qc_output),
                    "--summary-output",
                    str(summary_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(filter_breakdown_output.exists())
            self.assertTrue(annotation_qc_output.exists())
            self.assertTrue(summary_output.exists())

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_input = tmp_path / "does-not-exist.vcf.gz"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = qc.main(
                    [
                        "--filtered-vcf",
                        str(missing_input),
                        "--cohort-id",
                        "cohort",
                        "--stage",
                        "filtered",
                        "--variant-type",
                        "snp",
                        "--filter-breakdown-output",
                        str(tmp_path / "b.tsv"),
                        "--annotation-qc-output",
                        str(tmp_path / "a.tsv"),
                        "--summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_input), stderr.getvalue())

    def test_main_rejects_unknown_variant_type(self) -> None:
        stderr = io.StringIO()

        with self.assertRaises(SystemExit) as raised, contextlib.redirect_stderr(stderr):
            qc.main(
                [
                    "--filtered-vcf",
                    str(FIXTURES_DIR / "filtered_snp.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--stage",
                    "filtered",
                    "--variant-type",
                    "mnp",
                    "--filter-breakdown-output",
                    "b.tsv",
                    "--annotation-qc-output",
                    "a.tsv",
                    "--summary-output",
                    "s.txt",
                ]
            )

        self.assertEqual(raised.exception.code, 2)

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            filter_breakdown_output = tmp_path / "filter_breakdown.tsv"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--filtered-vcf",
                    str(FIXTURES_DIR / "filtered_indel_empty.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--stage",
                    "filtered",
                    "--variant-type",
                    "indel",
                    "--filter-breakdown-output",
                    str(filter_breakdown_output),
                    "--annotation-qc-output",
                    str(tmp_path / "a.tsv"),
                    "--summary-output",
                    str(tmp_path / "s.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(filter_breakdown_output.exists())


def _write_gzip_vcf(path: Path, lines: list[str]) -> None:
    """Write minimal gzip-compressed VCF-like text for a test fixture."""
    import gzip

    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    unittest.main()
