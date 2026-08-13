"""Unit tests for bin/summarize_variant_qc.py.

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
SCRIPT_PATH = REPO_ROOT / "bin" / "summarize_variant_qc.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_summarizer_module() -> types.ModuleType:
    """Load bin/summarize_variant_qc.py by path, without needing a package."""
    spec = importlib.util.spec_from_file_location("summarize_variant_qc", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses inspects sys.modules[cls.__module__] while the class
    # body executes, so the module must be registered before exec_module
    # runs, not just returned afterward.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


qc = _load_summarizer_module()


def _bcftools_stats_text(
    *,
    number_of_records: int = 2,
    samples: list[str] | None = None,
    missing_by_sample: dict[str, int] | None = None,
    include_tstv: bool = True,
    tstv_ratio: str = "0.00",
    include_psc: bool = True,
    omit_summary_keys: list[str] | None = None,
) -> str:
    """Build a minimal synthetic `bcftools stats`-shaped report for tests.

    Only the SN/TSTV/PSC rows the summarizer reads are included; real
    `bcftools stats` output also contains comment lines and many other
    record types (ID, SiS, AF, QUAL, ST, DP, PSI, HWE, VAF), which the
    "unknown record types" test adds explicitly.
    """
    samples = samples if samples is not None else ["sample_a", "sample_b"]
    missing_by_sample = missing_by_sample or {}
    omit_summary_keys = omit_summary_keys or []

    lines = ["# synthetic bcftools stats fixture for unit tests"]

    summary_fields = {
        "number of samples:": str(len(samples)),
        "number of records:": str(number_of_records),
        "number of SNPs:": str(number_of_records),
        "number of MNPs:": "0",
        "number of indels:": "0",
        "number of others:": "0",
        "number of multiallelic sites:": "0",
    }
    for key, value in summary_fields.items():
        if key not in omit_summary_keys:
            lines.append(f"SN\t0\t{key}\t{value}")

    if include_tstv:
        lines.append(
            f"TSTV\t0\t0\t{number_of_records}\t{tstv_ratio}\t0\t{number_of_records}\t{tstv_ratio}"
        )

    if include_psc:
        for sample in samples:
            missing = missing_by_sample.get(sample, 0)
            lines.append(f"PSC\t0\t{sample}\t1\t1\t0\t0\t1\t0\t8.0\t1\t0\t0\t{missing}")

    return "\n".join(lines) + "\n"


class ParseBcftoolsStatsTests(unittest.TestCase):
    """Tests for parse_bcftools_stats against real and synthetic reports."""

    def test_normal_snp_stats_from_real_bcftools_output(self) -> None:
        source = FIXTURES_DIR / "raw_snp.bcftools.stats.tsv"
        parsed = qc.parse_bcftools_stats(source.read_text(encoding="utf-8"), source)

        self.assertEqual(parsed.summary.number_of_samples, 2)
        self.assertEqual(parsed.summary.number_of_records, 2)
        self.assertEqual(parsed.summary.number_of_snps, 2)
        self.assertEqual(parsed.transition_transversion.transitions, 0)
        self.assertEqual(parsed.transition_transversion.transversions, 2)
        self.assertEqual(parsed.transition_transversion.ratio, "0.00")
        self.assertEqual([sample.sample for sample in parsed.per_sample], ["sample_a", "sample_b"])
        self.assertEqual(parsed.per_sample[0].average_depth, "8.0")
        self.assertEqual(parsed.per_sample[0].missing, 0)
        self.assertEqual(parsed.per_sample[1].average_depth, "5.5")

    def test_empty_indel_stats_has_zero_records_but_present_tstv_and_psc(self) -> None:
        # Real bcftools stats output for a 0-record VCF still emits TSTV
        # and PSC rows (all zero), rather than omitting them.
        source = FIXTURES_DIR / "empty_indel.bcftools.stats.tsv"
        parsed = qc.parse_bcftools_stats(source.read_text(encoding="utf-8"), source)

        self.assertEqual(parsed.summary.number_of_records, 0)
        self.assertEqual(parsed.transition_transversion.ratio, "0.00")
        self.assertEqual(len(parsed.per_sample), 2)
        self.assertTrue(all(sample.missing == 0 for sample in parsed.per_sample))

    def test_multiple_samples_preserve_encounter_order(self) -> None:
        text = _bcftools_stats_text(samples=["sample_b", "sample_a"])
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertEqual([sample.sample for sample in parsed.per_sample], ["sample_b", "sample_a"])

    def test_missing_genotypes_are_tracked_per_sample(self) -> None:
        text = _bcftools_stats_text(number_of_records=3, missing_by_sample={"sample_a": 1, "sample_b": 3})
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        missing_by_sample = {sample.sample: sample.missing for sample in parsed.per_sample}
        self.assertEqual(missing_by_sample, {"sample_a": 1, "sample_b": 3})

    def test_tstv_ratio_zero_is_kept_as_is(self) -> None:
        text = _bcftools_stats_text(tstv_ratio="0.00")
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertIsNotNone(parsed.transition_transversion)
        self.assertEqual(parsed.transition_transversion.ratio, "0.00")

    def test_missing_tstv_row_defaults_to_zero_counts(self) -> None:
        text = _bcftools_stats_text(include_tstv=False)
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertIsNone(parsed.transition_transversion)

    def test_missing_psc_rows_yields_no_samples(self) -> None:
        text = _bcftools_stats_text(include_psc=False)
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertEqual(parsed.per_sample, ())

    def test_unknown_record_types_and_comments_are_ignored(self) -> None:
        text = _bcftools_stats_text() + (
            "# a comment line with no tabs\n"
            "ID\t0\tcohort.vcf.gz\n"
            "SiS\t0\t1\t0\t0\t0\t0\t0\t0\t0\n"
            "UNKNOWN\t0\tsomething\tweird\n"
        )

        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertEqual(parsed.summary.number_of_records, 2)
        self.assertEqual(len(parsed.per_sample), 2)

    def test_missing_required_summary_metric_raises_with_file_and_cause(self) -> None:
        text = _bcftools_stats_text(omit_summary_keys=["number of samples:"])
        source = Path("/tmp/broken-summary.bcftools.stats.tsv")

        with self.assertRaises(qc.MalformedBcftoolsStatsError) as raised:
            qc.parse_bcftools_stats(text, source)

        message = str(raised.exception)
        self.assertIn(str(source), message)
        self.assertIn("number of samples:", message)

    def test_non_integer_summary_value_raises_with_file_and_cause(self) -> None:
        text = _bcftools_stats_text().replace(
            "SN\t0\tnumber of records:\t2",
            "SN\t0\tnumber of records:\tnot-a-number",
        )
        source = Path("/tmp/broken-value.bcftools.stats.tsv")

        with self.assertRaises(qc.MalformedBcftoolsStatsError) as raised:
            qc.parse_bcftools_stats(text, source)

        message = str(raised.exception)
        self.assertIn(str(source), message)
        self.assertIn("not-a-number", message)

    def test_malformed_psc_row_with_too_few_fields_raises(self) -> None:
        text = _bcftools_stats_text() + "PSC\t0\tsample_c\t1\n"
        source = Path("/tmp/broken-psc.bcftools.stats.tsv")

        with self.assertRaises(qc.MalformedBcftoolsStatsError) as raised:
            qc.parse_bcftools_stats(text, source)

        self.assertIn(str(source), str(raised.exception))


class BuildOutputTests(unittest.TestCase):
    """Tests for the output-contract builders: headers, metrics, and formatting."""

    def test_variant_qc_header_and_metric_order(self) -> None:
        self.assertEqual(
            qc.VARIANT_QC_HEADER,
            ("cohort_id", "stage", "variant_type", "metric", "value"),
        )

        source = FIXTURES_DIR / "raw_snp.bcftools.stats.tsv"
        parsed = qc.parse_bcftools_stats(source.read_text(encoding="utf-8"), source)
        rows = qc.build_variant_qc_rows("cohort", "raw", "snp", parsed)

        self.assertEqual(
            [row[3] for row in rows],
            [
                "number_of_samples",
                "number_of_records",
                "number_of_snps",
                "number_of_mnps",
                "number_of_indels",
                "number_of_others",
                "number_of_multiallelic_sites",
                "transitions",
                "transversions",
                "transition_transversion_ratio",
                "cohort_missing_genotypes",
                "cohort_total_genotypes",
                "cohort_missingness_rate",
                "sample_names",
            ],
        )

    def test_sample_qc_header(self) -> None:
        self.assertEqual(
            qc.SAMPLE_QC_HEADER,
            (
                "cohort_id",
                "stage",
                "variant_type",
                "sample",
                "reference_homozygous",
                "non_reference_homozygous",
                "heterozygous",
                "missing",
                "missingness_rate",
                "average_depth",
                "singletons",
            ),
        )

    def test_missingness_rate_uses_six_decimal_places(self) -> None:
        text = _bcftools_stats_text(number_of_records=3, missing_by_sample={"sample_a": 1})
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        rows = qc.build_sample_qc_rows("cohort", "raw", "snp", parsed)
        sample_a_row = next(row for row in rows if row[3] == "sample_a")

        self.assertEqual(sample_a_row[8], "0.333333")

    def test_empty_records_yields_na_missingness_not_zero(self) -> None:
        text = _bcftools_stats_text(number_of_records=0)
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        sample_rows = qc.build_sample_qc_rows("cohort", "raw", "snp", parsed)
        variant_rows = {row[3]: row[4] for row in qc.build_variant_qc_rows("cohort", "raw", "snp", parsed)}

        self.assertTrue(all(row[8] == "NA" for row in sample_rows))
        self.assertEqual(variant_rows["cohort_missingness_rate"], "NA")

    def test_missing_psc_rows_yields_empty_sample_qc_and_sample_names(self) -> None:
        text = _bcftools_stats_text(include_psc=False)
        parsed = qc.parse_bcftools_stats(text, Path("synthetic.tsv"))

        self.assertEqual(qc.build_sample_qc_rows("cohort", "raw", "snp", parsed), [])
        variant_rows = {row[3]: row[4] for row in qc.build_variant_qc_rows("cohort", "raw", "snp", parsed)}
        self.assertEqual(variant_rows["sample_names"], "")
        # cohort_total_genotypes comes from the SN section, independent
        # of how many PSC rows were actually present.
        self.assertEqual(variant_rows["cohort_total_genotypes"], "4")

    def test_raw_snp_fixture_matches_captured_baseline_values(self) -> None:
        # These values were captured from a real clean Docker smoke test
        # run of the (pre-refactor) AWK implementation; see the PR
        # description for the byte-level comparison across all 7
        # stage/type combinations.
        source = FIXTURES_DIR / "raw_snp.bcftools.stats.tsv"
        parsed = qc.parse_bcftools_stats(source.read_text(encoding="utf-8"), source)

        variant_qc = {row[3]: row[4] for row in qc.build_variant_qc_rows("cohort", "raw", "snp", parsed)}
        self.assertEqual(variant_qc["cohort_missing_genotypes"], "0")
        self.assertEqual(variant_qc["cohort_total_genotypes"], "4")
        self.assertEqual(variant_qc["cohort_missingness_rate"], "0.000000")
        self.assertEqual(variant_qc["sample_names"], "sample_a,sample_b")

        sample_qc = qc.build_sample_qc_rows("cohort", "raw", "snp", parsed)
        self.assertEqual(
            sample_qc[0],
            ["cohort", "raw", "snp", "sample_a", "1", "1", "0", "0", "0.000000", "8.0", "1"],
        )
        self.assertEqual(
            sample_qc[1],
            ["cohort", "raw", "snp", "sample_b", "1", "1", "0", "0", "0.000000", "5.5", "1"],
        )

        summary_text = qc.build_summary_text("cohort", "raw", "snp", parsed)
        self.assertIn("Missing genotypes: 0/4 (0.000000)", summary_text)


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_all_three_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            variant_qc_output = tmp_path / "variant_qc.tsv"
            sample_qc_output = tmp_path / "sample_qc.tsv"
            summary_output = tmp_path / "summary.txt"

            exit_code = qc.main(
                [
                    "--bcftools-stats",
                    str(FIXTURES_DIR / "raw_snp.bcftools.stats.tsv"),
                    "--cohort-id",
                    "cohort",
                    "--stage",
                    "raw",
                    "--variant-type",
                    "snp",
                    "--variant-qc-output",
                    str(variant_qc_output),
                    "--sample-qc-output",
                    str(sample_qc_output),
                    "--summary-output",
                    str(summary_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(variant_qc_output.exists())
            self.assertTrue(sample_qc_output.exists())
            self.assertTrue(summary_output.exists())

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_input = tmp_path / "does-not-exist.bcftools.stats.tsv"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = qc.main(
                    [
                        "--bcftools-stats",
                        str(missing_input),
                        "--cohort-id",
                        "cohort",
                        "--stage",
                        "raw",
                        "--variant-type",
                        "snp",
                        "--variant-qc-output",
                        str(tmp_path / "variant_qc.tsv"),
                        "--sample-qc-output",
                        str(tmp_path / "sample_qc.tsv"),
                        "--summary-output",
                        str(tmp_path / "summary.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_input), stderr.getvalue())

    def test_main_fails_clearly_for_malformed_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            broken_input = tmp_path / "broken.bcftools.stats.tsv"
            broken_input.write_text("SN\t0\tnumber of records:\t2\n", encoding="utf-8")
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = qc.main(
                    [
                        "--bcftools-stats",
                        str(broken_input),
                        "--cohort-id",
                        "cohort",
                        "--stage",
                        "raw",
                        "--variant-type",
                        "snp",
                        "--variant-qc-output",
                        str(tmp_path / "variant_qc.tsv"),
                        "--sample-qc-output",
                        str(tmp_path / "sample_qc.tsv"),
                        "--summary-output",
                        str(tmp_path / "summary.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(broken_input), stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            variant_qc_output = tmp_path / "variant_qc.tsv"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--bcftools-stats",
                    str(FIXTURES_DIR / "empty_indel.bcftools.stats.tsv"),
                    "--cohort-id",
                    "cohort",
                    "--stage",
                    "raw",
                    "--variant-type",
                    "indel",
                    "--variant-qc-output",
                    str(variant_qc_output),
                    "--sample-qc-output",
                    str(tmp_path / "sample_qc.tsv"),
                    "--summary-output",
                    str(tmp_path / "summary.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(variant_qc_output.exists())


if __name__ == "__main__":
    unittest.main()
