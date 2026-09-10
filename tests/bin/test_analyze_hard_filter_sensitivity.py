"""Tests for bin/analyze_hard_filter_sensitivity.py."""

from __future__ import annotations

import csv
import gzip
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "analyze_hard_filter_sensitivity.py"
SCENARIOS = REPO_ROOT / "conf" / "hard_filter_sensitivity_scenarios.json"


class SensitivityCliTests(unittest.TestCase):
    def test_current_audit_detects_cancelling_record_mismatches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "swapped.vcf.gz",
                [
                    "chr1\t1\t.\tA\tT\t50\tPASS\tQD=1",
                    "chr1\t2\t.\tA\tT\t50\tSNP_QD_LOW\tQD=10",
                ],
            )
            result, _, sensitivity, _ = _run(vcf, directory, "snp")
            self.assertEqual(0, result.returncode, result.stderr)
            current = {
                r["annotation"]: r for r in _read_tsv(sensitivity) if r["scenario"] == "current"
            }
            for annotation in ("QD", "ANY_FILTER"):
                self.assertEqual("0", current[annotation]["predicted_minus_observed"])
                self.assertEqual("2", current[annotation]["discordant_records"])

    def test_wrong_variant_type_filter_tag_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "wrong-type.vcf.gz",
                [
                    "chr1\t1\t.\tA\tAT\t50\tINDEL_QD_LOW\tQD=1",
                ],
            )
            result, distribution, sensitivity, summary = _run(vcf, directory, "snp")
            self.assertEqual(1, result.returncode)
            self.assertIn("unexpected FILTER", result.stderr)
            self.assertFalse(any(p.exists() for p in (distribution, sensitivity, summary)))

    def test_output_failure_does_not_leave_partial_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(directory / "empty.vcf.gz", [])
            (directory / "summary.txt").mkdir()
            result, distribution, sensitivity, _ = _run(vcf, directory, "snp")
            self.assertEqual(1, result.returncode)
            self.assertFalse(distribution.exists())
            self.assertFalse(sensitivity.exists())

    def test_snp_counts_current_thresholds_and_missing_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "snp.vcf.gz",
                [
                    "chr1\t1\t.\tA\tT\t50\tPASS\tQD=5;SOR=1;FS=10;MQ=60;MQRankSum=0;ReadPosRankSum=0",
                    "chr1\t2\t.\tA\tT\t20\tSNP_QD_LOW;SNP_SOR_HIGH\tQD=1;SOR=4;FS=70;MQ=30;MQRankSum=-13;ReadPosRankSum=-9",
                    "chr1\t3\t.\tA\tT\t.\tPASS\tQD=.;SOR=.;FS=.;MQ=.;MQRankSum=.;ReadPosRankSum=.",
                ],
            )
            result, distribution, sensitivity, summary = _run(vcf, directory, "snp")

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = _read_tsv(sensitivity)
            current = {row["annotation"]: row for row in rows if row["scenario"] == "current"}
            self.assertEqual(current["QD"]["hit_records"], "1")
            self.assertEqual(current["QD"]["missing_records"], "1")
            self.assertEqual(current["QD"]["predicted_minus_observed"], "0")
            self.assertEqual(current["FS"]["hit_records"], "1")
            self.assertEqual(current["ANY_FILTER"]["hit_records"], "1")
            self.assertEqual(current["ANY_FILTER"]["predicted_minus_observed"], "0")

            summary_rows = [
                row
                for row in _read_tsv(distribution)
                if row["row_type"] == "summary" and row["annotation"] == "QUAL"
            ]
            self.assertEqual(summary_rows[0]["present_records"], "2")
            self.assertEqual(summary_rows[0]["missing_records"], "1")
            self.assertIn("do not estimate accuracy", summary.read_text(encoding="utf-8"))

    def test_any_filter_denominator_excludes_wholly_unevaluable_records(self) -> None:
        # Two of these three records carry no filtered annotation at all,
        # so no threshold can ever hit them. Reporting them as evaluable
        # in the union row would restate an unevaluated record as one that
        # passed every filter -- the separation Issue #46 requires between
        # missing/evaluable and threshold-hit accounting.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "snp.vcf.gz",
                [
                    "chr1\t1\t.\tA\tT\t.\tPASS\tAC=1",
                    "chr1\t2\t.\tA\tT\t.\tPASS\tAC=1",
                    "chr1\t3\t.\tA\tT\t500\tPASS\t"
                    "QD=30;SOR=0.5;FS=0;MQ=60;MQRankSum=0;ReadPosRankSum=0",
                ],
            )
            result, _distribution, sensitivity, _summary = _run(vcf, directory, "snp")

            self.assertEqual(result.returncode, 0, result.stderr)
            any_filter = next(
                row
                for row in _read_tsv(sensitivity)
                if row["scenario"] == "current" and row["annotation"] == "ANY_FILTER"
            )
            self.assertEqual(any_filter["total_records"], "3")
            self.assertEqual(any_filter["present_records"], "1")
            self.assertEqual(any_filter["missing_records"], "2")
            # 0 of the 1 evaluable record is hit, which is a real rate --
            # distinct from the 0/3 that the total-record denominator
            # would have implied.
            self.assertEqual(any_filter["hit_records"], "0")
            self.assertEqual(any_filter["hit_rate_among_present"], "0.000000")
            self.assertEqual(any_filter["hit_rate_among_total"], "0.000000")

    def test_any_filter_present_plus_missing_equals_total(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "snp.vcf.gz",
                [
                    "chr1\t1\t.\tA\tT\t50\tPASS\tQD=5",
                    "chr1\t2\t.\tA\tT\t.\tPASS\tAC=1",
                ],
            )
            result, _distribution, sensitivity, _summary = _run(vcf, directory, "snp")

            self.assertEqual(result.returncode, 0, result.stderr)
            for row in _read_tsv(sensitivity):
                if row["annotation"] != "ANY_FILTER":
                    continue
                with self.subTest(scenario=row["scenario"]):
                    self.assertEqual(
                        int(row["present_records"]) + int(row["missing_records"]),
                        int(row["total_records"]),
                    )

    def test_boundary_values_do_not_hit_strict_current_comparisons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "boundary.vcf.gz",
                [
                    "chr1\t1\t.\tA\tAT\t30\tPASS\tQD=2;FS=200;ReadPosRankSum=-20",
                ],
            )
            result, _distribution, sensitivity, _summary = _run(vcf, directory, "indel")

            self.assertEqual(result.returncode, 0, result.stderr)
            current = [row for row in _read_tsv(sensitivity) if row["scenario"] == "current"]
            self.assertTrue(all(row["hit_records"] == "0" for row in current))

    def test_non_numeric_annotation_fails_before_outputs_are_written(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "broken.vcf.gz",
                ["chr1\t1\t.\tA\tT\t50\tPASS\tQD=not-a-number"],
            )
            result, distribution, sensitivity, summary = _run(vcf, directory, "snp")

            self.assertEqual(result.returncode, 1)
            self.assertIn("not numeric", result.stderr)
            self.assertTrue(all(not path.exists() for path in (distribution, sensitivity, summary)))

    def test_default_scenarios_are_versioned_and_complete(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(directory / "empty.vcf.gz", [])
            result, _distribution, sensitivity, _summary = _run(vcf, directory, "snp")

            self.assertEqual(result.returncode, 0, result.stderr)
            scenario_names = {row["scenario"] for row in _read_tsv(sensitivity)}
            self.assertEqual(scenario_names, {"lenient", "current", "stringent"})

    def test_current_scenario_matches_nextflow_defaults(self) -> None:
        payload = json.loads(SCENARIOS.read_text(encoding="utf-8"))
        current = payload["scenarios"]["current"]
        config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
        parameter_names = {
            "snp": {
                "QD": "snp_filter_qd_min",
                "QUAL": "snp_filter_qual_min",
                "SOR": "snp_filter_sor_max",
                "FS": "snp_filter_fs_max",
                "MQ": "snp_filter_mq_min",
                "MQRankSum": "snp_filter_mq_rank_sum_min",
                "ReadPosRankSum": "snp_filter_read_pos_rank_sum_min",
            },
            "indel": {
                "QD": "indel_filter_qd_min",
                "QUAL": "indel_filter_qual_min",
                "FS": "indel_filter_fs_max",
                "ReadPosRankSum": "indel_filter_read_pos_rank_sum_min",
            },
        }

        for variant_type, annotation_parameters in parameter_names.items():
            for annotation, parameter in annotation_parameters.items():
                with self.subTest(parameter=parameter):
                    match = re.search(
                        rf"^\s*{re.escape(parameter)}\s*=\s*(-?[0-9]+(?:\.[0-9]+)?)\s*$",
                        config,
                        re.MULTILINE,
                    )
                    self.assertIsNotNone(match)
                    self.assertEqual(float(match.group(1)), current[variant_type][annotation])

    def test_indel_distribution_still_reports_unfiltered_annotations(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "indel.vcf.gz",
                [
                    "chr1\t1\t.\tA\tAT\t50\tPASS\tQD=5;SOR=2;FS=10;MQ=60;MQRankSum=0;ReadPosRankSum=0"
                ],
            )
            result, distribution, _sensitivity, _summary = _run(vcf, directory, "indel")

            self.assertEqual(result.returncode, 0, result.stderr)
            summaries = {
                row["annotation"]: row
                for row in _read_tsv(distribution)
                if row["row_type"] == "summary"
            }
            self.assertEqual(
                set(summaries),
                {"QD", "QUAL", "SOR", "FS", "MQ", "MQRankSum", "ReadPosRankSum"},
            )
            self.assertEqual(summaries["SOR"]["present_records"], "1")


def _write_vcf(path: Path, rows: list[str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            handle.write(row + "\n")
    return path


def _run(vcf: Path, directory: Path, variant_type: str):
    distribution = directory / "distribution.tsv"
    sensitivity = directory / "sensitivity.tsv"
    summary = directory / "summary.txt"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--filtered-vcf",
            str(vcf),
            "--cohort-id",
            "cohort",
            "--variant-type",
            variant_type,
            "--scenario-config",
            str(SCENARIOS),
            "--distribution-output",
            str(distribution),
            "--sensitivity-output",
            str(sensitivity),
            "--summary-output",
            str(summary),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    return result, distribution, sensitivity, summary


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


if __name__ == "__main__":
    unittest.main()
