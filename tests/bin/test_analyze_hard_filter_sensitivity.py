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
            current = {
                row["annotation"]: row for row in rows if row["scenario"] == "current"
            }
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
            self.assertIn(
                "do not estimate accuracy", summary.read_text(encoding="utf-8")
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
            current = [
                row for row in _read_tsv(sensitivity) if row["scenario"] == "current"
            ]
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
            self.assertTrue(
                all(not path.exists() for path in (distribution, sensitivity, summary))
            )

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
                    self.assertEqual(
                        float(match.group(1)), current[variant_type][annotation]
                    )

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
