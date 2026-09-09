"""Tests for bin/prepare_joint_genotyping_benchmark.py."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "prepare_joint_genotyping_benchmark.py"
CONTAINER = "broadinstitute/gatk:4.6.2.0@sha256:" + "a" * 64
COMMIT = "b" * 40


class BenchmarkPreparationTests(unittest.TestCase):
    def test_51_inputs_produce_sample_map_interval_candidates_and_pending_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            result, sample_map, interval_plan, evidence_path = _run(
                directory, fai, manifest
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with sample_map.open(encoding="utf-8", newline="") as handle:
                sample_rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(len(sample_rows), 51)
            self.assertEqual(sample_rows[0][0], "sample_000")
            self.assertEqual(len(sample_rows[0]), 3)
            self.assertTrue(sample_rows[0][1].endswith("sample_000.g.vcf.gz"))
            self.assertTrue(sample_rows[0][2].endswith("sample_000.g.vcf.gz.tbi"))

            interval_rows = _read_tsv(interval_plan)
            plan_ids = {row["plan_id"] for row in interval_rows}
            self.assertEqual(plan_ids, {"baseline_per_contig", "candidate_split_group"})
            candidate = [
                row
                for row in interval_rows
                if row["plan_id"] == "candidate_split_group"
            ]
            self.assertTrue(any("window" in row["strategy"] for row in candidate))
            small = next(
                row
                for row in candidate
                if row["interval_group_id"] == "small_scaffolds"
            )
            self.assertEqual(json.loads(small["intervals_json"]), ["small1", "small2"])

            evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
            self.assertEqual(evidence["status"], "PENDING_REAL_BENCHMARK")
            self.assertEqual(evidence["inputs"]["sample_count"], 51)
            self.assertEqual(len(evidence["experiments"]), 5)
            self.assertTrue(
                all(
                    experiment["measurements"]["peak_rss_bytes"] is None
                    for experiment in evidence["experiments"]
                )
            )
            self.assertEqual(evidence["decision"]["full_327_sample_gate"], "PENDING")

    def test_50_inputs_are_rejected_before_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 50)
            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("need at least 51", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_sample_count_must_exceed_selected_batch_size(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            result, *outputs = _run(
                directory,
                fai,
                manifest,
                extra_args=["--batch-size", "51"],
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("need at least 52", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_manifest_and_header_sample_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            rows = list(
                csv.DictReader(manifest.read_text(encoding="utf-8").splitlines())
            )
            rows[0]["sample_id"] = "wrong"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("does not match manifest", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_reference_dictionary_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            fai.write_text("chr1\t40000001\nsmall1\t500\n", encoding="utf-8")

            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("does not match reference FAI", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))


def _build_inputs(directory: Path, sample_count: int) -> tuple[Path, Path]:
    contigs = (("chr1", 40_000_001), ("small1", 500), ("small2", 700))
    fai = directory / "reference.fa.fai"
    fai.write_text(
        "".join(f"{name}\t{length}\t0\t0\t0\n" for name, length in contigs),
        encoding="utf-8",
    )
    manifest = directory / "gvcfs.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("sample_id", "gvcf", "gvcf_index"))
        for index in range(sample_count):
            sample = f"sample_{index:03d}"
            gvcf = directory / f"{sample}.g.vcf.gz"
            with gzip.open(gvcf, "wt", encoding="utf-8") as output:
                output.write("##fileformat=VCFv4.2\n")
                for name, length in contigs:
                    output.write(f"##contig=<ID={name},length={length}>\n")
                output.write(
                    f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
                )
            index_path = Path(str(gvcf) + ".tbi")
            index_path.write_bytes(f"index-{sample}".encode())
            writer.writerow((sample, gvcf, index_path))
    return fai, manifest


def _run(
    directory: Path,
    fai: Path,
    manifest: Path,
    extra_args: list[str] | None = None,
):
    sample_map = directory / "sample_name_map.tsv"
    interval_plan = directory / "interval_plan.tsv"
    evidence = directory / "evidence.json"
    arguments = [
        sys.executable,
        str(SCRIPT),
        "--gvcf-manifest",
        str(manifest),
        "--reference-fai",
        str(fai),
        "--pipeline-commit",
        COMMIT,
        "--gatk-container",
        CONTAINER,
        "--sample-name-map-output",
        str(sample_map),
        "--interval-plan-output",
        str(interval_plan),
        "--evidence-template-output",
        str(evidence),
    ]
    if extra_args:
        arguments.extend(extra_args)
    result = subprocess.run(
        arguments,
        capture_output=True,
        text=True,
        check=False,
    )
    return result, sample_map, interval_plan, evidence


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


if __name__ == "__main__":
    unittest.main()
