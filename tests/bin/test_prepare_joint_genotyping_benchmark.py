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
    def test_sample_map_uses_literal_tabs_without_csv_quoting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / 'quoted"directory'
            directory.mkdir()
            fai, manifest = _build_inputs(directory, 51)
            result, sample_map, _, _ = _run(directory, fai, manifest)
            self.assertEqual(0, result.returncode, result.stderr)
            first = sample_map.read_text().splitlines()[0].split("\t")
            self.assertEqual(str(directory / "sample_000.g.vcf.gz"), first[1])
            self.assertEqual(str(directory / "sample_000.g.vcf.gz.tbi"), first[2])

    def test_manifest_sample_names_are_not_silently_trimmed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            manifest.write_text(manifest.read_text().replace("sample_000,", " sample_000,"))
            result, *outputs = _run(directory, fai, manifest)
            self.assertEqual(1, result.returncode)
            self.assertFalse(any(p.exists() for p in outputs))

    def test_extra_manifest_column_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            lines = manifest.read_text().splitlines()
            lines[1] += ",unaccounted-source"
            manifest.write_text("\n".join(lines) + "\n")
            result, *outputs = _run(directory, fai, manifest)
            self.assertEqual(1, result.returncode)
            self.assertFalse(any(p.exists() for p in outputs))

    def test_fake_container_digest_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            result, *outputs = _run(
                directory, fai, manifest, ["--gatk-container", "gatk@sha256:fake"]
            )
            self.assertEqual(1, result.returncode)
            self.assertFalse(any(p.exists() for p in outputs))

    def test_grouped_scaffolds_preserve_reference_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            contigs = (("small1", 500), ("chr1", 40000001), ("small2", 700))
            fai, manifest = _build_inputs(directory, 51, contigs)
            result, _, intervals, _ = _run(directory, fai, manifest)
            self.assertEqual(0, result.returncode, result.stderr)
            names = []
            for row in _read_tsv(intervals):
                if row["plan_id"] != "candidate_split_group":
                    continue
                for interval in json.loads(row["intervals_json"]):
                    name = interval.split(":")[0]
                    if not names or names[-1] != name:
                        names.append(name)
            self.assertEqual([c[0] for c in contigs], names)

    def test_51_inputs_produce_sample_map_interval_candidates_and_pending_evidence(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            result, sample_map, interval_plan, evidence_path = _run(directory, fai, manifest)

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
            self.assertEqual(
                plan_ids, {"baseline_per_contig", "candidate_split_only", "candidate_split_group"}
            )
            candidate = [row for row in interval_rows if row["plan_id"] == "candidate_split_group"]
            self.assertTrue(any("window" in row["strategy"] for row in candidate))
            small = next(row for row in candidate if row["interval_group_id"] == "small_scaffolds")
            self.assertEqual(json.loads(small["intervals_json"]), ["small1", "small2"])

            evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
            self.assertEqual(evidence["status"], "PENDING_REAL_BENCHMARK")
            self.assertEqual(evidence["inputs"]["sample_count"], 51)
            self.assertEqual(len(evidence["experiments"]), 6)
            self.assertFalse(evidence["input_validation"]["benchmark_ready"])
            self.assertFalse(evidence["input_validation"]["lineage_verified"])
            experiments = {e["experiment_id"]: e for e in evidence["experiments"]}
            for experiment in experiments.values():
                if experiment["compare_to"] is None:
                    continue
                parent = experiments[experiment["compare_to"]]
                factors = ("gvcf_input_mode", "interval_plan_id", "reblock_gvcfs", "consolidate")
                self.assertEqual(1, sum(experiment[key] != parent[key] for key in factors))
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
            rows = list(csv.DictReader(manifest.read_text(encoding="utf-8").splitlines()))
            rows[0]["sample_id"] = "wrong"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("does not match manifest", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_index_belonging_to_another_sample_is_rejected(self) -> None:
        # Every other check passes for this manifest: both files exist,
        # no path is reused, and each header sample still matches its own
        # row. Only the gVCF/index pairing is wrong, and it would reach
        # GenomicsDBImport through the sample-name-map's third column.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            rows = list(csv.DictReader(manifest.read_text(encoding="utf-8").splitlines()))
            shifted = [row["gvcf_index"] for row in rows[1:]] + [rows[0]["gvcf_index"]]
            for row, index_path in zip(rows, shifted, strict=True):
                row["gvcf_index"] = index_path
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("must be the gVCF's own index", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_small_scaffold_bound_above_the_window_size_is_rejected(self) -> None:
        # Inverted bounds would group every contig, chromosome-scale ones
        # included, into the single "small_scaffolds" task, so E2 would
        # silently measure a one-interval plan instead of the documented
        # chromosome split.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)

            result, *outputs = _run(
                directory,
                fai,
                manifest,
                ["--window-size-bp", "1000000", "--small-scaffold-max-bp", "50000000"],
            )

            self.assertEqual(result.returncode, 1)
            self.assertIn("--small-scaffold-max-bp must be smaller", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))

    def test_candidate_plan_tiles_every_contig_exactly_once(self) -> None:
        # The split/group candidate must stay a partition of the
        # reference: a gap would drop variants from the benchmark and an
        # overlap would import the same region twice.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)

            result, _sample_map, interval_plan, _evidence = _run(directory, fai, manifest)
            self.assertEqual(result.returncode, 0, result.stderr)

            lengths = {}
            for line in fai.read_text(encoding="utf-8").splitlines():
                name, length = line.split("\t")[:2]
                lengths[name] = int(length)

            covered: dict[str, list[tuple[int, int]]] = {name: [] for name in lengths}
            for row in _read_tsv(interval_plan):
                if row["plan_id"] != "candidate_split_group":
                    continue
                for interval in json.loads(row["intervals_json"]):
                    if ":" in interval:
                        contig, span = interval.split(":")
                        start, end = span.split("-")
                        covered[contig].append((int(start), int(end)))
                    else:
                        covered[interval].append((1, lengths[interval]))

            for name, length in lengths.items():
                with self.subTest(contig=name):
                    spans = sorted(covered[name])
                    self.assertEqual(1, spans[0][0])
                    self.assertEqual(length, spans[-1][1])
                    self.assertEqual(sum(end - start + 1 for start, end in spans), length)
                    for earlier, later in zip(spans, spans[1:], strict=False):
                        self.assertEqual(earlier[1] + 1, later[0])

    def test_reference_dictionary_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fai, manifest = _build_inputs(directory, 51)
            fai.write_text("chr1\t40000001\nsmall1\t500\n", encoding="utf-8")

            result, *outputs = _run(directory, fai, manifest)

            self.assertEqual(result.returncode, 1)
            self.assertIn("does not match reference FAI", result.stderr)
            self.assertTrue(all(not output.exists() for output in outputs))


def _build_inputs(directory: Path, sample_count: int, contigs=None) -> tuple[Path, Path]:
    if contigs is None:
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
                output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
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
