#!/usr/bin/env python3
"""Bounded synthetic failure/resume and corrupt/missing-cache regression.

Uses only bundled reads. All writes and the single deliberate cache-output
deletion occur under a new empty --output-dir; never point it at a real run.
"""

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin"))
from resume_cache_inventory import snapshot, verify  # noqa: E402


def payload(directory: Path) -> dict:
    result = {}
    for relative in (
        "variants/raw/cohort.raw.vcf.gz",
        "gs_panel/cohort.gs_panel.genotype_matrix.tsv.gz",
    ):
        with gzip.open(directory / relative, "rt") as handle:
            result[relative] = [line.rstrip("\n") for line in handle if not line.startswith("##")]
    return result


def run_pipeline(
    base: Path, name: str, samples: Path, failure_config: Path | None = None, resume: bool = False
) -> tuple[int, float, Path]:
    launch = base / name
    launch.mkdir(exist_ok=True)
    trace = launch / ("resume.trace.txt" if resume else "trace.txt")
    command = [
        "nextflow",
        "run",
        str(REPO),
        "-profile",
        "test,docker",
        "-work-dir",
        str(launch / "work"),
        "--input",
        str(samples),
        "--outdir",
        str(launch / "results"),
        "-with-trace",
        str(trace),
    ]
    if failure_config:
        command += ["-c", str(failure_config)]
    if resume:
        command += ["-resume"]
    start = time.monotonic()
    with (launch / ("resume.log" if resume else "run.log")).open("w") as log:
        result = subprocess.run(
            command,
            cwd=launch,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "NXF_VER": "26.04.6"},
            check=False,
        )
    return result.returncode, time.monotonic() - start, trace


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    base = args.output_dir.resolve()
    if base.exists() and any(base.iterdir()):
        parser.error("output-dir must be new or empty; existing runs are never altered")
    base.mkdir(parents=True, exist_ok=True)
    # Mutated-input regression below must touch only disposable copies.
    shutil.copytree(REPO / "tests/data/reads", base / "reads")
    samples = base / "samples.csv"
    samples.write_text(
        (REPO / "tests/data/samplesheet.csv")
        .read_text()
        .replace("tests/data/reads/", str(base / "reads") + "/")
    )
    failure = base / "stop.config"
    failure.write_text(
        "process { withName: BUILD_GS_PANEL { beforeScript = 'exit 73'; errorStrategy = 'terminate' } }\n"
    )
    reference_code, reference_time, _ = run_pipeline(base, "baseline", samples)
    if reference_code:
        parser.exit(1, "baseline failed; inspect baseline/run.log\n")
    failed_code, failed_time, _ = run_pipeline(base, "interrupted", samples, failure_config=failure)
    launch = base / "interrupted"
    if failed_code == 0 or (launch / "results/provenance/cohort.run_manifest.json").exists():
        parser.exit(1, "controlled failure was not observed or produced a completion manifest\n")
    inventory = snapshot(launch / "work")
    if verify(launch / "work", inventory):
        parser.exit(1, "stopped successful-task cache is already inconsistent\n")
    code, resume_time, trace = run_pipeline(base, "interrupted", samples, resume=True)
    if code or payload(launch / "results") != payload(base / "baseline/results"):
        parser.exit(1, "resume failed or changed biological payload\n")
    with trace.open() as handle:
        cached = sum(
            row.get("status") == "CACHED" for row in csv.DictReader(handle, delimiter="\t")
        )
    if cached == 0:
        parser.exit(1, "resume did not exercise a cached task\n")
    inventory = snapshot(launch / "work")
    candidates = [
        row
        for row in inventory["files"]
        if row["path"].endswith("cohort.gs_panel.genotype_matrix.tsv.gz")
    ]
    if len(candidates) != 1:
        parser.exit(1, "could not uniquely identify the synthetic matrix cache output\n")
    target = launch / "work" / candidates[0]["path"]
    if target.is_symlink() or not target.resolve().is_relative_to(base):
        parser.exit(1, "cache mutation target escaped the isolated test directory\n")
    with target.open("r+b") as handle:
        first = handle.read(1)
        handle.seek(0)
        handle.write(bytes([first[0] ^ 1]))
    if not any(row["reason"] == "changed" for row in verify(launch / "work", inventory)):
        parser.exit(1, "same-length corruption was not detected\n")
    target.unlink()  # Only the just-created, isolated synthetic test cache file.
    if not any(row["reason"] == "missing" for row in verify(launch / "work", inventory)):
        parser.exit(1, "missing cache output was not detected\n")
    code, missing_time, _ = run_pipeline(base, "interrupted", samples, resume=True)
    if code or payload(launch / "results") != payload(base / "baseline/results"):
        parser.exit(1, "missing-cache re-execution failed or changed biological payload\n")
    plan_path = launch / "results/provenance/cohort.joint_plan.json"
    original_binding = json.loads(plan_path.read_text())["input_summary_sha256"]
    raw_input = sorted((base / "reads").glob("*.fastq.gz"))[0]
    old_stat = raw_input.stat()
    # Change gzip's timestamp header, leaving the FASTQ reads and compressed
    # size unchanged; restore filesystem mtime to defeat metadata-only caching.
    with raw_input.open("r+b") as handle:
        if handle.read(2) != b"\x1f\x8b":
            parser.exit(1, "synthetic input is not gzip\n")
        handle.seek(4)
        value = handle.read(1)
        handle.seek(4)
        handle.write(bytes([value[0] ^ 1]))
    os.utime(raw_input, ns=(old_stat.st_atime_ns, old_stat.st_mtime_ns))
    code, changed_input_time, trace = run_pipeline(base, "interrupted", samples, resume=True)
    if code or payload(launch / "results") != payload(base / "baseline/results"):
        parser.exit(1, "same-size changed-input resume failed or changed biological payload\n")
    if json.loads(plan_path.read_text())["input_summary_sha256"] == original_binding:
        parser.exit(1, "changed FASTQ bytes retained the old approval input binding\n")
    with trace.open() as handle:
        fastp_rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if ":FASTP " in row.get("name", "")
        ]
    if not fastp_rows or any(row.get("status") != "COMPLETED" for row in fastp_rows):
        parser.exit(1, "scientific processing reused metadata-only cache for changed FASTQ bytes\n")
    (base / "resume_evidence.json").write_text(
        json.dumps(
            {
                "scope": "bundled_synthetic_only",
                "controlled_failure_exit": failed_code,
                "cached_tasks_after_resume": cached,
                "biological_payload_equivalent": True,
                "same_length_corruption_rejected": True,
                "missing_output_reexecuted": True,
                "same_size_mtime_changed_input_rebound_and_reexecuted": True,
                "wall_seconds": {
                    "baseline": reference_time,
                    "failed": failed_time,
                    "resume": resume_time,
                    "missing_cache_resume": missing_time,
                    "changed_input_resume": changed_input_time,
                },
                "real_51_sample_validation": "not_run",
                "production_sla": "not_established",
            },
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
