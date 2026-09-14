#!/usr/bin/env python3
"""Measure how one fixed interval's cost moves with cohort size, at 13/26/51 samples.

Auxiliary to the E0-E4 suite, never an input to its one-factor decisions. Its
only purpose is to reduce the sample-count uncertainty in a 327-sample
projection, because every E0-E4 figure was measured at exactly 51 samples and
says nothing about that axis.

Subsets are nested (13 subset of 26 subset of 51) and taken as a deterministic
prefix of the validated manifest's own sample order, so no sample is chosen for
its result. Reference, GATK image, interval, heaps, CPUs, batch size and
concurrency are identical across the three levels; only cohort size moves.

Batching regimes differ by construction and must not be fitted as one line.
With batch-size 50, 13 and 26 samples import as a single batch while 51 imports
as two (50 + 1). GenomicsDBImport is therefore reported per regime. GenotypeGVCFs
reads a finished workspace and has no batching, so its three points share one
regime and are the ones a trend may be fitted to.

No FASTQ or gVCF is generated. Existing validated gVCFs are read only.
"""

import argparse
import datetime
import json
import shutil
import subprocess
import threading
import time
from pathlib import Path

from benchmark_tools import GATK, ToolRunner, batches_from_log, filesystem_metrics, joint_integrity
from execute_experiments import (
    GENOMICSDB_JAVA_HEAP_MIB,
    GENOMICSDB_MEMORY_GIB,
    GENOTYPE_MEMORY_GIB,
    GENOTYPE_NATIVE_RESERVE_GIB,
)
from run_generation import memory
from stage_generation import sha256

SUBSET_SIZES = (13, 26, 51)
BATCH_SIZE = 50


class HostMonitor:
    """Sample host memory and swap for the whole measurement, as the suite does."""

    def __init__(self, path):
        self.path = path
        self.done = threading.Event()
        self.thread = None
        self.samples = []

    def __enter__(self):
        self.thread = threading.Thread(target=self.run, daemon=True)
        self.thread.start()
        return self

    def run(self):
        with self.path.open("w") as handle:
            handle.write("timestamp\tmem_available_bytes\tswap_bytes\n")
            while not self.done.is_set():
                available, swap = memory()
                self.samples.append((available, swap))
                handle.write(
                    f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{available}\t{swap}\n"
                )
                handle.flush()
                self.done.wait(5)

    def __exit__(self, *_):
        self.done.set()
        if self.thread is not None:
            self.thread.join(timeout=30)
        return False

    def summary(self):
        if not self.samples:
            return {}
        memories = [row[0] for row in self.samples]
        swaps = [row[1] for row in self.samples]
        return {
            "monitor_samples": len(self.samples),
            "mem_available_min_bytes": min(memories),
            "swap_start_bytes": swaps[0],
            "swap_peak_bytes": max(swaps),
            "swap_end_bytes": swaps[-1],
            "swap_delta_bytes": swaps[-1] - swaps[0],
        }


def measure_level(runner, root, level_dir, samples, reference_name, contigs, interval, monitor):
    level_dir.mkdir(parents=True, exist_ok=False)
    sample_map = level_dir / "sample_name_map.tsv"
    sample_map.write_text(
        "".join(
            f"{row['sample']}\t/input/{row['gvcf_name']}\t/input/{row['gvcf_name']}.tbi\n"
            for row in samples
        )
    )
    heap = runner.java_heap_gib(GENOTYPE_MEMORY_GIB, GENOTYPE_NATIVE_RESERVE_GIB)
    started = time.monotonic()
    imported = runner.run(
        level_dir,
        "GenomicsDBImport",
        [
            "gatk",
            "--java-options",
            f"-Xmx{GENOMICSDB_JAVA_HEAP_MIB}m",
            "GenomicsDBImport",
            "--sample-name-map",
            "/bench/" + str(sample_map.relative_to(root)),
            "--validate-sample-name-map",
            "true",
            "--genomicsdb-workspace-path",
            "workspace",
            "--intervals",
            interval,
            "--reader-threads",
            "8",
            "--batch-size",
            str(BATCH_SIZE),
            "--consolidate",
            "false",
            "--tmp-dir",
            ".",
        ],
        memory_gib=GENOMICSDB_MEMORY_GIB,
    )
    batching = batches_from_log(
        (level_dir / "GenomicsDBImport.log").read_text(), len(samples), BATCH_SIZE
    )
    workspace = filesystem_metrics(level_dir / "workspace")
    genotyped = runner.run(
        level_dir,
        "GenotypeGVCFs",
        [
            "gatk",
            "--java-options",
            f"-Xmx{heap}g",
            "GenotypeGVCFs",
            "--reference",
            "/reference/" + reference_name,
            "--variant",
            "gendb://workspace",
            "--intervals",
            interval,
            "--sample-ploidy",
            "2",
            "--only-output-calls-starting-in-intervals",
            "true",
            "--output",
            "joint.vcf.gz",
            "--create-output-variant-index",
            "true",
        ],
        memory_gib=GENOTYPE_MEMORY_GIB,
    )
    integrity = joint_integrity(
        level_dir / "joint.vcf.gz", [row["sample"] for row in samples], contigs
    )
    return {
        "sample_count": len(samples),
        "batch_size": BATCH_SIZE,
        "batching_regime": ("single_batch" if len(samples) <= BATCH_SIZE else "multi_batch"),
        "batches_observed": batching["observed"],
        "batches_completed": batching["completed"],
        "reader_initialization_fell_back_to_serial": batching[
            "reader_initialization_fell_back_to_serial"
        ],
        "import_wall_seconds": imported["wall_seconds"],
        "import_peak_rss_bytes": imported["peak_rss_bytes"],
        "import_exit_code": imported["exit_code"],
        "import_attempt": imported["attempt"],
        "workspace_bytes": workspace["bytes"],
        "workspace_file_count": workspace["file_count"],
        "genotype_wall_seconds": genotyped["wall_seconds"],
        "genotype_peak_rss_bytes": genotyped["peak_rss_bytes"],
        "genotype_exit_code": genotyped["exit_code"],
        "genotype_attempt": genotyped["attempt"],
        "output_sample_count": integrity["sample_count"],
        "output_sample_order_sha256": integrity["sample_order_sha256"],
        "output_variant_count": integrity["variant_count"],
        "output_called_alleles": integrity["called_alleles"],
        "output_accounting_sha256": integrity["accounting_sha256"],
        "sample_order_valid": integrity["sample_order_valid"],
        "contig_order_valid": integrity["contig_order_valid"],
        "allele_accounting_valid": integrity["allele_accounting_valid"],
        "level_wall_seconds": time.monotonic() - started,
        "host_during_level": monitor.summary(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validated-cohort", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--interval",
        required=True,
        help="one fixed 20 Mb chromosome window, identical at every cohort size",
    )
    args = parser.parse_args()

    validated = json.loads(args.validated_cohort.read_text())
    evidence = validated["evidence"]
    if not evidence.get("benchmark_ready") or not evidence.get("lineage_verified"):
        raise ValueError("Validated cohort gate is not open")
    samples = evidence["samples"]
    if len(samples) != 51 or len({row["sample"] for row in samples}) != 51:
        raise ValueError("Exactly 51 unique validated samples are required")

    input_dir = Path(validated["gvcf_directory"])
    for row in samples:
        gvcf = input_dir / row["gvcf_name"]
        if (
            sha256(gvcf) != row["gvcf_sha256"]
            or sha256(Path(str(gvcf) + ".tbi")) != row["gvcf_index_sha256"]
        ):
            raise ValueError("Validated gVCF or index changed before measurement")

    reference = validated["reference_paths"]
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
    from prepare_joint_genotyping_benchmark import read_reference_fai

    contigs = [(row.name, row.length) for row in read_reference_fai(Path(reference["fai"]))]

    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    stop = threading.Event()
    runner = ToolRunner(root, input_dir, Path(reference["fasta"]).parent, stop)

    record = {
        "purpose": (
            "auxiliary sample-count scaling measurement; not an input to any E0-E4 "
            "one-factor decision"
        ),
        "cohort_id": evidence["cohort_id"],
        "production_sha": evidence["production_sha"],
        "gatk_container": GATK,
        "reference_accession": evidence["reference_accession"],
        "interval": args.interval,
        "subset_sizes": list(SUBSET_SIZES),
        "subsets_are_nested": True,
        "subset_selection": (
            "deterministic prefix of the validated manifest's own sample order; "
            "no sample chosen for its result"
        ),
        "fixed_across_levels": {
            "genomicsdb_cpus": 8,
            "genomicsdb_memory_gib": GENOMICSDB_MEMORY_GIB,
            "genomicsdb_java_heap_mib": GENOMICSDB_JAVA_HEAP_MIB,
            "genotype_memory_gib": GENOTYPE_MEMORY_GIB,
            "genotype_java_heap_gib": runner.java_heap_gib(
                GENOTYPE_MEMORY_GIB, GENOTYPE_NATIVE_RESERVE_GIB
            ),
            "batch_size": BATCH_SIZE,
            "sample_ploidy": 2,
            "only_output_calls_starting_in_intervals": True,
            "concurrent_tasks": 1,
            "reader_threads_requested": 8,
        },
        "concurrency_note": (
            "every level runs alone at concurrency 1, so wall times are comparable "
            "to each other but not to the E0-E4 suite, which ran three tasks at once"
        ),
        "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "levels": [],
    }
    (root / "sample_scaling.private.json").write_text(json.dumps(record, indent=2) + "\n")

    ordered = list(samples)
    for size in SUBSET_SIZES:
        if shutil.disk_usage(root).free < 300_000_000_000:
            raise RuntimeError("Storage headroom guard failed before the next level")
        with HostMonitor(root / f"host_{size}.tsv") as monitor:
            level = measure_level(
                runner,
                root,
                root / f"samples_{size:03d}",
                ordered[:size],
                Path(reference["fasta"]).name,
                contigs,
                args.interval,
                monitor,
            )
        record["levels"].append(level)
        record["finished_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (root / "sample_scaling.private.json").write_text(json.dumps(record, indent=2) + "\n")
        print(
            f"level {size}: import {level['import_peak_rss_bytes'] / 2**30:.2f} GiB / "
            f"{level['import_wall_seconds']:.0f}s ({level['batching_regime']}), "
            f"genotype {level['genotype_peak_rss_bytes'] / 2**30:.2f} GiB / "
            f"{level['genotype_wall_seconds']:.0f}s",
            flush=True,
        )

    nested = [{row["sample"] for row in ordered[:size]} for size in SUBSET_SIZES]
    record["nesting_verified"] = all(nested[i] < nested[i + 1] for i in range(len(nested) - 1))
    record["finished_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (root / "sample_scaling.private.json").write_text(json.dumps(record, indent=2) + "\n")
    subprocess.run(["sync"], check=False)
    print("sample scaling measurement complete")


if __name__ == "__main__":
    main()
