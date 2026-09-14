"""Execute the fixed one-factor Issue45 experiments after independent validation.

This module does not release RAM or grant execution permission. The controller
must own the resource guard and service restoration for the entire call.
"""

import concurrent.futures
import csv
import datetime
import io
import json
import sys
import time
from pathlib import Path

from benchmark_tools import (
    GATK,
    ToolRunner,
    batches_from_log,
    compare_callsets,
    filesystem_metrics,
    joint_integrity,
)
from stage_generation import sha256
from validate_generated_cohort import bcftools_records, read_header, verify_bgzf

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from prepare_joint_genotyping_benchmark import (  # noqa: E402
    INTERVAL_PLAN_HEADER,
    build_interval_plan_rows,
    read_reference_fai,
)

# Resource policy, fixed identically for every experiment so that only the
# declared one-factor change differs between a pair being compared.
#
# Issue #45 first real run: GenotypeGVCFs was OOM-killed (exit 247,
# OOMKilled=true) on the largest contig at 51 samples with a 15 GiB Java heap
# inside a 16 GiB container, after its progress meter collapsed from ~1.2M to
# ~3 records/minute -- a GC death spiral, so the live set genuinely approached
# the heap ceiling rather than the container merely clipping native overhead.
# Per this repository's standing methodology (nextflow.config, Issues #30/#33)
# a reading at its own cgroup ceiling is not a true peak, so GenotypeGVCFs was
# re-measured alone at a deliberately generous ceiling and given its own tier
# below. GenomicsDBImport keeps 16 GiB: its measured peak was 1.67 GiB, so that
# allocation is already ample and is left unchanged. This is a resource
# allocation change only -- no scientific parameter, threshold, ploidy,
# reference or interval semantics is affected.
MAXIMUM_CONCURRENT_TASKS = 3
GENOMICSDB_MEMORY_GIB = 16
GENOMICSDB_JAVA_HEAP_MIB = 13107
GENOTYPE_MEMORY_GIB = 32
# GenomicsDB's native TileDB buffers, plus JVM metaspace/stacks/code cache, all
# sit outside the Java heap, so the heap is derived this far below the container
# limit rather than being written as its own independent number.
GENOTYPE_NATIVE_RESERVE_GIB = 6
GATHER_MEMORY_GIB = 8
GATHER_NATIVE_RESERVE_GIB = 1
REBLOCK_MEMORY_GIB = 16
REBLOCK_NATIVE_RESERVE_GIB = 3
COMPLETED = "COMPLETED_AWAITING_COMPARATIVE_REVIEW"

EXPERIMENTS = (
    ("E0", None, "baseline_per_contig", False, False, False),
    ("E1", "E0", "baseline_per_contig", True, False, False),
    ("E2a", "E1", "candidate_split_only", True, False, False),
    ("E2b", "E2a", "candidate_split_group", True, False, False),
    ("E3", "E2b", "candidate_split_group", True, True, False),
    ("E4", "E2b", "candidate_split_group", True, False, True),
)


def validate_tiling(groups, contigs):
    """Independently verify exact coverage/order of the serialized execution plan."""
    lengths = dict(contigs)
    expected_index, expected_start = 0, 1
    for group in groups:
        for interval in group["intervals"]:
            if ":" in interval:
                name, coordinates = interval.rsplit(":", 1)
                start, end = (int(value) for value in coordinates.split("-"))
            else:
                name, start, end = interval, 1, lengths[interval]
            if expected_index >= len(contigs) or name != contigs[expected_index][0]:
                raise ValueError("Interval plan reordered or repeated a reference contig")
            if start != expected_start or not start <= end <= lengths[name]:
                raise ValueError("Interval plan has a gap, overlap, or invalid endpoint")
            if end == lengths[name]:
                expected_index += 1
                expected_start = 1
            else:
                expected_start = end + 1
    if expected_index != len(contigs) or expected_start != 1:
        raise ValueError("Interval plan does not cover the complete reference")


def parallel_map(function, values, workers=MAXIMUM_CONCURRENT_TASKS):
    # Keep only three intervals/samples active. A failed future prevents queued
    # tasks from starting; active tasks retain their own measurements.
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    futures = [pool.submit(function, value) for value in values]
    positions = {future: index for index, future in enumerate(futures)}
    results = [None] * len(futures)
    try:
        for future in concurrent.futures.as_completed(futures):
            results[positions[future]] = future.result()
        return results
    finally:
        pool.shutdown(wait=True, cancel_futures=True)


def process_summary(records):
    return {
        "task_count": len(records),
        "wall_seconds_sum": sum(row["wall_seconds"] for row in records),
        "peak_rss_bytes_max": max((row["peak_rss_bytes"] for row in records), default=0),
        "user_cpu_seconds_sum": sum(row["user_cpu_seconds"] for row in records),
        "system_cpu_seconds_sum": sum(row["system_cpu_seconds"] for row in records),
        "retry_count": sum(row["retry_count"] for row in records),
    }


# The suite deliberately continues past a Reblock failure, because E4's
# consolidate comparison is against E2b and does not depend on E3. A resumed run
# has to reach the same conclusion from the retained record instead of treating
# that already-accepted outcome as a reason to refuse to continue.
CONTINUES_AFTER_FAILURE = ("E3",)


def completed_experiment(root, experiment):
    """Reuse an already-finished experiment, but never overwrite retained failure evidence."""
    directory = root / experiment
    record = directory / "experiment_result.private.json"
    if not directory.exists():
        return None
    if not record.exists():
        raise ValueError(
            f"{experiment} directory exists without a result; retain it and use a new run"
        )
    result = json.loads(record.read_text())
    if result.get("status") != COMPLETED and experiment not in CONTINUES_AFTER_FAILURE:
        raise ValueError(
            f"{experiment} holds retained {result.get('status')} evidence; "
            "diagnose it and start a new run directory rather than overwriting it"
        )
    return result


def run_experiment(
    runner, directory, configuration, groups, samples, reference_name, contigs, batch_size=50
):
    directory.mkdir(exist_ok=False)
    started = time.monotonic()
    config = {
        **configuration,
        "sample_count": len(samples),
        "batch_size": batch_size,
        "interval_count": len(groups),
        "intervals": groups,
        "gatk_container": GATK,
        "maximum_concurrent_tasks": MAXIMUM_CONCURRENT_TASKS,
        "genomicsdb_cpus": 8,
        "genomicsdb_memory_bytes": GENOMICSDB_MEMORY_GIB * 1024**3,
        "genomicsdb_java_heap_mib": GENOMICSDB_JAVA_HEAP_MIB,
        "genotype_memory_bytes": runner.effective_memory_gib(GENOTYPE_MEMORY_GIB) * 1024**3,
        "genotype_java_heap_gib": runner.java_heap_gib(
            GENOTYPE_MEMORY_GIB, GENOTYPE_NATIVE_RESERVE_GIB
        ),
        "sample_ploidy": 2,
        "reader_threads_requested": 8,
        "only_output_calls_starting_in_intervals": True,
        "reblock_policy": {
            "gq_bands": [20, 100],
            "keep_all_alts": True,
            "floor_blocks": False,
            "drop_low_quals": False,
        },
        "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    (directory / "experiment_config.json").write_text(json.dumps(config, indent=2) + "\n")
    source_paths = {row["sample"]: "/input/" + row["gvcf_name"] for row in samples}
    result = {"configuration": config, "status": "RUNNING", "tasks": []}
    try:
        if config["reblock"]:
            reblocked = directory / "reblocked"
            reblocked.mkdir()

            def reblock(sample):
                sid = sample["sample"]
                target = reblocked / f"{sid}.g.vcf.gz"
                tool_dir = directory / "reblock_tasks" / sid
                container_target = "/bench/" + str(target.relative_to(runner.root))
                metric = runner.run(
                    tool_dir,
                    "ReblockGVCF",
                    [
                        "gatk",
                        "--java-options",
                        f"-Xmx{runner.java_heap_gib(REBLOCK_MEMORY_GIB, REBLOCK_NATIVE_RESERVE_GIB)}g",
                        "ReblockGVCF",
                        "--reference",
                        "/reference/" + reference_name,
                        "--variant",
                        source_paths[sid],
                        "--output",
                        container_target,
                        "--gvcf-gq-bands",
                        "20",
                        "--gvcf-gq-bands",
                        "100",
                        "--keep-all-alts",
                        "true",
                        "--floor-blocks",
                        "false",
                        "--drop-low-quals",
                        "false",
                        "--create-output-variant-index",
                        "true",
                    ],
                    cpus=4,
                    memory_gib=REBLOCK_MEMORY_GIB,
                )
                return {
                    "sample": sid,
                    "metrics": metric,
                    "input_bytes": sample["gvcf_bytes"],
                    "output_bytes": target.stat().st_size,
                    "gvcf_sha256": sha256(target),
                    "index_sha256": sha256(Path(str(target) + ".tbi")),
                }

            derived = parallel_map(reblock, samples)
            validation_dir = directory / "reblock_validation"
            validation_dir.mkdir()
            first_header = None
            for row in derived:
                target = reblocked / f"{row['sample']}.g.vcf.gz"
                header = read_header(target, row["sample"], contigs)
                fingerprint = header["shared_header_definitions_sha256"]
                if first_header is not None and fingerprint != first_header:
                    raise ValueError("Reblocked sample headers are heterogeneous")
                first_header = fingerprint
                verify_bgzf(target)
                sequential = bcftools_records(target, contigs, False, validation_dir)
                indexed = bcftools_records(target, contigs, True, validation_dir)
                if sequential != indexed:
                    raise ValueError("Reblocked gVCF/index traversal mismatch")
                row.update(sequential)
                source_paths[row["sample"]] = "/bench/" + str(target.relative_to(runner.root))
            result["reblock_samples"] = derived
            result["reblock"] = {
                **process_summary([row["metrics"] for row in derived]),
                "input_bytes": sum(row["input_bytes"] for row in derived),
                "output_bytes": sum(row["output_bytes"] for row in derived),
            }
        sample_map = directory / "sample_name_map.tsv"
        sample_map.write_text(
            "".join(f"{sid}\t{path}\t{path}.tbi\n" for sid, path in source_paths.items())
        )

        def interval_task(group):
            task = directory / "intervals" / group["id"]
            intervals = [
                argument
                for interval in group["intervals"]
                for argument in ("--intervals", interval)
            ]
            if config["sample_name_map"]:
                inputs = [
                    "--sample-name-map",
                    "/bench/" + str(sample_map.relative_to(runner.root)),
                    "--validate-sample-name-map",
                    "true",
                ]
            else:
                inputs = [
                    argument for path in source_paths.values() for argument in ("--variant", path)
                ]
            imported = runner.run(
                task,
                "GenomicsDBImport",
                [
                    "gatk",
                    "--java-options",
                    f"-Xmx{GENOMICSDB_JAVA_HEAP_MIB}m",
                    "GenomicsDBImport",
                    *inputs,
                    "--genomicsdb-workspace-path",
                    "workspace",
                    *intervals,
                    "--reader-threads",
                    "8",
                    "--batch-size",
                    str(batch_size),
                    "--consolidate",
                    str(config["consolidate"]).lower(),
                    "--tmp-dir",
                    ".",
                ],
            )
            batch = batches_from_log(
                (task / "GenomicsDBImport.log").read_text(), len(samples), batch_size
            )
            workspace = filesystem_metrics(task / "workspace")
            genotyped = runner.run(
                task,
                "GenotypeGVCFs",
                [
                    "gatk",
                    "--java-options",
                    f"-Xmx{runner.java_heap_gib(GENOTYPE_MEMORY_GIB, GENOTYPE_NATIVE_RESERVE_GIB)}g",
                    "GenotypeGVCFs",
                    "--reference",
                    "/reference/" + reference_name,
                    "--variant",
                    "gendb://workspace",
                    *intervals,
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
            return {
                "group": group,
                "import": imported,
                "genotype": genotyped,
                "batching": batch,
                "workspace": workspace,
            }

        result["tasks"] = parallel_map(interval_task, groups)
        gathered_dir = directory / "gather"
        arguments = [
            argument
            for group in groups
            for argument in (
                "--INPUT",
                "/bench/"
                + str(
                    (directory / "intervals" / group["id"] / "joint.vcf.gz").relative_to(
                        runner.root
                    )
                ),
            )
        ]
        result["gather"] = runner.run(
            gathered_dir,
            "GatherVcfs",
            [
                "gatk",
                "--java-options",
                f"-Xmx{runner.java_heap_gib(GATHER_MEMORY_GIB, GATHER_NATIVE_RESERVE_GIB)}g",
                "GatherVcfs",
                *arguments,
                "--OUTPUT",
                "cohort.raw.vcf.gz",
            ],
            cpus=4,
            memory_gib=GATHER_MEMORY_GIB,
        )
        result["gather_index"] = runner.run(
            gathered_dir,
            "IndexFeatureFile",
            [
                "gatk",
                "--java-options",
                f"-Xmx{runner.java_heap_gib(GATHER_MEMORY_GIB, GATHER_NATIVE_RESERVE_GIB)}g",
                "IndexFeatureFile",
                "--input",
                "cohort.raw.vcf.gz",
            ],
            cpus=4,
            memory_gib=GATHER_MEMORY_GIB,
        )
        result["integrity"] = joint_integrity(
            gathered_dir / "cohort.raw.vcf.gz", [row["sample"] for row in samples], contigs
        )
        result["import"] = process_summary([row["import"] for row in result["tasks"]])
        result["genotype"] = process_summary([row["genotype"] for row in result["tasks"]])
        result["workspace"] = {
            key: sum(row["workspace"][key] for row in result["tasks"])
            for key in ("bytes", "file_count")
        }
        result["status"] = COMPLETED
    except BaseException as error:
        result.update(status="FAILED", error_type=type(error).__name__, error_private=str(error))
        raise
    finally:
        result["wall_seconds_including_validation"] = time.monotonic() - started
        result["finished_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (directory / "experiment_result.private.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
    return result


def execute_suite(root, validated, stop_event):
    evidence = validated["evidence"]
    if not evidence.get("benchmark_ready") or not evidence.get("lineage_verified"):
        raise ValueError("Independent cohort validation has not opened the execution gate")
    samples = evidence["samples"]
    if (
        evidence["sample_count"] != 51
        or len(samples) != 51
        or len({r["sample"] for r in samples}) != 51
    ):
        raise ValueError("Exactly 51 unique validated samples are required")
    reference = validated["reference_paths"]
    for role, path in reference.items():
        if sha256(Path(path)) != evidence["reference"][role]["sha256"]:
            raise ValueError("Validated reference changed")
    input_dir = Path(validated["gvcf_directory"])
    for sample in samples:
        gvcf = input_dir / sample["gvcf_name"]
        if (
            sha256(gvcf) != sample["gvcf_sha256"]
            or sha256(Path(str(gvcf) + ".tbi")) != sample["gvcf_index_sha256"]
        ):
            raise ValueError("Validated gVCF/index changed before benchmark")
    ordered = read_reference_fai(Path(reference["fai"]))
    contigs = [(row.name, row.length) for row in ordered]
    rows = build_interval_plan_rows(ordered, 20_000_000, 1_000_000)
    plan_path = root / "interval_plan.tsv"
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(INTERVAL_PLAN_HEADER)
    writer.writerows(rows)
    if plan_path.exists():
        # A resumed run must reproduce the identical plan the earlier experiments used.
        if plan_path.read_text() != buffer.getvalue():
            raise ValueError("Interval plan changed since this run directory was started")
    else:
        plan_path.write_text(buffer.getvalue())
    plans = {}
    for row in rows:
        plans.setdefault(row[0], []).append(
            {
                "id": f"interval_{int(row[2]):04d}",
                "intervals": json.loads(row[3]),
                "total_bp": int(row[4]),
            }
        )
    for groups in plans.values():
        validate_tiling(groups, contigs)
    runner = ToolRunner(root, input_dir, Path(reference["fasta"]).parent, stop_event)
    results = {}
    for experiment, compare_to, plan, sample_map, reblock, consolidate in EXPERIMENTS:
        earlier = completed_experiment(root, experiment)
        if earlier is not None:
            results[experiment] = earlier
            continue
        if stop_event.is_set():
            raise RuntimeError("Resource guard stopped the suite")
        configuration = {
            "experiment_id": experiment,
            "compare_to": compare_to,
            "plan": plan,
            "sample_name_map": sample_map,
            "reblock": reblock,
            "consolidate": consolidate,
        }
        try:
            result = run_experiment(
                runner,
                root / experiment,
                configuration,
                plans[plan],
                samples,
                Path(reference["fasta"]).name,
                contigs,
            )
        except Exception:
            if experiment not in CONTINUES_AFTER_FAILURE:
                raise
            # A real Reblock failure does not block the independent consolidate comparison.
            result = json.loads((root / experiment / "experiment_result.private.json").read_text())
        if compare_to and "integrity" in result:
            baseline = results[compare_to]["integrity"]
            result["comparison"] = {
                key: result["integrity"][key] == baseline[key]
                for key in (
                    "sample_order",
                    "ordered_contigs",
                    "variant_count",
                    "record_sha256",
                    "variant_gt_sha256",
                    "accounting_sha256",
                    "shared_header_sha256",
                )
            }
            result["record_difference_audit"] = compare_callsets(
                root / compare_to / "gather/cohort.raw.vcf.gz",
                root / experiment / "gather/cohort.raw.vcf.gz",
                contigs,
            )
        results[experiment] = result
        (root / "suite_results.private.json").write_text(json.dumps(results, indent=2) + "\n")
    for sample in samples:
        path = input_dir / sample["gvcf_name"]
        if (
            sha256(path) != sample["gvcf_sha256"]
            or sha256(Path(str(path) + ".tbi")) != sample["gvcf_index_sha256"]
        ):
            raise ValueError("Source inputs changed during benchmark")
    return results
