#!/usr/bin/env python3
"""Turn a private Issue45 benchmark run into the sanitized evidence committed to Git.

Absolute paths, usernames and host directory layout never reach the public
record. Sample identities are public accessions and stay as they are. Nothing
here re-derives a measurement: it selects, renames and aggregates what the run
already recorded, so a disagreement with the private record is a bug here.
"""

import argparse
import json
from pathlib import Path

EXPERIMENT_ORDER = ("E0", "E1", "E2a", "E2b", "E3", "E4")
COMPARISON_KEYS = (
    "sample_order",
    "ordered_contigs",
    "variant_count",
    "record_sha256",
    "variant_gt_sha256",
    "accounting_sha256",
    "shared_header_sha256",
)


def gib(value):
    return round(value / 1024**3, 3)


def process_block(summary):
    return {
        "task_count": summary["task_count"],
        "wall_seconds_sum": round(summary["wall_seconds_sum"], 1),
        "peak_rss_gib_max": gib(summary["peak_rss_bytes_max"]),
        "user_cpu_seconds_sum": round(summary["user_cpu_seconds_sum"], 1),
        "system_cpu_seconds_sum": round(summary["system_cpu_seconds_sum"], 1),
        "retry_count": summary["retry_count"],
    }


def task_block(task):
    group = task["group"]
    return {
        "interval_id": group["id"],
        "interval_count": len(group["intervals"]),
        "total_bp": group["total_bp"],
        "import_wall_seconds": round(task["import"]["wall_seconds"], 1),
        "import_peak_rss_gib": gib(task["import"]["peak_rss_bytes"]),
        "import_attempt": task["import"]["attempt"],
        "genotype_wall_seconds": round(task["genotype"]["wall_seconds"], 1),
        "genotype_peak_rss_gib": gib(task["genotype"]["peak_rss_bytes"]),
        "genotype_attempt": task["genotype"]["attempt"],
        "workspace_bytes": task["workspace"]["bytes"],
        "workspace_file_count": task["workspace"]["file_count"],
        "batches_observed": task["batching"]["observed"],
        "batches_completed": task["batching"]["completed"],
        "reader_initialization_fell_back_to_serial": task["batching"][
            "reader_initialization_fell_back_to_serial"
        ],
    }


def experiment_block(experiment, result):
    config = result["configuration"]
    block = {
        "experiment_id": experiment,
        "compare_to": config["compare_to"],
        "status": result["status"],
        "plan": config["plan"],
        "sample_count": config["sample_count"],
        "batch_size": config["batch_size"],
        "interval_count": config["interval_count"],
        "sample_name_map": config["sample_name_map"],
        "reblock": config["reblock"],
        "consolidate": config["consolidate"],
        "genomicsdb_cpus": config["genomicsdb_cpus"],
        "genomicsdb_memory_gib": gib(config["genomicsdb_memory_bytes"]),
        "genotype_memory_gib": gib(config["genotype_memory_bytes"]),
        "genotype_java_heap_gib": config["genotype_java_heap_gib"],
        "maximum_concurrent_tasks": config["maximum_concurrent_tasks"],
        "only_output_calls_starting_in_intervals": config[
            "only_output_calls_starting_in_intervals"
        ],
        "wall_seconds_including_validation": round(result["wall_seconds_including_validation"], 1),
    }
    if result["status"] != "COMPLETED_AWAITING_COMPARATIVE_REVIEW":
        block["error_type"] = result.get("error_type")
        return block
    integrity = result["integrity"]
    block.update(
        import_summary=process_block(result["import"]),
        genotype_summary=process_block(result["genotype"]),
        gather_wall_seconds=round(result["gather"]["wall_seconds"], 1),
        gather_peak_rss_gib=gib(result["gather"]["peak_rss_bytes"]),
        gather_index_wall_seconds=round(result["gather_index"]["wall_seconds"], 1),
        gather_index_peak_rss_gib=gib(result["gather_index"]["peak_rss_bytes"]),
        workspace_bytes=result["workspace"]["bytes"],
        workspace_file_count=result["workspace"]["file_count"],
        output_sample_count=integrity["sample_count"],
        output_sample_order_sha256=integrity["sample_order_sha256"],
        output_variant_count=integrity["variant_count"],
        output_called_alleles=integrity["called_alleles"],
        output_record_sha256=integrity["record_sha256"],
        output_variant_gt_sha256=integrity["variant_gt_sha256"],
        output_accounting_sha256=integrity["accounting_sha256"],
        output_shared_header_sha256=integrity["shared_header_sha256"],
        output_per_contig_variant_counts=integrity["per_contig_variant_counts"],
        sample_order_valid=integrity["sample_order_valid"],
        contig_order_valid=integrity["contig_order_valid"],
        allele_accounting_valid=integrity["allele_accounting_valid"],
        tasks=[task_block(task) for task in result["tasks"]],
    )
    if "reblock" in result:
        block["reblock_summary"] = {
            **process_block(result["reblock"]),
            "input_bytes": result["reblock"]["input_bytes"],
            "output_bytes": result["reblock"]["output_bytes"],
            "size_ratio": round(
                result["reblock"]["output_bytes"] / result["reblock"]["input_bytes"], 4
            ),
        }
    if "comparison" in result:
        block["comparison_to_baseline"] = {
            key: result["comparison"][key] for key in COMPARISON_KEYS
        }
    if "record_difference_audit" in result:
        block["record_difference_audit"] = result["record_difference_audit"]
    return block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--validated-cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    lineage = json.loads((args.run_dir / "execution_lineage.json").read_text())
    cohort = json.loads(args.validated_cohort.read_text())["evidence"]
    host = [
        line.split("\t")
        for line in (args.run_dir / "host.tsv").read_text().splitlines()[1:]
        if line.strip()
    ]
    memory_values = [int(row[1]) for row in host]
    swap_values = [int(row[2]) for row in host]
    storage_values = [int(row[3]) for row in host]

    experiments = {}
    for experiment in EXPERIMENT_ORDER:
        record = args.run_dir / experiment / "experiment_result.private.json"
        if record.exists():
            experiments[experiment] = experiment_block(experiment, json.loads(record.read_text()))

    evidence = {
        "schema_version": 1,
        "run_id": args.run_dir.name,
        "cohort": {
            "cohort_id": cohort["cohort_id"],
            "sample_count": cohort["sample_count"],
            "unique_biosamples": cohort["unique_biosamples"],
            "bioproject": sorted({row["bioproject"] for row in cohort["samples"]}),
            "reference_accession": cohort["reference_accession"],
            "production_sha": cohort["production_sha"],
            "nextflow_version": cohort["nextflow_version"],
            "gatk_version": cohort["gatk_version"],
            "gatk_container": cohort["gatk_container"],
            "sample_ploidy": cohort["sample_ploidy"],
            "lineage_verified": cohort["lineage_verified"],
            "benchmark_ready": cohort["benchmark_ready"],
            "input_gvcf_bytes": sum(row["gvcf_bytes"] for row in cohort["samples"]),
            "input_gvcf_records": sum(row["record_count"] for row in cohort["samples"]),
        },
        "execution_lineage": lineage,
        "host_during_run": {
            "monitor_samples": len(host),
            "mem_available_min_gib": gib(min(memory_values)) if memory_values else None,
            "swap_start_gib": gib(swap_values[0]) if swap_values else None,
            "swap_peak_gib": gib(max(swap_values)) if swap_values else None,
            "swap_end_gib": gib(swap_values[-1]) if swap_values else None,
            "swap_delta_gib": gib(swap_values[-1] - swap_values[0]) if swap_values else None,
            "storage_free_min_bytes": min(storage_values) if storage_values else None,
        },
        "experiments": experiments,
    }
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=False) + "\n")
    print(f"wrote {args.output} for {len(experiments)} experiments")


if __name__ == "__main__":
    main()
