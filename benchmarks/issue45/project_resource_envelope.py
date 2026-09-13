#!/usr/bin/env python3
"""Project the 327-sample Joint Genotyping resource envelope from committed evidence.

Every output keeps five things apart, because mixing them is how a projection
comes to look more certain than its inputs: what was observed, the empirical
trend fitted to those observations, the assumptions needed to carry the trend
to 327 samples, the projected range that follows, and the uncertainty left
over. The script decides nothing; the gate is a separate, reviewed judgement.

Sample count is not used alone as the predictor of GenotypeGVCFs memory. Along
the 13/26/51 prefix series cohort composition changes together with sample
count (seven highly divergent samples enter between 26 and 51), so an N-only
fit would attribute composition to N. Memory is modelled on genotype cells,
samples times output variants, which E0's eleven chromosomes test
independently at a fixed sample count. How many variants 327 samples produce
is unmeasured and is carried as explicit scenarios rather than one number.

Reads only committed evidence, so the projection can be re-derived by anyone.
"""

import argparse
import csv
import json
import math
from pathlib import Path

TARGET_SAMPLES = 327
BATCH_SIZE = 50
GIB = 2**30


def linear_fit(xs, ys):
    """Least squares y = a + b*x, with R^2 and the largest absolute residual."""
    n = len(xs)
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / sxx
    intercept = mean_y - slope * mean_x
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    total = sum((y - mean_y) ** 2 for y in ys)
    r2 = 1 - sum(r * r for r in residuals) / total if total else 1.0
    return {
        "intercept": intercept,
        "slope": slope,
        "r2": r2,
        "max_abs_residual": max(abs(r) for r in residuals),
    }


def power_fit(xs, ys):
    """Least squares in log space, y = c * x^k."""
    fit = linear_fit([math.log(x) for x in xs], [math.log(y) for y in ys])
    return {
        "coefficient": math.exp(fit["intercept"]),
        "exponent": fit["slope"],
        "r2_log": fit["r2"],
    }


def harmonic(n):
    """Watterson's a_n = sum_{i=1}^{n-1} 1/i."""
    return sum(1 / i for i in range(1, n))


def ratio_range(ns, ys, target=TARGET_SAMPLES):
    """Target/last-observed ratio under three N-trend forms; the range spans them."""
    last_n, last_y = ns[-1], ys[-1]
    linear = linear_fit(ns, ys)
    power = power_fit(ns, ys)
    upper_pair = last_y + (ys[-1] - ys[-2]) / (ns[-1] - ns[-2]) * (target - last_n)
    ratios = {
        "linear": (linear["intercept"] + linear["slope"] * target) / last_y,
        "power": power["coefficient"] * target ** power["exponent"] / last_y,
        "upper_pair_slope": upper_pair / last_y,
    }
    return {
        "linear_fit": linear,
        "power_fit": power,
        "ratio_to_last_observed": ratios,
        "ratio_min": min(ratios.values()),
        "ratio_max": max(ratios.values()),
    }


def variant_scenarios(v51, sites_first_in_14_26, sites_first_in_27_51):
    """Variant count for the measured window at 327 samples, one entry per assumption."""
    added = TARGET_SAMPLES - 51
    rate_typical = sites_first_in_14_26 / 13
    rate_divergent_mix = sites_first_in_27_51 / 25
    return {
        "S0_no_new_sites": {
            "variants": v51,
            "assumption": "the 276 additional samples add no variant site; a floor, since sites are practically never removed by adding samples",
        },
        "S1_neutral_harmonic": {
            "variants": v51 * harmonic(TARGET_SAMPLES) / harmonic(51),
            "assumption": "discovery follows Watterson's a_n from 51 to 327 samples, i.e. the 51 samples are representative of a panmictic 327; contradicted in part by the observed divergent subgroup",
        },
        "S2_typical_rate_constant": {
            "variants": v51 + added * rate_typical,
            "assumption": f"every added sample contributes {rate_typical:.0f} new sites, the rate observed for positions 14-26, without the decline discovery normally shows",
        },
        "S3_divergent_mix_rate_constant": {
            "variants": v51 + added * rate_divergent_mix,
            "assumption": f"every added sample contributes {rate_divergent_mix:.0f} new sites, the rate observed for positions 27-51 including seven divergent samples, without decline; a deliberately pessimistic bound",
        },
    }


def build_envelope(results, scaling, observed, fastq_327_bytes):
    experiments = results["experiments"]
    levels = sorted(scaling["levels"], key=lambda level: level["sample_count"])
    ns = [level["sample_count"] for level in levels]
    variants = [level["output_variant_count"] for level in levels]
    cells_m = [n * v / 1e6 for n, v in zip(ns, variants)]
    genotype_rss = [level["genotype_peak_rss_bytes"] / GIB for level in levels]

    # E0: one task per contig at a fixed 51 samples; chromosomes only, never capped.
    e0 = experiments["E0"]
    e0_counts = e0["output_per_contig_variant_counts"]
    # Counts are keyed in dictionary order and E0 runs one task per contig in that order.
    if e0["plan"] != "baseline_per_contig" or len(e0_counts) != len(e0["tasks"]):
        raise ValueError("E0 tasks do not pair one-to-one with contigs in dictionary order")
    chromosomes = [
        {
            "contig": contig,
            "bp": task["total_bp"],
            "variants": e0_counts[contig],
            "genotype_peak_rss_gib": task["genotype_peak_rss_gib"],
            "genotype_wall_seconds": task["genotype_wall_seconds"],
        }
        for contig, task in zip(e0_counts, e0["tasks"])
        if task["total_bp"] > 20_000_000
    ]
    e0_cells_fit = linear_fit(
        [51 * row["variants"] / 1e6 for row in chromosomes],
        [row["genotype_peak_rss_gib"] for row in chromosomes],
    )
    e0_length_fit = linear_fit(
        [row["bp"] / 1e6 for row in chromosomes],
        [row["genotype_peak_rss_gib"] for row in chromosomes],
    )
    e0_wall_length = linear_fit(
        [row["bp"] / 1e6 for row in chromosomes],
        [row["genotype_wall_seconds"] for row in chromosomes],
    )
    e0_wall_variants = linear_fit(
        [row["variants"] / 1e6 for row in chromosomes],
        [row["genotype_wall_seconds"] for row in chromosomes],
    )
    pair_slopes = [
        (genotype_rss[i + 1] - genotype_rss[i]) / (cells_m[i + 1] - cells_m[i])
        for i in range(len(levels) - 1)
    ]
    slopes = [e0_cells_fit["slope"], *pair_slopes]
    slope_low, slope_high = min(slopes), max(slopes)
    max_observed_cells_m = max(51 * row["variants"] / 1e6 for row in chromosomes)
    max_observed_cells_rss = max(row["genotype_peak_rss_gib"] for row in chromosomes)

    e2a = experiments["E2a"]
    genotype_wall = ratio_range(ns, [level["genotype_wall_seconds"] for level in levels])
    import_wall = ratio_range(ns, [level["import_wall_seconds"] for level in levels])
    workspace = ratio_range(ns, [level["workspace_bytes"] for level in levels])
    output_bytes = [level["output_vcf_bytes"] for level in levels]
    bytes_per_cell = [b / (c * 1e6) for b, c in zip(output_bytes, cells_m)]
    bytes_per_cell_slope = math.log(bytes_per_cell[-1] / bytes_per_cell[-2]) / math.log(
        cells_m[-1] / cells_m[-2]
    )
    one_batch_files = {
        level["workspace_file_count"] for level in levels if len(level["batches_completed"]) == 1
    }
    two_batch_files = {
        level["workspace_file_count"] for level in levels if len(level["batches_completed"]) == 2
    }
    if len(one_batch_files) != 1 or len(two_batch_files) != 1:
        raise ValueError("Workspace file count is not a single value per batching regime")
    files_per_batch = two_batch_files.pop() - next(iter(one_batch_files))
    files_fixed = next(iter(one_batch_files)) - files_per_batch
    batches_327 = math.ceil(TARGET_SAMPLES / BATCH_SIZE)

    composition = scaling["composition"]["sites_by_first_non_reference_carrier"]
    window = levels[-1]
    window_bp = 20_000_000
    scenarios = variant_scenarios(
        window["output_variant_count"],
        composition["positions_14_26"],
        composition["positions_27_51"],
    )
    window_cells_51 = cells_m[-1]
    window_rss_51 = genotype_rss[-1]
    window_density_51 = window["output_variant_count"] / (window_bp / 1e6)
    chromosome_variants = sum(row["variants"] for row in chromosomes)
    chromosome_bp = sum(row["bp"] for row in chromosomes)

    host = observed["host"]
    generation = observed["generation_51"]
    tier_gib = e2a["genotype_memory_gib"]
    concurrency = e2a["maximum_concurrent_tasks"]
    task_hours = {
        "import": [
            e2a["import_summary"]["wall_seconds_sum"] * r / 3600
            for r in (import_wall["ratio_min"], import_wall["ratio_max"])
        ],
        "genotype": [
            e2a["genotype_summary"]["wall_seconds_sum"] * r / 3600
            for r in (genotype_wall["ratio_min"], genotype_wall["ratio_max"])
        ],
    }
    scheduling_overhead = e2a["wall_seconds_including_validation"] / (
        (e2a["import_summary"]["wall_seconds_sum"] + e2a["genotype_summary"]["wall_seconds_sum"])
        / concurrency
    )
    fastq_ratio = fastq_327_bytes / generation["fastq_bytes"]
    sample_ratio = TARGET_SAMPLES / 51
    workspace_gb = [
        e2a["workspace_bytes"] * workspace["ratio_min"] / 1e9,
        e2a["workspace_bytes"] * max(workspace["ratio_max"], sample_ratio) / 1e9,
    ]

    projected_scenarios = {}
    for name, scenario in scenarios.items():
        growth = scenario["variants"] / window["output_variant_count"]
        cells_327 = TARGET_SAMPLES * scenario["variants"] / 1e6
        rss = [
            window_rss_51 + slope * (cells_327 - window_cells_51)
            for slope in (slope_low, slope_high)
        ]
        density_327 = window_density_51 * growth
        max_window_bp = max_observed_cells_m * 1e6 / (TARGET_SAMPLES * density_327) * 1e6
        chromosome_tasks = sum(math.ceil(row["bp"] / max_window_bp) for row in chromosomes)
        cells_ratio = TARGET_SAMPLES * growth / 51
        callset_gb = [
            observed["e2a_outputs"]["gathered_callset_bytes"] * cells_ratio * factor / 1e9
            for factor in (cells_ratio**bytes_per_cell_slope, 1.0)
        ]
        # The production plan is one task per contig; anchor at the longest one's own E0 reading.
        longest = max(chromosomes, key=lambda row: row["bp"])
        longest_rss = [
            longest["genotype_peak_rss_gib"]
            + slope * (TARGET_SAMPLES * growth - 51) * longest["variants"] / 1e6
            for slope in (slope_low, slope_high)
        ]
        storage_tb = [
            (fastq_327_bytes + generation["storage_consumed_bytes_net"] * ratio) / 1e12
            + workspace / 1e3
            + 2 * callset / 1e3
            for ratio, workspace, callset in zip(
                (sample_ratio, fastq_ratio), workspace_gb, callset_gb
            )
        ]
        projected_scenarios[name] = {
            "window_variants": round(scenario["variants"]),
            "variant_growth_vs_51": round(growth, 3),
            "window_genotype_cells_millions": round(cells_327, 1),
            "genotype_peak_rss_gib_20mb_window": [round(x, 1) for x in rss],
            "fits_current_genotype_tier": rss[1] <= tier_gib,
            "concurrent_tasks_within_launch_gate_at_projected_rss": math.floor(
                host["benchmark_launch_memory_gate_bytes"] / GIB / rss[1]
            ),
            "max_window_bp_within_observed_cells": round(max_window_bp),
            "chromosome_tasks_at_that_window": chromosome_tasks,
            "gathered_callset_gb": [round(x, 1) for x in callset_gb],
            "callset_gb_with_per_interval_copies": [round(2 * x, 1) for x in callset_gb],
            "production_per_contig_longest_genotype_rss_gib": [round(x, 1) for x in longest_rss],
            "full_run_retained_storage_tb": [round(x, 2) for x in storage_tb],
        }

    return {
        "schema_version": 1,
        "target_samples": TARGET_SAMPLES,
        "decides_nothing": "inputs to the reviewed 327 gate; see gate_327_decision.json",
        "observed": {
            "scaling_window": scaling["interval"],
            "scaling_levels": [
                {
                    "samples": n,
                    "output_variants": v,
                    "genotype_cells_millions": round(c, 3),
                    "genotype_peak_rss_gib": round(r, 3),
                    "genotype_wall_seconds": round(level["genotype_wall_seconds"], 1),
                    "import_wall_seconds": round(level["import_wall_seconds"], 1),
                    "import_peak_rss_gib": round(level["import_peak_rss_bytes"] / GIB, 3),
                    "workspace_bytes": level["workspace_bytes"],
                    "workspace_file_count": level["workspace_file_count"],
                    "output_vcf_bytes": level["output_vcf_bytes"],
                    "batches_completed": level["batches_completed"],
                }
                for n, v, c, r, level in zip(ns, variants, cells_m, genotype_rss, levels)
            ],
            "e0_chromosomes_at_51": chromosomes,
            "e2a_adopted_plan_at_51": {
                "interval_count": e2a["interval_count"],
                "import_wall_seconds_sum": e2a["import_summary"]["wall_seconds_sum"],
                "genotype_wall_seconds_sum": e2a["genotype_summary"]["wall_seconds_sum"],
                "elapsed_seconds_including_validation": e2a["wall_seconds_including_validation"],
                "import_peak_rss_gib_max": e2a["import_summary"]["peak_rss_gib_max"],
                "genotype_peak_rss_gib_max": e2a["genotype_summary"]["peak_rss_gib_max"],
                "workspace_bytes": e2a["workspace_bytes"],
                "workspace_file_count": e2a["workspace_file_count"],
                "gathered_callset_bytes": observed["e2a_outputs"]["gathered_callset_bytes"],
                "genotype_tier_gib": tier_gib,
                "concurrent_tasks": concurrency,
            },
            "chromosome_variant_density_per_mb_at_51": round(
                chromosome_variants / (chromosome_bp / 1e6)
            ),
            "scaling_window_variant_density_per_mb_at_51": round(window_density_51),
            "composition": scaling["composition"],
            "input_gvcf_51": observed["input_gvcf_51"],
            "fastq_bytes": {"selected_51": generation["fastq_bytes"], "all_327": fastq_327_bytes},
            "generation_51": generation,
            "host": host,
            "production_resource_contract": observed["production_resource_contract"],
        },
        "empirical_trends": {
            "genotype_rss_vs_cells": {
                "e0_eleven_chromosomes_fixed_51_samples": e0_cells_fit,
                "e0_same_points_vs_length_for_comparison": e0_length_fit,
                "scaling_pair_slopes_gib_per_million_cells": {
                    "13_to_26": pair_slopes[0],
                    "26_to_51": pair_slopes[1],
                },
                "slope_range_gib_per_million_cells": [slope_low, slope_high],
                "largest_observed_task_cells_millions": max_observed_cells_m,
                "largest_observed_task_rss_gib": max_observed_cells_rss,
                "reading": "three independent estimates, one at fixed N across contigs and two along N, agree within 0.16-0.19 GiB per million genotype cells",
            },
            "genotype_rss_vs_n_alone_not_used": {
                "power_fit": power_fit(ns, genotype_rss),
                "reason": "composition changes with N along the prefix series, so an N-only exponent is not a property of sample count",
            },
            "genotype_wall_vs_n": genotype_wall,
            "genotype_wall_at_fixed_n": {
                "vs_length_r2": e0_wall_length["r2"],
                "vs_variants_r2": e0_wall_variants["r2"],
                "reading": "wall time follows interval length more closely than variant count, so it is less composition-sensitive than memory",
            },
            "import_wall_vs_n": import_wall,
            "import_peak_rss_gib": {
                "observed": [round(level["import_peak_rss_bytes"] / GIB, 3) for level in levels],
                "reading": "no trend with N; readers are bounded by the batch size",
            },
            "workspace_bytes_vs_n": workspace,
            "workspace_files_per_array": {
                "fixed": files_fixed,
                "per_batch": files_per_batch,
                "reading": "47 files with one batch and 87 with two",
            },
            "output_bytes_per_cell": {
                "observed": [round(x, 2) for x in bytes_per_cell],
                "log_log_slope_upper_pair": bytes_per_cell_slope,
            },
        },
        "assumptions": {
            "A1_cells_model_extrapolates": f"GenotypeGVCFs peak RSS stays linear in genotype cells beyond the largest observed task ({max_observed_cells_m:.0f} M cells)",
            "A2_variant_scenarios": {name: s["assumption"] for name, s in scenarios.items()},
            "A3_window_is_representative": "the scaling window's variant growth applies genome-wide; at 51 samples it is denser than the chromosome average, which makes window-level memory conservative",
            "A4_seven_batches": f"{batches_327} import batches behave like the observed one- and two-batch regimes, per sample and per batch",
            "A5_same_host_and_io": f"the same host, storage and {concurrency}-task concurrency as the suite",
            "A6_upstream_scales_with_fastq_bytes": "generation storage and CPU scale with FASTQ bytes (ratio from the published 327 sizes) or with sample count; both are shown",
            "A7_retained_intermediates": "production has no cleanup configured, so a full run retains intermediates the way the 51-sample generation did",
        },
        "projected": {
            "per_variant_scenario": projected_scenarios,
            "task_hours": {
                "import": [round(x, 1) for x in task_hours["import"]],
                "genotype": [round(x, 1) for x in task_hours["genotype"]],
                "elapsed_hours_at_3_tasks": [
                    round((i + g) / 3 * scheduling_overhead, 1)
                    for i, g in zip(task_hours["import"], task_hours["genotype"])
                ],
                "elapsed_hours_at_2_tasks": [
                    round((i + g) / 2 * scheduling_overhead, 1)
                    for i, g in zip(task_hours["import"], task_hours["genotype"])
                ],
                "elapsed_hours_at_1_task": [
                    round((i + g) * scheduling_overhead, 1)
                    for i, g in zip(task_hours["import"], task_hours["genotype"])
                ],
                "scheduling_overhead_factor_observed_e2a": round(scheduling_overhead, 3),
            },
            "workspace": {
                "bytes_gb": [round(x, 1) for x in workspace_gb],
                "files_per_array": files_fixed + files_per_batch * batches_327,
                "files_e2a_plan": e2a["interval_count"]
                * (files_fixed + files_per_batch * batches_327),
            },
            "import_batches_per_array": batches_327,
            "input_gvcf_gb": [
                round(observed["input_gvcf_51"]["bytes_total"] * sample_ratio / 1e9, 1),
                round(observed["input_gvcf_51"]["bytes_total"] * fastq_ratio / 1e9, 1),
            ],
            "generation": {
                "retained_storage_tb": [
                    round(generation["storage_consumed_bytes_net"] * sample_ratio / 1e12, 2),
                    round(generation["storage_consumed_bytes_net"] * fastq_ratio / 1e12, 2),
                ],
                "cpu_hours": [
                    round(generation["cpu_hours_total"] * sample_ratio),
                    round(generation["cpu_hours_total"] * fastq_ratio),
                ],
                "elapsed_hours_same_allocation": [
                    round(generation["elapsed_seconds"] * sample_ratio / 3600),
                    round(generation["elapsed_seconds"] * fastq_ratio / 3600),
                ],
            },
            "fastq_tb": round(fastq_327_bytes / 1e12, 2),
            "data_volume_tb": {
                "total": round(host["data_volume_total_bytes"] / 1e12, 2),
                "free_now": round(host["data_volume_free_bytes_after_scaling"] / 1e12, 2),
            },
        },
        "uncertainty": [
            "Variant count at 327 samples is unmeasured; S0-S3 span a factor of about eight and dominate the memory and callset ranges.",
            f"A1 extrapolates the cells model beyond its largest observation ({max_observed_cells_m:.0f} M cells); even S0 needs about three times that for a 20 Mb window.",
            "Peak RSS is a resident-set reading under a 26 GiB heap, not a measured live set; the E0 OOM shows the live set does track it near a ceiling.",
            "The 51-sample cohort is an accession-order selection, not a random sample of the 327; the divergent subgroup's share of the 327 is unknown.",
            "Import at seven batches is extrapolated from one- and two-batch regimes.",
            "Wall-time ranges come from one window and the three N-trend forms; composition could add to genotype time, though less than to memory.",
            "Upstream generation figures are one 51-sample run on this host with every intermediate retained; they are context for a full run, not a Joint Genotyping measurement.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-results", type=Path, required=True)
    parser.add_argument("--sample-scaling", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--candidate-cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.candidate_cohort.open() as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["library_strategy"] == "WGS" and row["library_layout"] == "PAIRED"
        ]
    if len(rows) != TARGET_SAMPLES:
        raise ValueError("Candidate cohort no longer has 327 WGS paired-end runs")
    fastq_327 = sum(int(size) for row in rows for size in row["fastq_bytes"].split(";"))
    envelope = build_envelope(
        json.loads(args.benchmark_results.read_text()),
        json.loads(args.sample_scaling.read_text()),
        json.loads(args.observed.read_text()),
        fastq_327,
    )
    args.output.write_text(json.dumps(envelope, indent=2) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
