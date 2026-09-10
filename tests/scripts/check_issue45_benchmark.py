#!/usr/bin/env python3
"""Exercise real pinned benchmark commands on two synthetic samples only."""

import argparse
import json
import shutil
import sys
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks/issue45"))
from benchmark_tools import ToolRunner, compare_callsets  # noqa: E402
from execute_experiments import parallel_map, run_experiment, validate_tiling  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generation-run", type=Path, required=True)
    parser.add_argument("--reference-directory", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs, reference, gvcfs = (root / name for name in ("input", "reference", "gvcfs"))
    for directory in (inputs, reference, gvcfs):
        directory.mkdir()
    source = args.reference_directory.resolve()
    for name in ("synthetic.fa", "synthetic.fa.fai"):
        shutil.copyfile(source / name, reference / name)
    shutil.copyfile(source / "prebuilt_dict/synthetic.dict", reference / "synthetic.dict")
    ids = ["sample_a", "sample_b"]
    for sid in ids:
        for suffix in (".markdup.bam", ".markdup.bam.bai"):
            candidates = {
                p.resolve() for p in (args.generation_run / "generation-work").rglob(sid + suffix)
            }
            if len(candidates) != 1:
                raise ValueError("Expected one synthetic BAM/index identity")
            shutil.copyfile(candidates.pop(), inputs / (sid + suffix))
    stop = threading.Event()
    producer = ToolRunner(root, inputs, reference, stop)

    def generate(sid):
        producer.run(
            root / "generation" / sid,
            "HaplotypeCaller",
            [
                "gatk",
                "--java-options",
                "-Xmx3g",
                "HaplotypeCaller",
                "--reference",
                "/reference/synthetic.fa",
                "--input",
                f"/input/{sid}.markdup.bam",
                "--output",
                f"/bench/gvcfs/{sid}.g.vcf.gz",
                "--emit-ref-confidence",
                "GVCF",
                "--sample-ploidy",
                "2",
                "--native-pair-hmm-threads",
                "4",
                "--create-output-variant-index",
                "true",
            ],
            cpus=4,
            memory_gib=4,
        )

    parallel_map(generate, ids, workers=2)
    samples = [
        {
            "sample": sid,
            "gvcf_name": sid + ".g.vcf.gz",
            "gvcf_bytes": (gvcfs / (sid + ".g.vcf.gz")).stat().st_size,
        }
        for sid in ids
    ]
    contigs = [("chrSynthetic1", 5000), ("chrSynthetic2", 5000)]
    baseline = [
        {"id": f"interval_{i:04d}", "intervals": [name], "total_bp": length}
        for i, (name, length) in enumerate(contigs)
    ]
    split = [
        {
            "id": f"interval_{i:04d}_{j}",
            "intervals": [f"{name}:{start}-{end}"],
            "total_bp": end - start + 1,
        }
        for i, (name, _) in enumerate(contigs)
        for j, (start, end) in enumerate(((1, 2500), (2501, 5000)))
    ]
    grouped = [
        {"id": "interval_0000", "intervals": [name for name, _ in contigs], "total_bp": 10000}
    ]
    for groups in (baseline, split, grouped):
        validate_tiling(groups, contigs)
    runner = ToolRunner(root, gvcfs, reference, stop)
    results = {}
    # The two equally sized synthetic contigs cannot represent both long
    # chromosomes and small scaffolds. Test splitting and grouping separately
    # against E1 here; the real suite retains E2a -> E2b dependency and one factor.
    cases = (
        ("baseline", baseline, False, False, False),
        ("sample_map", baseline, True, False, False),
        ("split", split, True, False, False),
        ("group", grouped, True, False, False),
        ("reblock", grouped, True, True, False),
        ("consolidate", grouped, True, False, True),
    )
    for name, groups, sample_map, reblock, consolidate in cases:
        config = {
            "experiment_id": "SYNTHETIC_" + name,
            "compare_to": None,
            "sample_name_map": sample_map,
            "reblock": reblock,
            "consolidate": consolidate,
        }
        results[name] = run_experiment(
            runner, root / name, config, groups, samples, "synthetic.fa", contigs, batch_size=1
        )
        if name not in ("baseline", "reblock"):
            for key in (
                "sample_order",
                "ordered_contigs",
                "variant_gt_sha256",
                "accounting_sha256",
            ):
                if results[name]["integrity"][key] != results["baseline"]["integrity"][key]:
                    raise ValueError(f"Synthetic {name} changes {key}")
            audit = compare_callsets(
                root / "baseline/gather/cohort.raw.vcf.gz",
                root / name / "gather/cohort.raw.vcf.gz",
                contigs,
            )
            if not audit[
                "all_differences_explained_by_high_qd_jitter_with_unchanged_current_filter"
            ]:
                raise ValueError(f"Unexplained synthetic {name} record difference: {audit}")
            results[name]["record_difference_audit"] = audit
    (root / "synthetic_results.private.json").write_text(json.dumps(results, indent=2) + "\n")
    full_workflow = args.generation_run / "full/results/variants/raw/cohort.raw.vcf.gz"
    baseline_audit = compare_callsets(
        full_workflow, root / "baseline/gather/cohort.raw.vcf.gz", contigs
    )
    if any(
        baseline_audit[key]
        for key in ("different_records", "left_only_variants", "right_only_variants")
    ):
        raise ValueError(
            "Targeted baseline differs from the complete synthetic production workflow"
        )
    print(
        "PASS: actual pinned tools exercise variant arguments, sample map, splitting, grouping, Reblock and consolidate"
    )
    print("Synthetic execution only; this does not open any real cohort/lineage gate")


if __name__ == "__main__":
    main()
