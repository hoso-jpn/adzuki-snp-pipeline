#!/usr/bin/env python3
"""Plan bounded research execution; no plan certifies accuracy or a service SLA."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

from hash_input_fastqs import INPUT_PROVENANCE_COLUMNS

SYNTHETIC_REFERENCE_SHA256 = "716fcd6b16040b4220d45ae93eeb45d3b9ffce5c8a413d869c921669aa53f889"
OBSERVED_REFERENCE_SHA256 = "e9838db1b048b54b21534285aaf95eae64cb3019d85af270b44e20b3d545f383"
# Historical observations, not guaranteed capacities for a different cohort.
OBSERVED_FASTQ_BYTES = 207_358_941_216
MAX_INTERVALS = 10_000
_CHECKSUM = re.compile(r"sha256:[0-9a-f]{64}\Z")


def bind_input_provenance(samples: list[dict], provenance_files: list[Path]) -> list[dict]:
    """Reuse measured FASTQ hashes, pairing by sample/read group rather than order."""
    provenance = {}
    for path in provenance_files:
        for line in path.read_text().splitlines():
            values = line.split("\t")
            if len(values) != len(INPUT_PROVENANCE_COLUMNS):
                raise ValueError("invalid input provenance row")
            row = dict(zip(INPUT_PROVENANCE_COLUMNS, values, strict=True))
            row.pop("rank")
            key = (row["sample_id"], row["read_group_id"])
            if not all(key) or key in provenance:
                raise ValueError("duplicate or invalid input provenance identity")
            provenance[key] = row
    expected = [(row["sample_id"], row["read_group_id"]) for row in samples]
    if len(set(expected)) != len(expected) or set(expected) != set(provenance):
        raise ValueError("input provenance must match each planned read group exactly")
    return [{**row, **provenance[(row["sample_id"], row["read_group_id"])]} for row in samples]


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def make_intervals(contigs: list[tuple[str, int]], window: int) -> list[dict]:
    if window < 1 or window > 20_000_000:
        raise ValueError(
            "window must be between 1 and 20000000 bp; 20 Mb is not a universal safe size"
        )
    names = set()
    for contig, length in contigs:
        if not contig or contig in names or length < 1 or any(c in contig for c in ":\t\r\n"):
            raise ValueError("invalid/duplicate contig or non-positive length")
        names.add(contig)
    if not names:
        raise ValueError("no reference intervals")
    count = sum((length + window - 1) // window for _, length in contigs)
    if count > MAX_INTERVALS:
        raise ValueError(f"interval count {count} exceeds the planning limit {MAX_INTERVALS}")
    result = []
    for contig, length in contigs:
        for start in range(1, length + 1, window):
            rank = len(result)
            result.append(
                {
                    "id": f"interval_{rank + 1:06d}",
                    "rank": rank,
                    "contig": contig,
                    "start": start,
                    "end": min(start + window - 1, length),
                }
            )
    return result


def build_plan(
    *,
    contigs: list[tuple[str, int]],
    samples: list[dict],
    reference: dict,
    pipeline_sha: str,
    ploidy: int,
    window: int,
    task_memory_gib: float,
    concurrency: int,
    memory_budget_gib: float,
    free_bytes: int,
    storage_multiplier: float = 6,
    review: dict | None = None,
) -> dict:
    if not samples or any(row.get("bytes", -1) < 0 or not row.get("sample_id") for row in samples):
        raise ValueError("invalid input-size summary")
    if any(
        not isinstance(row.get(field), str) or not _CHECKSUM.fullmatch(row[field])
        for row in samples
        for field in ("fastq_1_checksum", "fastq_2_checksum")
    ):
        raise ValueError("planned inputs require measured FASTQ SHA-256 checksums")
    identities = [(row["sample_id"], row["read_group_id"]) for row in samples]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate planned read group")
    if (
        task_memory_gib <= 0
        or concurrency < 1
        or memory_budget_gib <= 0
        or storage_multiplier < 1
        or not all(
            math.isfinite(v) for v in (task_memory_gib, memory_budget_gib, storage_multiplier)
        )
    ):
        raise ValueError("invalid resource budget")
    intervals = make_intervals(contigs, window)
    sample_ids = sorted({row["sample_id"] for row in samples})
    input_bytes = sum(row["bytes"] for row in samples)
    summary_sha = digest(sorted(samples, key=lambda row: (row["sample_id"], row["read_group_id"])))
    review_binding = {
        "purpose": "public_scale_validation",
        "pipeline_git_commit": pipeline_sha,
        "reference_fasta_sha256": reference["fasta"]["sha256"],
        "reference_bundle_fingerprint": reference["fingerprint"],
        "input_summary_sha256": summary_sha,
        "input_binding": "measured_FASTQ_SHA256_and_read_group_metadata",
        "sample_count": len(sample_ids),
        "ploidy": ploidy,
        "window_bp": window,
        "task_memory_gib": task_memory_gib,
        "concurrency": concurrency,
        "memory_budget_gib": memory_budget_gib,
    }
    if review is not None and not isinstance(review, dict):
        raise ValueError("review must be a JSON object")
    review_valid = (
        review is not None
        and all(review.get(k) == v for k, v in review_binding.items())
        and bool(review.get("reviewer"))
        and bool(review.get("rationale"))
    )
    reasons = []
    if len(sample_ids) > 51:
        reasons.append("more_than_51_samples_unvalidated; 327_sample_NO_GO_unchanged")
    elif len(sample_ids) > 20 and not review_valid:
        reasons.append("21_to_51_samples_require_input_bound_public_validation_review")
    if (
        reference["fasta"]["sha256"]
        not in {
            SYNTHETIC_REFERENCE_SHA256,
            OBSERVED_REFERENCE_SHA256,
        }
        and not review_valid
    ):
        reasons.append("unobserved_reference_requires_input_bound_validation_review")
    if input_bytes > OBSERVED_FASTQ_BYTES and not review_valid:
        reasons.append("input_larger_than_observed_requires_input_bound_validation_review")
    if sum(length for _, length in contigs) > 500_000_000:
        reasons.append("reference_larger_than_research_envelope")
    if ploidy != 2 and sum(length for _, length in contigs) > 1_000_000:
        reasons.append("non_diploid_real_reference_not_validated")
    if task_memory_gib * concurrency > memory_budget_gib:
        reasons.append("genotype_concurrency_exceeds_launch_memory_budget")
    # The multiplier is an explicit planning assumption, not a measured upper bound.
    required_bytes = math.ceil(input_bytes * storage_multiplier)
    if required_bytes > free_bytes:
        reasons.append("insufficient_snapshot_free_space")
    return {
        "schema_version": 1,
        "contract": "joint_execution_plan_v1",
        "status": "blocked" if reasons else "eligible_for_research_execution",
        "commercial_validation": "not_established",
        "blocking_reasons": reasons,
        "pipeline_git_commit": pipeline_sha,
        "reference_bundle_fingerprint": reference["fingerprint"],
        "sample_count": len(sample_ids),
        "canonical_sample_order": sample_ids,
        "input_summary_sha256": summary_sha,
        "intervals": intervals,
        "resources": {
            "genotype_task_memory_gib": task_memory_gib,
            "max_concurrency": concurrency,
            "launch_memory_budget_gib": memory_budget_gib,
            "heap_fraction": 0.8,
            "interval_count": len(intervals),
            "interval_count_limit": MAX_INTERVALS,
            "input_bytes_observed": input_bytes,
            "free_bytes_snapshot": free_bytes,
            "retained_bytes_estimated": required_bytes,
            "storage_multiplier_assumption": storage_multiplier,
        },
        "evidence": {
            "full_pipeline_samples_observed": 20,
            "targeted_joint_samples_observed": 51,
            "new_production_profile_real_replay": "not_completed",
            "sample327": "NO_GO",
            "automatic_cleanup": False,
            "historical_reference_fasta_sha256": OBSERVED_REFERENCE_SHA256,
            "historical_51_sample_fastq_bytes": OBSERVED_FASTQ_BYTES,
        },
        "review_template": {**review_binding, "reviewer": "", "rationale": ""},
        "review_sha256": digest(review) if review else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fai", "samples", "reference-manifest", "output", "intervals"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--pipeline-sha", required=True)
    parser.add_argument("--ploidy", type=int, required=True)
    parser.add_argument("--window", type=int, default=20_000_000)
    parser.add_argument("--task-memory-gib", type=float, default=16)
    parser.add_argument("--concurrency", type=int, default=3)
    parser.add_argument("--memory-budget-gib", type=float, default=110)
    parser.add_argument("--output-free-bytes", type=int, required=True)
    parser.add_argument("--review", type=Path)
    parser.add_argument("--input-provenance", type=Path, nargs="+", required=True)
    parser.add_argument("--require-eligible", action="store_true")
    args = parser.parse_args()
    try:
        contigs = [
            (row.split("\t")[0], int(row.split("\t")[1]))
            for row in args.fai.read_text().splitlines()
        ]
        samples = bind_input_provenance(
            [json.loads(row) for row in args.samples.read_text().splitlines()],
            args.input_provenance,
        )
        plan = build_plan(
            contigs=contigs,
            samples=samples,
            reference=json.loads(args.reference_manifest.read_text()),
            pipeline_sha=args.pipeline_sha,
            ploidy=args.ploidy,
            window=args.window,
            task_memory_gib=args.task_memory_gib,
            concurrency=args.concurrency,
            memory_budget_gib=args.memory_budget_gib,
            free_bytes=min(shutil.disk_usage(".").free, args.output_free_bytes),
            review=json.loads(args.review.read_text()) if args.review else None,
        )
        args.output.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
        args.intervals.write_text(
            "".join(
                f"{r['id']}\t{r['rank']}\t{r['contig']}\t{r['start']}\t{r['end']}\n"
                for r in plan["intervals"]
            )
        )
        if args.require_eligible and plan["blocking_reasons"]:
            parser.exit(2, "run plan blocked: " + "; ".join(plan["blocking_reasons"]) + "\n")
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"run plan rejected: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
