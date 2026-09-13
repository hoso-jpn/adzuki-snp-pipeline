#!/usr/bin/env python3
"""Classify every differing record between two Issue45 callsets against a fixed gate.

GATK 4.6.2.0 QualByDepth replaces a raw QD at or above MAX_QD_BEFORE_FIXING
(35.0) with IDEAL_HIGH_QD (30.0) + nextGaussian() * JITTER_SIGMA (3.0).
Re-partitioning the genome changes that random stream, so QD alone may differ
between two runs whose every other field agrees. The inclusive comparison is
taken from the pinned jar's own bytecode, not from documentation.

Why reconstructing raw QD from the final VCF is imperfect: GenotypeGVCFsEngine
passes likelihoods=null when it annotates QD, so the depth QualByDepth.getDepth()
used is the genotype AD, or DP when AD is absent -- exactly what `raw_qd` recomputes.
The discrepancy is not a read-likelihood depth this reconstruction cannot see.
It is that GenotypeGVCFs performs genotype cleanup *after* INFO annotation,
migrating MIN_DP into DP and synthesising AD from DP where AD is absent, so the
AD/DP finally serialized need not match the genotype state at the moment QD was
annotated. Serialized QUAL precision adds further reconstruction error.

The comparison helper classifies only raw QD >= 35.01, holding a 0.01 margin for
that reconstruction error. Records inside the margin are unclassified by
construction rather than scientifically unexplained, and this script tests that
claim per record instead of asserting it.

A record whose QD differs while QUAL, AD, GT and every other field agree can only
have been jittered on at least one side, since QD is otherwise deterministic.
When neither output QD equals the deterministic reconstruction, jitter on both
sides is the natural reading -- but the deterministic value is itself back-computed
from the final VCF, so this is strong circumstantial evidence, not proof. Such
records are therefore reported as consistent with the reconstruction limit and
never as directly explained.

This script is reporting-only. It never runs a benchmark tool and is excluded
from the execution lineage a resumed suite verifies.
"""

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter

# Verified against the pinned image's own bytecode rather than assumed. In
# gatk-package-4.6.2.0-local.jar, QualByDepth.fixTooHighQD compiles to
# `dcmpg; ifge` against 35.0d, so the jitter branch is taken when raw QD is
# greater than OR EQUAL TO 35.0, and QD is emitted through String.format("%.2f").
# The jitter is IDEAL_HIGH_QD + Random.nextGaussian() * JITTER_SIGMA.
GATK_MAX_QD_BEFORE_FIXING = 35.0
GATK_IDEAL_HIGH_QD = 30.0
GATK_JITTER_SIGMA = 3.0
HELPER_MARGIN = 35.01
JITTER_PLAUSIBLE_SIGMA = 5.0
DETERMINISTIC_TOLERANCE = 0.011  # QD is printed to 2 decimals


def positions(handle):
    current, group = None, {}
    for line in handle:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos = (f[0], int(f[1]))
        if current is not None and pos != current:
            yield current, group
            group = {}
        current = pos
        group[(f[3], f[4])] = f
    if current is not None:
        yield current, group


def info_of(fields):
    return dict(
        (item.split("=", 1) if "=" in item else (item, None)) for item in fields[7].split(";")
    )


def plausible_jitter(value):
    return abs(value - GATK_IDEAL_HIGH_QD) <= JITTER_PLAUSIBLE_SIGMA * GATK_JITTER_SIGMA


def classify(left, right, raw):
    """Return (category, rollup) for one QD-only differing record."""
    try:
        lq, rq = float(info_of(left)["QD"]), float(info_of(right)["QD"])
    except (KeyError, ValueError):
        return "non_finite_or_missing_annotation", "unresolved"
    if not (math.isfinite(lq) and math.isfinite(rq)):
        return "non_finite_or_missing_annotation", "unresolved"
    if raw is None:
        return "cannot_reconstruct_raw_qd", "unresolved"
    if not math.isfinite(raw):
        return "non_finite_or_missing_annotation", "unresolved"
    if raw >= HELPER_MARGIN:
        return "directly_explained_high_qd_jitter", "directly_explained"
    # Inside the helper's conservative margin, or below GATK's own threshold.
    jittered_both = (
        abs(lq - raw) > DETERMINISTIC_TOLERANCE and abs(rq - raw) > DETERMINISTIC_TOLERANCE
    )
    if (
        raw >= GATK_MAX_QD_BEFORE_FIXING
        and jittered_both
        and plausible_jitter(lq)
        and plausible_jitter(rq)
    ):
        # At or above GATK's own threshold but inside the helper's reconstruction
        # margin, with neither output equal to the deterministic value. The
        # comparison mirrors the pinned bytecode's `ifge`: a reconstruction
        # landing exactly on 35.0 is a jitter case for GATK, and treating it
        # otherwise would report a difference GATK's own predicate explains.
        return "reconstruction_margin_at_or_above_gatk_threshold", "reconstruction_limit_consistent"
    if jittered_both and plausible_jitter(lq) and plausible_jitter(rq):
        return "near_threshold_below_gatk_threshold", "unresolved"
    return "other_unresolved", "unresolved"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--left-label", default="E1")
    parser.add_argument("--right-label", default="E2a")
    parser.add_argument("--helpers", required=True)
    parser.add_argument("--per-record-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()

    sys.path.insert(0, args.helpers)
    from benchmark_tools import raw_qd

    categories, rollups = Counter(), Counter()
    non_qd_fields = Counter()
    filter_flips = []
    qd_min_seen = math.inf
    raw_band = []
    total_differing = 0

    with (
        gzip.open(args.left, "rt") as a,
        gzip.open(args.right, "rt") as b,
        open(args.per_record_tsv, "w", newline="") as out,
    ):
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                f"{args.left_label}_qd",
                f"{args.right_label}_qd",
                f"{args.left_label}_filter",
                f"{args.right_label}_filter",
                "reconstructed_raw_qd",
                "classification",
                "rollup",
                "qd_lt2_filter_flip",
            ]
        )
        left, right = positions(a), positions(b)
        x, y = next(left, None), next(right, None)
        while x is not None and y is not None:
            if x[0] != y[0]:
                if x[0] < y[0]:
                    x = next(left, None)
                else:
                    y = next(right, None)
                continue
            for key in sorted(x[1].keys() & y[1].keys()):
                lf, rf = x[1][key], y[1][key]
                if lf == rf:
                    continue
                total_differing += 1
                li, ri = info_of(lf), info_of(rf)
                differing = {
                    k for k in li.keys() | ri.keys() if k not in li or k not in ri or li[k] != ri[k]
                }
                qd_only = differing == {"QD"} and lf[:7] == rf[:7] and lf[8:] == rf[8:]
                if not qd_only:
                    for field in differing - {"QD"}:
                        non_qd_fields[field] += 1
                    for index, name in ((5, "QUAL"), (6, "FILTER"), (8, "FORMAT")):
                        if lf[index] != rf[index]:
                            non_qd_fields[name] += 1
                    if lf[9:] != rf[9:]:
                        non_qd_fields["SAMPLE_COLUMNS"] += 1
                    writer.writerow(
                        [
                            lf[0],
                            lf[1],
                            lf[3],
                            lf[4],
                            li.get("QD"),
                            ri.get("QD"),
                            lf[6],
                            rf[6],
                            "",
                            "not_qd_only",
                            "unresolved",
                            "",
                        ]
                    )
                    categories["not_qd_only"] += 1
                    rollups["unresolved"] += 1
                    continue
                try:
                    raw = raw_qd(lf)
                except (KeyError, ValueError, TypeError, IndexError):
                    raw = None
                category, rollup = classify(lf, rf, raw)
                categories[category] += 1
                rollups[rollup] += 1
                try:
                    lq, rq = float(li["QD"]), float(ri["QD"])
                    qd_min_seen = min(qd_min_seen, lq, rq)
                    flip = (lq < 2.0) != (rq < 2.0)
                except (KeyError, ValueError):
                    flip = False
                if flip:
                    filter_flips.append(
                        {
                            "chrom": lf[0],
                            "pos": lf[1],
                            "ref": lf[3],
                            "alt": lf[4],
                            "left_qd": li.get("QD"),
                            "right_qd": ri.get("QD"),
                            "left_filter": lf[6],
                            "right_filter": rf[6],
                        }
                    )
                if raw is not None and math.isfinite(raw) and raw < HELPER_MARGIN:
                    raw_band.append(raw)
                    writer.writerow(
                        [
                            lf[0],
                            lf[1],
                            lf[3],
                            lf[4],
                            li.get("QD"),
                            ri.get("QD"),
                            lf[6],
                            rf[6],
                            f"{raw:.6f}",
                            category,
                            rollup,
                            str(flip).lower(),
                        ]
                    )
            x, y = next(left, None), next(right, None)

    unresolved = rollups["unresolved"]
    if filter_flips or non_qd_fields:
        gate = "REJECT"
    elif unresolved == 0:
        gate = "ADOPT_CANDIDATE"
    else:
        gate = "DEFER"

    summary = {
        "comparison": f"{args.right_label} against {args.left_label}",
        "total_differing_records": total_differing,
        "categories": dict(categories),
        "rollup": {
            "directly_explained": rollups["directly_explained"],
            "reconstruction_limit_consistent": rollups["reconstruction_limit_consistent"],
            "unresolved": rollups["unresolved"],
        },
        "rollup_percent_of_differing": {
            key: (round(rollups[key] / total_differing * 100, 6) if total_differing else 0.0)
            for key in (
                "directly_explained",
                "reconstruction_limit_consistent",
                "unresolved",
            )
        },
        "scientific_semantics": {
            "qd_lt2_filter_membership_flips": len(filter_flips),
            "filter_flip_examples": filter_flips[:10],
            "non_qd_field_differences": dict(non_qd_fields),
            "minimum_qd_observed_among_differing_records": (
                None if qd_min_seen is math.inf else round(qd_min_seen, 2)
            ),
            "qd_lt2_filter_threshold": 2.0,
        },
        "reconstruction_margin_band": {
            "count": len(raw_band),
            "min_raw_qd": round(min(raw_band), 6) if raw_band else None,
            "max_raw_qd": round(max(raw_band), 6) if raw_band else None,
            "gatk_threshold": GATK_MAX_QD_BEFORE_FIXING,
            "helper_margin": HELPER_MARGIN,
            "all_at_or_above_gatk_threshold": (
                bool(raw_band) and min(raw_band) >= GATK_MAX_QD_BEFORE_FIXING
            ),
        },
        "gate": gate,
    }
    with open(args.summary_json, "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
