#!/usr/bin/env python3
"""Issue #64: how a DP/GQ genotype mask would change a real GS panel, over a threshold grid.

This is a **sensitivity analysis**, not an accuracy measurement: there is no
truth set, so nothing here says which threshold gives better genotypes. It
records what each threshold pair would remove and what that does to
missingness and allele counts, so a later adoption decision starts from
measured consequences instead of a number chosen in advance.

One streaming pass over the GS-eligible PASS VCF covers the whole grid. For
every call the build would encode, each configured field's status comes
from the shared evaluator in bin/gs_genotype_quality.py (evaluated with a
zero threshold, so a readable value is `pass` and everything else keeps its
own status). The value's rank against the sorted thresholds then says, for
every grid point at once, whether the call survives: it is kept at
(DP threshold i, GQ threshold j) exactly when both ranks reach i and j. Per
row the calls are folded into a small rank-by-rank table and suffix-summed,
so the work per row does not multiply by the number of grid points.

Calls whose DP or GQ cannot be read are kept here (the `unevaluated`
policy), and counted separately per field, so each grid point's masking is
driven by the thresholds alone. `--cross-check` compares one grid point with
the accounting of a real `build_gs_panel.py --genotype-quality-mask` run
under the equivalent policy, which ties this analysis to the production
evaluator rather than to a second copy of its rule.

Memory is one row plus per-sample and per-grid-point counters and
fixed-size histograms.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "bin"))

import gs_genotype_quality as quality  # noqa: E402
from build_gs_panel import FIXED_COLUMN_COUNT, classify_genotype  # noqa: E402

MISSINGNESS_BIN_EDGES = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
DELTA_AF_BIN_EDGES = (0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)
MAF_CUTOFFS = (0.01, 0.05)


def _bin(value: float, edges: tuple[float, ...]) -> int:
    """Index of the half-open-left bin: [0], (e0, e1], ..., with exact 0 in bin 0."""
    if value <= edges[0]:
        return 0
    return min(bisect.bisect_left(edges, value), len(edges) - 1)


def _bin_labels(edges: tuple[float, ...]) -> list[str]:
    return [f"=={edges[0]}"] + [f"({edges[k - 1]},{edges[k]}]" for k in range(1, len(edges))]


class Grid:
    def __init__(self, dp_thresholds: list[int], gq_thresholds: list[int]) -> None:
        self.dp = sorted(set(dp_thresholds))
        self.gq = sorted(set(gq_thresholds))
        # Index 0 is "no threshold on this field".
        self.ni = len(self.dp) + 1
        self.nj = len(self.gq) + 1

    def rank(self, thresholds: list[int], status: str, value: str | None) -> int:
        """How many thresholds (from the lowest) this call passes; unreadable passes all."""
        if status != quality.PASS:
            return len(thresholds)
        return bisect.bisect_right(thresholds, int(value))

    def label(self, i: int, j: int) -> str:
        dp = "none" if i == 0 else str(self.dp[i - 1])
        gq = "none" if j == 0 else str(self.gq[j - 1])
        return f"dp>={dp},gq>={gq}"


def _suffix(table: list[list[int]], ni: int, nj: int) -> list[list[int]]:
    """out[i][j] = sum of table[a][b] for a >= i and b >= j."""
    out = [[0] * (nj + 1) for _ in range(ni + 1)]
    for a in range(ni - 1, -1, -1):
        row, below, current = table[a], out[a + 1], out[a]
        for b in range(nj - 1, -1, -1):
            current[b] = row[b] + current[b + 1] + below[b] - below[b + 1]
    return out


def analyse(vcf: Path, grid: Grid) -> dict[str, object]:
    ni, nj = grid.ni, grid.nj
    status_policy = quality.GenotypeQualityPolicy(
        min_dp=0,
        min_gq=0,
        missing_format_field=quality.UNEVALUATED,
        missing_value=quality.UNEVALUATED,
        malformed_value=quality.UNEVALUATED,
    )
    points = [(i, j) for i in range(ni) for j in range(nj)]

    samples: tuple[str, ...] = ()
    variants = 0
    cells = 0
    shape_counts = {"missing": 0, "non_diploid": 0, "non_biallelic_index": 0}
    evaluated = 0
    field_status = {key: dict.fromkeys(quality.FIELD_STATUSES, 0) for key in ("DP", "GQ")}
    dp_values_hist: dict[int, int] = {}
    gq_values_hist: dict[int, int] = {}
    # Cohort corner table by original dosage.
    corner = {dosage: [[0] * nj for _ in range(ni)] for dosage in ("-1", "0", "1")}
    sample_corner: list[list[list[int]]] = []
    sample_original_missing: list[int] = []
    per_point = {
        point: {
            "variant_missingness_hist": [0] * len(MISSINGNESS_BIN_EDGES),
            "variants_fully_missing": 0,
            "variants_ac_changed": 0,
            "variants_an_changed": 0,
            "variants_polymorphic_to_monomorphic": 0,
            "variants_an_zero": 0,
            "an_sum": 0,
            "delta_af_hist": [0] * len(DELTA_AF_BIN_EDGES),
            "maf_crossings": {
                str(c): {"above_to_below": 0, "below_to_above": 0} for c in MAF_CUTOFFS
            },
        }
        for point in points
    }

    with gzip.open(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = tuple(fields[FIXED_COLUMN_COUNT:])
                sample_corner = [[[0] * nj for _ in range(ni)] for _ in samples]
                sample_original_missing = [0] * len(samples)
                continue
            variants += 1
            keys = fields[8].split(":")
            gt_index = keys.index("GT")
            evaluator = quality.RowQualityEvaluator(status_policy, fields[8])
            dp_index = keys.index("DP") if "DP" in keys else None
            gq_index = keys.index("GQ") if "GQ" in keys else None

            row_corner_n = [[0] * nj for _ in range(ni)]
            row_corner_ac = [[0] * nj for _ in range(ni)]
            fixed_an = fixed_ac = original_missing = 0
            for position, sample_field in enumerate(fields[FIXED_COLUMN_COUNT:]):
                cells += 1
                subfields = sample_field.split(":")
                gt = subfields[gt_index]
                cell = classify_genotype(gt)
                if cell.category != "standard":
                    shape_counts[cell.category] += 1
                    original_missing += 1
                    sample_original_missing[position] += 1
                    ac, an = quality.allele_counts([gt])
                    fixed_ac += ac
                    fixed_an += an
                    continue
                evaluated += 1
                verdict = evaluator.evaluate(subfields)
                dp_status, gq_status = verdict.statuses
                field_status["DP"][dp_status] += 1
                field_status["GQ"][gq_status] += 1
                dp_value = subfields[dp_index] if dp_status == quality.PASS else None
                gq_value = subfields[gq_index] if gq_status == quality.PASS else None
                if dp_value is not None:
                    v = int(dp_value)
                    dp_values_hist[min(v, 100)] = dp_values_hist.get(min(v, 100), 0) + 1
                if gq_value is not None:
                    v = int(gq_value)
                    gq_values_hist[min(v, 99)] = gq_values_hist.get(min(v, 99), 0) + 1
                a = grid.rank(grid.dp, dp_status, dp_value)
                b = grid.rank(grid.gq, gq_status, gq_value)
                alt = 1 if cell.dosage == "0" else (2 if cell.dosage == "1" else 0)
                row_corner_n[a][b] += 1
                row_corner_ac[a][b] += alt
                corner[cell.dosage][a][b] += 1
                sample_corner[position][a][b] += 1

            kept_n = _suffix(row_corner_n, ni, nj)
            kept_ac = _suffix(row_corner_ac, ni, nj)
            sample_count = len(samples)
            ac_pre = fixed_ac + kept_ac[0][0]
            an_pre = fixed_an + 2 * kept_n[0][0]
            af_pre = ac_pre / an_pre if an_pre else None
            for i, j in points:
                stats = per_point[(i, j)]
                kept = kept_n[i][j]
                missing = original_missing + (kept_n[0][0] - kept)
                stats["variant_missingness_hist"][
                    _bin(missing / sample_count, MISSINGNESS_BIN_EDGES)
                ] += 1
                if missing == sample_count:
                    stats["variants_fully_missing"] += 1
                ac_post = fixed_ac + kept_ac[i][j]
                an_post = fixed_an + 2 * kept
                stats["an_sum"] += an_post
                if ac_post != ac_pre:
                    stats["variants_ac_changed"] += 1
                if an_post != an_pre:
                    stats["variants_an_changed"] += 1
                if an_post == 0:
                    stats["variants_an_zero"] += 1
                if 0 < ac_pre < an_pre and (ac_post == 0 or ac_post == an_post):
                    stats["variants_polymorphic_to_monomorphic"] += 1
                if af_pre is not None and an_post:
                    af_post = ac_post / an_post
                    stats["delta_af_hist"][_bin(abs(af_post - af_pre), DELTA_AF_BIN_EDGES)] += 1
                    maf_pre, maf_post = min(af_pre, 1 - af_pre), min(af_post, 1 - af_post)
                    for cutoff in MAF_CUTOFFS:
                        crossing = stats["maf_crossings"][str(cutoff)]
                        if maf_pre >= cutoff > maf_post:
                            crossing["above_to_below"] += 1
                        elif maf_pre < cutoff <= maf_post:
                            crossing["below_to_above"] += 1

    cohort_kept = {dosage: _suffix(table, ni, nj) for dosage, table in corner.items()}
    sample_kept = [_suffix(table, ni, nj) for table in sample_corner]
    original_missing_total = sum(shape_counts.values())
    grid_results = []
    for i, j in points:
        kept = {dosage: cohort_kept[dosage][i][j] for dosage in corner}
        masked_total = evaluated - sum(kept.values())
        # Reason split: masked by DP alone / GQ alone / both, from the corner tables.
        dp_fail = sum(
            corner[d][a][b]
            for d in corner
            for a in range(ni)
            for b in range(nj)
            if a < i and b >= j
        )
        gq_fail = sum(
            corner[d][a][b]
            for d in corner
            for a in range(ni)
            for b in range(nj)
            if a >= i and b < j
        )
        both_fail = sum(
            corner[d][a][b] for d in corner for a in range(ni) for b in range(nj) if a < i and b < j
        )
        stats = per_point[(i, j)]
        sample_missing = [
            (sample_original_missing[s] + (sample_kept[s][0][0] - sample_kept[s][i][j])) / variants
            for s in range(len(samples))
        ]
        ordered = sorted(sample_missing)
        grid_results.append(
            {
                "policy": grid.label(i, j),
                "min_dp": None if i == 0 else grid.dp[i - 1],
                "min_gq": None if j == 0 else grid.gq[j - 1],
                "masked_calls": masked_total,
                "masked_low_dp_only": dp_fail,
                "masked_low_gq_only": gq_fail,
                "masked_low_dp_and_low_gq": both_fail,
                "masked_by_original_dosage": {
                    "hom_ref": evaluated_by(corner, "-1") - kept["-1"],
                    "het": evaluated_by(corner, "0") - kept["0"],
                    "hom_alt": evaluated_by(corner, "1") - kept["1"],
                },
                "post_mask_missing_cells": original_missing_total + masked_total,
                "post_mask_missingness": (original_missing_total + masked_total) / cells,
                "sample_missingness": {
                    "min": ordered[0],
                    "median": ordered[len(ordered) // 2],
                    "max": ordered[-1],
                    "per_sample": dict(zip(samples, (round(value, 6) for value in sample_missing))),
                },
                "variant_missingness_histogram": dict(
                    zip(_bin_labels(MISSINGNESS_BIN_EDGES), stats["variant_missingness_hist"])
                ),
                "variants_fully_missing": stats["variants_fully_missing"],
                "variants_ac_changed": stats["variants_ac_changed"],
                "variants_an_changed": stats["variants_an_changed"],
                "variants_an_zero": stats["variants_an_zero"],
                "variants_polymorphic_to_monomorphic": stats["variants_polymorphic_to_monomorphic"],
                "mean_an": stats["an_sum"] / variants if variants else None,
                "abs_delta_af_histogram": dict(
                    zip(_bin_labels(DELTA_AF_BIN_EDGES), stats["delta_af_hist"])
                ),
                "maf_threshold_crossings": stats["maf_crossings"],
            }
        )

    return {
        "variants": variants,
        "samples": len(samples),
        "total_genotype_cells": cells,
        "originally_missing_or_non_standard": {
            "missing_calls": shape_counts["missing"],
            "non_diploid_calls": shape_counts["non_diploid"],
            "non_biallelic_index_calls": shape_counts["non_biallelic_index"],
        },
        "evaluated_calls": evaluated,
        "field_status_threshold_independent": field_status,
        "field_status_note": (
            "statuses come from the shared evaluator at threshold 0, so 'pass' means a readable "
            "value and below_threshold is always 0 here; per-threshold masking is in 'grid'"
        ),
        "dp_value_histogram_capped_100": dict(sorted(dp_values_hist.items())),
        "gq_value_histogram_capped_99": dict(sorted(gq_values_hist.items())),
        "unreadable_policy_in_this_analysis": "unevaluated (kept and counted per field)",
        "grid": grid_results,
    }


def evaluated_by(corner: dict[str, list[list[int]]], dosage: str) -> int:
    return sum(sum(row) for row in corner[dosage])


def cross_check(
    result: dict[str, object], accounting: Path, min_dp: int | None, min_gq: int | None
) -> dict[str, object]:
    """Compare one grid point with a real build's accounting under the same policy."""
    reported = {}
    for line in accounting.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        reported[parts[1]] = parts[2]
    point = next(p for p in result["grid"] if p["min_dp"] == min_dp and p["min_gq"] == min_gq)
    checks = {
        "total_genotype_cells": (
            result["total_genotype_cells"],
            int(reported["total_genotype_cells"]),
        ),
        "quality_evaluated_calls": (
            result["evaluated_calls"],
            int(reported["quality_evaluated_calls"]),
        ),
        "quality_masked_calls": (point["masked_calls"], int(reported["quality_masked_calls"])),
        "quality_masked_hom_ref_calls": (
            point["masked_by_original_dosage"]["hom_ref"],
            int(reported["quality_masked_hom_ref_calls"]),
        ),
        "quality_masked_het_calls": (
            point["masked_by_original_dosage"]["het"],
            int(reported["quality_masked_het_calls"]),
        ),
        "quality_masked_hom_alt_calls": (
            point["masked_by_original_dosage"]["hom_alt"],
            int(reported["quality_masked_hom_alt_calls"]),
        ),
        "masked_low_dp_only": (
            point["masked_low_dp_only"],
            int(reported.get("quality_masked_reason.dp_low.gq_none", 0)),
        ),
        "masked_low_gq_only": (
            point["masked_low_gq_only"],
            int(reported.get("quality_masked_reason.dp_none.gq_low", 0)),
        ),
        "masked_low_dp_and_low_gq": (
            point["masked_low_dp_and_low_gq"],
            int(reported.get("quality_masked_reason.dp_low.gq_low", 0)),
        ),
    }
    return {
        "policy": point["policy"],
        "checks": {
            name: {"sensitivity": a, "build": b, "equal": a == b} for name, (a, b) in checks.items()
        },
        "all_equal": all(a == b for a, b in checks.values()),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gs-pass-vcf", required=True, type=Path)
    parser.add_argument(
        "--dp-thresholds", required=True, help="comma-separated, e.g. 3,5,8,10,15,20"
    )
    parser.add_argument("--gq-thresholds", required=True, help="comma-separated, e.g. 10,20,30")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cross-check-accounting", type=Path)
    parser.add_argument("--cross-check-min-dp", type=int)
    parser.add_argument("--cross-check-min-gq", type=int)
    args = parser.parse_args(argv)

    grid = Grid(
        [int(x) for x in args.dp_thresholds.split(",")],
        [int(x) for x in args.gq_thresholds.split(",")],
    )
    result = analyse(args.gs_pass_vcf, grid)
    if args.cross_check_accounting:
        result["cross_check_with_build_gs_panel"] = cross_check(
            result, args.cross_check_accounting, args.cross_check_min_dp, args.cross_check_min_gq
        )
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    check = result.get("cross_check_with_build_gs_panel")
    return 0 if check is None or check["all_equal"] else 2


if __name__ == "__main__":
    sys.exit(main())
