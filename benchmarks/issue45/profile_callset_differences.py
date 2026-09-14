#!/usr/bin/env python3
"""Profile every difference between two joint callsets, field by field.

Written for E3, where ReblockGVCFs deliberately changes the reference-confidence
representation of its inputs, so a byte-identical callset is not the expectation
and its absence is not a failure. The opposite error matters just as much:
waving every difference through as "expected from reblocking" would hide a real
semantic change inside a representational one.

This script therefore decides nothing. It measures the difference profile at
field granularity and separates the categories that carry scientific meaning --
the variant set, per-sample genotypes, FILTER membership -- from annotation-only
changes, then reports which annotations actually moved and by how much, so each
can be judged against what the tool is documented to do rather than by assertion.
"""

import argparse
import gzip
import json
import math
from collections import Counter, defaultdict


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


def genotypes(fields):
    keys = fields[8].split(":")
    index = keys.index("GT") if "GT" in keys else None
    if index is None:
        return None
    return [sample.split(":")[index] for sample in fields[9:]]


def called(gt):
    return gt not in (".", "./.", ".|.")


def numeric_delta(left, right, store, field):
    try:
        a, b = float(left), float(right)
    except (TypeError, ValueError):
        store[field].append(None)
        return
    if math.isfinite(a) and math.isfinite(b):
        store[field].append(b - a)
    else:
        store[field].append(None)


def summarise(values):
    numbers = [v for v in values if v is not None]
    out = {"count": len(values), "unevaluable": len(values) - len(numbers)}
    if numbers:
        numbers.sort()
        out.update(
            min=round(numbers[0], 6),
            max=round(numbers[-1], 6),
            median=round(numbers[len(numbers) // 2], 6),
            mean=round(sum(numbers) / len(numbers), 6),
            increased=sum(1 for n in numbers if n > 0),
            decreased=sum(1 for n in numbers if n < 0),
            unchanged_value=sum(1 for n in numbers if n == 0),
        )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--left-label", default="left")
    parser.add_argument("--right-label", default="right")
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--examples", type=int, default=8)
    args = parser.parse_args()

    counts = Counter()
    info_fields, format_fields, column_fields = Counter(), Counter(), Counter()
    info_deltas = defaultdict(list)
    gt_transitions = Counter()
    genotype_change_records = 0
    filter_change_records = 0
    variant_set_examples, genotype_examples, filter_examples = [], [], []

    with gzip.open(args.left, "rt") as a, gzip.open(args.right, "rt") as b:
        left, right = positions(a), positions(b)
        x, y = next(left, None), next(right, None)
        while x is not None or y is not None:
            if y is None or (x is not None and x[0] < y[0]):
                for key, fields in x[1].items():
                    counts["left_only_variants"] += 1
                    if len(variant_set_examples) < args.examples:
                        variant_set_examples.append(
                            {
                                "side": args.left_label,
                                "chrom": fields[0],
                                "pos": fields[1],
                                "ref": fields[3],
                                "alt": fields[4],
                                "qual": fields[5],
                                "filter": fields[6],
                            }
                        )
                x = next(left, None)
                continue
            if x is None or y[0] < x[0]:
                for key, fields in y[1].items():
                    counts["right_only_variants"] += 1
                    if len(variant_set_examples) < args.examples:
                        variant_set_examples.append(
                            {
                                "side": args.right_label,
                                "chrom": fields[0],
                                "pos": fields[1],
                                "ref": fields[3],
                                "alt": fields[4],
                                "qual": fields[5],
                                "filter": fields[6],
                            }
                        )
                y = next(right, None)
                continue
            for key in x[1].keys() - y[1].keys():
                counts["left_only_variants"] += 1
            for key in y[1].keys() - x[1].keys():
                counts["right_only_variants"] += 1
            for key in x[1].keys() & y[1].keys():
                lf, rf = x[1][key], y[1][key]
                counts["shared_variants"] += 1
                if lf == rf:
                    counts["identical_records"] += 1
                    continue
                counts["different_records"] += 1

                if lf[6] != rf[6]:
                    filter_change_records += 1
                    column_fields["FILTER"] += 1
                    if len(filter_examples) < args.examples:
                        filter_examples.append(
                            {
                                "chrom": lf[0],
                                "pos": lf[1],
                                "ref": lf[3],
                                "alt": lf[4],
                                args.left_label: lf[6],
                                args.right_label: rf[6],
                            }
                        )
                for index, name in ((2, "ID"), (5, "QUAL"), (8, "FORMAT_order")):
                    if lf[index] != rf[index]:
                        column_fields[name] += 1

                lg, rg = genotypes(lf), genotypes(rf)
                if lg is None or rg is None or len(lg) != len(rg):
                    counts["genotype_unevaluable"] += 1
                elif lg != rg:
                    genotype_change_records += 1
                    for a_gt, b_gt in zip(lg, rg):
                        if a_gt != b_gt:
                            gt_transitions[f"{a_gt}->{b_gt}"] += 1
                    if len(genotype_examples) < args.examples:
                        genotype_examples.append(
                            {
                                "chrom": lf[0],
                                "pos": lf[1],
                                "ref": lf[3],
                                "alt": lf[4],
                                "changed_samples": sum(1 for p, q in zip(lg, rg) if p != q),
                                "left_called": sum(1 for g in lg if called(g)),
                                "right_called": sum(1 for g in rg if called(g)),
                            }
                        )

                li, ri = info_of(lf), info_of(rf)
                for field in li.keys() | ri.keys():
                    lv, rv = li.get(field), ri.get(field)
                    if lv == rv and field in li and field in ri:
                        continue
                    if field not in li:
                        info_fields[field + " (added)"] += 1
                    elif field not in ri:
                        info_fields[field + " (removed)"] += 1
                    else:
                        info_fields[field] += 1
                        numeric_delta(lv, rv, info_deltas, field)

                lk, rk = lf[8].split(":"), rf[8].split(":")
                for ls, rs in zip(lf[9:], rf[9:]):
                    lm = dict(zip(lk, ls.split(":")))
                    rm = dict(zip(rk, rs.split(":")))
                    for field in lm.keys() | rm.keys():
                        if lm.get(field) != rm.get(field):
                            format_fields[field] += 1
            x, y = next(left, None), next(right, None)

    semantic = {
        "variant_set_changed": counts["left_only_variants"] + counts["right_only_variants"],
        "records_with_genotype_change": genotype_change_records,
        "records_with_filter_change": filter_change_records,
        "genotype_unevaluable": counts["genotype_unevaluable"],
    }
    summary = {
        "comparison": f"{args.right_label} against {args.left_label}",
        "totals": {
            "shared_variants": counts["shared_variants"],
            "identical_records": counts["identical_records"],
            "different_records": counts["different_records"],
            "left_only_variants": counts["left_only_variants"],
            "right_only_variants": counts["right_only_variants"],
        },
        "scientific_semantic_categories": semantic,
        "no_semantic_difference_observed": all(value == 0 for value in semantic.values()),
        "annotation_differences": {
            "info_fields": dict(info_fields.most_common()),
            "format_fields": dict(format_fields.most_common()),
            "column_fields": dict(column_fields.most_common()),
        },
        "info_numeric_deltas": {
            field: summarise(values) for field, values in sorted(info_deltas.items())
        },
        "genotype_transitions": dict(gt_transitions.most_common(20)),
        "examples": {
            "variant_set": variant_set_examples,
            "genotype_change": genotype_examples,
            "filter_change": filter_examples,
        },
        "interpretation_note": (
            "This profile decides nothing. A difference in the variant set, in any "
            "per-sample genotype, or in FILTER membership is a scientific-semantic "
            "difference regardless of its cause. Annotation-only differences still "
            "require a stated reason per field; none is excused merely for appearing "
            "in a reblocked comparison."
        ),
    }
    with open(args.summary_json, "w") as handle:
        json.dump(summary, handle, indent=2)
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "totals",
                    "scientific_semantic_categories",
                    "no_semantic_difference_observed",
                    "annotation_differences",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
