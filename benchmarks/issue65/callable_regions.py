#!/usr/bin/env python3
"""Per-sample and cohort callable regions from existing alignments (Issue #65).

Input is a gzip-compressed stream from

    samtools depth -a -Q <mapq> -q <baseq> -r <region> -f <bam list>

(one column per sample, every position of the region, samtools' default read
flag filter: unmapped, secondary, QC-fail and duplicate reads are not counted).
It is read twice: first for each sample's depth distribution over the region,
then to label positions. Outputs, for that region:

* `cohort_callable.bed` / `cohort_non_callable.bed` -- a partition of the region;
* `median_depth_*.bed` -- a partition of the region by the cohort median depth;
* `callable_summary.json` -- the rule, per-sample figures, bases and hashes.

Rule, recorded in the output:

* a sample is callable at a position when `min_depth <= depth <= max_depth`,
  with `max_depth = max_depth_factor x` the sample's median depth over the
  region (excess depth marks collapsed repeats and duplications);
* the cohort is callable at a position when at least `min_sample_fraction` of
  the samples are callable there. N reference bases have depth 0, so they are
  never callable when `min_depth >= 1`.

"Callable" states that the alignments meet these depth and quality conditions.
It is not a statement that genotypes there are accurate. These are benchmark
stratification parameters, not pipeline defaults.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

DEPTH_BIN_EDGES = (3, 6, 10, 15, 25)
HIST_CAP = 10000


def depth_label(value: float) -> str:
    lower = 0
    for edge in DEPTH_BIN_EDGES:
        if value < edge:
            return f"median_depth_{lower:02d}_{edge:02d}"
        lower = edge
    return f"median_depth_{lower:02d}_plus"


def all_depth_labels() -> list[str]:
    labels, lower = [], 0
    for edge in DEPTH_BIN_EDGES:
        labels.append(f"median_depth_{lower:02d}_{edge:02d}")
        lower = edge
    return [*labels, f"median_depth_{lower:02d}_plus"]


def _quantile(hist: list[int], q: float) -> int:
    target, running = q * sum(hist), 0
    for depth, count in enumerate(hist):
        running += count
        if count and running >= target:
            return depth
    return 0


def _rows(path: Path, contig: str, start: int, end: int, samples: int):
    expected = start
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != contig:
                raise SystemExit(f"unexpected contig {fields[0]!r}")
            pos = int(fields[1])
            if pos != expected:
                raise SystemExit(
                    f"depth stream is not contiguous at {pos} (expected {expected}); use samtools depth -a"
                )
            if len(fields) - 2 != samples:
                raise SystemExit(f"{len(fields) - 2} depth columns at {pos}, expected {samples}")
            expected += 1
            yield pos, [int(x) for x in fields[2:]]
    if expected != end + 1:
        raise SystemExit(f"depth stream ended at {expected - 1}, expected {end}")


def build(
    *,
    depth: Path,
    samples: list[str],
    region: str,
    min_depth: int,
    max_depth_factor: float,
    min_sample_fraction: float,
    mapq: int,
    baseq: int,
    output: Path,
) -> dict[str, object]:
    if min_depth < 1:
        raise SystemExit("--min-depth must be at least 1")
    contig, span = region.split(":")
    start, end = (int(x) for x in span.replace(",", "").split("-"))
    region_bases = end - start + 1

    hists = [[0] * (HIST_CAP + 1) for _ in samples]
    for _pos, depths in _rows(depth, contig, start, end, len(samples)):
        for i, value in enumerate(depths):
            hists[i][min(value, HIST_CAP)] += 1
    medians = [_quantile(h, 0.5) for h in hists]
    max_depths = [max_depth_factor * m for m in medians]
    need = min_sample_fraction * len(samples)

    output.mkdir(parents=True, exist_ok=False)
    labels = ["cohort_callable", "cohort_non_callable", *all_depth_labels()]
    handles = {label: (output / f"{label}.bed").open("w") for label in labels}
    bases = dict.fromkeys(labels, 0)
    sample_callable = [0] * len(samples)
    runs: dict[str, list] = {}

    def close_run(key: str) -> None:
        label, first, last = runs[key]
        handles[label].write(f"{contig}\t{first - 1}\t{last}\n")
        bases[label] += last - first + 1

    try:
        for pos, depths in _rows(depth, contig, start, end, len(samples)):
            callable_count = 0
            for i, value in enumerate(depths):
                if min_depth <= value <= max_depths[i]:
                    callable_count += 1
                    sample_callable[i] += 1
            ordered = sorted(depths)
            mid = len(ordered) // 2
            median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2
            for key, label in (
                (
                    "callable",
                    "cohort_callable" if callable_count >= need else "cohort_non_callable",
                ),
                ("depth", depth_label(median)),
            ):
                run = runs.get(key)
                if run is not None and run[0] == label and run[2] == pos - 1:
                    run[2] = pos
                else:
                    if run is not None:
                        close_run(key)
                    runs[key] = [label, pos, pos]
        for key in runs:
            close_run(key)
    finally:
        for handle in handles.values():
            handle.close()
    if bases["cohort_callable"] + bases["cohort_non_callable"] != region_bases:
        raise SystemExit("callable labels do not partition the region")
    if sum(bases[label] for label in all_depth_labels()) != region_bases:
        raise SystemExit("depth labels do not partition the region")

    summary = {
        "region": region,
        "region_bases": region_bases,
        "samples_in_depth_column_order": samples,
        "rule": {
            "depth_source": "samtools depth -a (default flag filter: UNMAP,SECONDARY,QCFAIL,DUP)",
            "mapping_quality_min": mapq,
            "base_quality_min": baseq,
            "min_depth": min_depth,
            "max_depth_factor": max_depth_factor,
            "max_depth": "max_depth_factor x the sample's median depth over this region",
            "sample_level": "min_depth <= depth <= max_depth",
            "cohort_level": f"callable in at least {min_sample_fraction} of the {len(samples)} samples",
            "reference_n": "N bases have depth 0 and are never callable",
            "median_depth_bins": all_depth_labels(),
            "not_a_claim": "callable describes alignment depth and quality, not genotype accuracy",
            "status": "benchmark stratification parameters; not pipeline defaults",
        },
        "bases": bases,
        "samples": {
            name: {
                "region_median_depth": medians[i],
                "max_depth": max_depths[i],
                "callable_bases": sample_callable[i],
                "callable_fraction": round(sample_callable[i] / region_bases, 6),
                "depth_q10_q50_q90": [_quantile(hists[i], q) for q in (0.1, 0.5, 0.9)],
                "zero_depth_bases": hists[i][0],
            }
            for i, name in enumerate(samples)
        },
        "files": {
            f"{label}.bed": hashlib.sha256((output / f"{label}.bed").read_bytes()).hexdigest()
            for label in labels
        },
    }
    (output / "callable_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--depth", required=True, type=Path, help="gzip samtools depth -a output")
    parser.add_argument("--samples", required=True, help="comma-separated, in depth column order")
    parser.add_argument("--region", required=True, help="contig:start-end, 1-based inclusive")
    parser.add_argument("--min-depth", type=int, required=True)
    parser.add_argument("--max-depth-factor", type=float, required=True)
    parser.add_argument("--min-sample-fraction", type=float, required=True)
    parser.add_argument(
        "--mapq", type=int, required=True, help="as passed to samtools depth -Q (recorded)"
    )
    parser.add_argument(
        "--baseq", type=int, required=True, help="as passed to samtools depth -q (recorded)"
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    build(
        depth=args.depth,
        samples=args.samples.split(","),
        region=args.region,
        min_depth=args.min_depth,
        max_depth_factor=args.max_depth_factor,
        min_sample_fraction=args.min_sample_fraction,
        mapq=args.mapq,
        baseq=args.baseq,
        output=args.output_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
