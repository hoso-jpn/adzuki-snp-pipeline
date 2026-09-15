#!/usr/bin/env python3
"""Empirical self-alignment mappability for the exact reference (Issue #65).

No precomputed mask from another organism or assembly is used. Reads are cut
from the reference itself -- length `--read-length`, one every `--step` bases,
error-free -- aligned back with the pipeline's own pinned BWA-MEM2 container,
and each `step`-sized segment is labelled by the read that starts at it:

* `high_mappability`  primary alignment to its own origin with MAPQ >= threshold;
* `low_mappability`   unmapped, placed elsewhere, or MAPQ below the threshold;
* `mappability_undefined` no read starts there (the read would contain N, or
  would run past the contig end).

The three files partition the reference. This is a property of the reference
under these exact read and aligner settings, not of any sample's reads: real
reads have errors, other lengths and pairing, which this does not model.

    mappability.py reads --reference-fasta ref.fna --read-length 150 --step 50 | gzip > reads.fq.gz
    bwa-mem2 mem -t 32 ref.fna reads.fq.gz | samtools view -F 0x900 | \
        mappability.py summarize --reference-fai ref.fna.fai --read-length 150 --step 50 --min-mapq 30 --output-dir out
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from build_region_assets import read_contigs


def reads(fasta: Path, read_length: int, step: int, out) -> None:
    quality = "I" * read_length
    for name, seq in read_contigs(fasta):
        upper = seq.upper()
        for start in range(0, len(upper) - read_length + 1, step):
            chunk = upper[start : start + read_length]
            if "N" in chunk:
                continue
            out.write(f"@{name}:{start}\n{chunk}\n+\n{quality}\n")


def summarize(
    fai: Path, read_length: int, step: int, min_mapq: int, sam, output: Path
) -> dict[str, object]:
    contigs = []
    for line in Path(fai).read_text().splitlines():
        if line:
            name, length = line.split("\t")[:2]
            contigs.append((name, int(length)))
    labels = {
        name: bytearray((length + step - 1) // step) for name, length in contigs
    }  # 0 undefined, 1 high, 2 low
    reads_seen = high = 0
    for line in sam:
        if line.startswith("@"):
            continue
        fields = line.split("\t", 6)
        flag = int(fields[1])
        if flag & 0x900:
            continue
        origin_contig, origin_start = fields[0].rsplit(":", 1)
        origin_start = int(origin_start)
        reads_seen += 1
        mapped = not flag & 0x4
        exact = mapped and fields[2] == origin_contig and int(fields[3]) - 1 == origin_start
        good = exact and int(fields[4]) >= min_mapq
        labels[origin_contig][origin_start // step] = 1 if good else 2
        high += good
    output.mkdir(parents=True, exist_ok=False)
    names = {1: "high_mappability", 2: "low_mappability", 0: "mappability_undefined"}
    bases = dict.fromkeys(names.values(), 0)
    handles = {code: (output / f"{name}.bed").open("w") for code, name in names.items()}
    try:
        for name, length in contigs:
            array = labels[name]
            run_code, run_start = array[0] if array else 0, 0
            for index in range(1, len(array) + 1):
                code = array[index] if index < len(array) else None
                if code != run_code:
                    end = min(index * step, length)
                    handles[run_code].write(f"{name}\t{run_start * step}\t{end}\n")
                    bases[names[run_code]] += end - run_start * step
                    run_code, run_start = code, index
    finally:
        for handle in handles.values():
            handle.close()
    total = sum(length for _, length in contigs)
    if sum(bases.values()) != total:
        raise SystemExit("mappability labels do not partition the reference")
    manifest = {
        "reads_aligned": reads_seen,
        "reads_high": high,
        "parameters": {
            "read_length": read_length,
            "step": step,
            "min_mapq": min_mapq,
            "reads": "error-free, single-end, forward strand, skipped when they contain N",
            "rule": "primary alignment at its own origin with MAPQ >= min_mapq",
        },
        "reference_fai_sha256": hashlib.sha256(Path(fai).read_bytes()).hexdigest(),
        "assets": {
            name: {
                "file": f"{name}.bed",
                "bases": bases[name],
                "fraction_of_reference": round(bases[name] / total, 6),
                "sha256": hashlib.sha256((output / f"{name}.bed").read_bytes()).hexdigest(),
                "kind": "partition:mappability",
            }
            for name in names.values()
        },
    }
    (output / "mappability.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("reads")
    r.add_argument("--reference-fasta", required=True, type=Path)
    r.add_argument("--read-length", type=int, default=150)
    r.add_argument("--step", type=int, default=50)
    s = sub.add_parser("summarize")
    s.add_argument("--reference-fai", required=True, type=Path)
    s.add_argument("--read-length", type=int, default=150)
    s.add_argument("--step", type=int, default=50)
    s.add_argument("--min-mapq", type=int, default=30)
    s.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.command == "reads":
        reads(args.reference_fasta, args.read_length, args.step, sys.stdout)
    else:
        summarize(
            args.reference_fai,
            args.read_length,
            args.step,
            args.min_mapq,
            sys.stdin,
            args.output_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
