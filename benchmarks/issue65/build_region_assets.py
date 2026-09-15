#!/usr/bin/env python3
"""Derive sequence-based region assets from the reference FASTA itself (Issue #65).

Everything here is computed from the exact reference the pipeline aligns to,
so coordinates cannot drift from another assembly version:

* `n_regions.bed`         runs of N (assembly gaps and any other unknown bases);
* `repeat_windowmasker.bed` lower-case runs. NCBI documents that eukaryotic
  RefSeq/GenBank genomic FASTA marks repeats identified by WindowMasker in
  lower case (genomes FTP README, "Repetitive sequences"), so for this RefSeq
  assembly these runs are NCBI's WindowMasker repeat annotation carried in the
  same file -- not a mask borrowed from another organism.
* `homopolymer_ge{N}.bed` single-base runs of at least N (upper-cased), a
  sequence-derived low-complexity tag;
* `gc_*.bed`              a partition of the reference into fixed windows by GC
  fraction of their non-N bases, with windows that are mostly N set apart as
  `gc_undefined`.

When NCBI's `genomic_gaps.txt` is given, every listed gap must be an N run in
the FASTA, or the tool fails: that check is what shows the gap annotation and
the sequence are the same coordinate system.

Memory is one contig's sequence at a time (the longest is 65 Mb here).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

GC_EDGES = (0.25, 0.30, 0.35, 0.40, 0.45)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_contigs(fasta: Path):
    name, parts = None, []
    with Path(fasta).open(encoding="ascii") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def gc_label(fraction: float | None) -> str:
    if fraction is None:
        return "gc_undefined"
    lower = 0.0
    for edge in GC_EDGES:
        if fraction < edge:
            return f"gc_{int(lower * 100):02d}_{int(edge * 100):02d}"
        lower = edge
    return f"gc_{int(lower * 100):02d}_100"


def all_gc_labels() -> list[str]:
    labels, lower = [], 0.0
    for edge in GC_EDGES:
        labels.append(f"gc_{int(lower * 100):02d}_{int(edge * 100):02d}")
        lower = edge
    return [*labels, f"gc_{int(lower * 100):02d}_100", "gc_undefined"]


def build(
    fasta: Path, output: Path, window: int, homopolymer_min: int, gaps: Path | None
) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=False)
    n_re = re.compile(r"[Nn]+")
    lower_re = re.compile(r"[acgtn]+")
    homo_re = re.compile(r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % ((homopolymer_min,) * 4))
    handles = {
        "n_regions": (output / "n_regions.bed").open("w"),
        "repeat_windowmasker": (output / "repeat_windowmasker.bed").open("w"),
        f"homopolymer_ge{homopolymer_min}": (output / f"homopolymer_ge{homopolymer_min}.bed").open(
            "w"
        ),
        **{label: (output / f"{label}.bed").open("w") for label in all_gc_labels()},
    }
    bases = dict.fromkeys(handles, 0)
    contigs: list[tuple[str, int]] = []
    n_runs: dict[str, list[tuple[int, int]]] = {}
    try:
        for name, seq in read_contigs(fasta):
            contigs.append((name, len(seq)))
            upper = seq.upper()
            for key, regex, text in (
                ("n_regions", n_re, seq),
                ("repeat_windowmasker", lower_re, seq),
                (f"homopolymer_ge{homopolymer_min}", homo_re, upper),
            ):
                for match in regex.finditer(text):
                    handles[key].write(f"{name}\t{match.start()}\t{match.end()}\n")
                    bases[key] += match.end() - match.start()
                    if key == "n_regions":
                        n_runs.setdefault(name, []).append((match.start(), match.end()))
            runs: list[list] = []
            for start in range(0, len(seq), window):
                chunk = upper[start : start + window]
                n = chunk.count("N")
                defined = len(chunk) - n
                fraction = (
                    None
                    if defined < len(chunk) / 2
                    else (chunk.count("G") + chunk.count("C")) / defined
                )
                label = gc_label(fraction)
                end = start + len(chunk)
                if runs and runs[-1][0] == label and runs[-1][2] == start:
                    runs[-1][2] = end
                else:
                    runs.append([label, start, end])
            for label, start, end in runs:
                handles[label].write(f"{name}\t{start}\t{end}\n")
                bases[label] += end - start
    finally:
        for handle in handles.values():
            handle.close()

    gap_check = None
    if gaps is not None:
        opener = gzip.open if str(gaps).endswith(".gz") else open
        listed = unmatched = 0
        with opener(gaps, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                contig, start1, stop1 = line.split("\t")[:3]
                listed += 1
                start0, end0 = int(start1) - 1, int(stop1)
                if not any(s <= start0 and end0 <= e for s, e in n_runs.get(contig, [])):
                    unmatched += 1
        if unmatched:
            raise SystemExit(f"{unmatched} of {listed} NCBI gaps are not N runs in the FASTA")
        gap_check = {"ncbi_gaps_listed": listed, "all_inside_fasta_n_runs": True}

    genome = sum(length for _, length in contigs)
    gc_total = sum(bases[label] for label in all_gc_labels())
    if gc_total != genome:
        raise SystemExit("GC windows do not partition the reference")
    manifest = {
        "reference_fasta_sha256": sha256(fasta),
        "reference_bases": genome,
        "contigs": len(contigs),
        "parameters": {
            "gc_window_bp": window,
            "gc_edges": list(GC_EDGES),
            "homopolymer_min": homopolymer_min,
            "gc_undefined_rule": "more than half of the window is N",
        },
        "assets": {
            key: {
                "file": f"{key}.bed",
                "sha256": sha256(output / f"{key}.bed"),
                "bases": bases[key],
                "fraction_of_reference": round(bases[key] / genome, 6),
                "kind": "partition:gc" if key.startswith("gc_") else "tag",
            }
            for key in handles
        },
        "ncbi_gap_cross_check": gap_check,
        "sources": {
            "repeat_windowmasker": "lower-case runs of the RefSeq genomic FASTA; NCBI genomes FTP README: "
            "'Repetitive sequences in eukaryotic genome assembly sequence files, as identified by "
            "WindowMasker, have been masked to lower-case.'",
        },
    }
    (output / "region_assets.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gc-window", type=int, default=1000)
    parser.add_argument("--homopolymer-min", type=int, default=10)
    parser.add_argument("--ncbi-gaps", type=Path)
    args = parser.parse_args(argv)
    build(
        args.reference_fasta, args.output_dir, args.gc_window, args.homopolymer_min, args.ncbi_gaps
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
