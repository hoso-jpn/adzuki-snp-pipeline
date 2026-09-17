#!/usr/bin/env python3
"""Write the Issue #65 synthetic truth fixture (reference, truth, query, regions).

The fixture is small enough to account for by hand. Each record exists to
exercise one accounting rule, named in the comment beside it. The expected
counts live in the unit test, derived from these comments, not from running
the engine -- so the test checks the engine rather than restating it.

Run from the repository root:
    python3 benchmarks/issue65/make_synthetic_fixture.py tests/bin/fixtures/issue65
"""

from __future__ import annotations

import random
import sys
from pathlib import Path


def reference_sequences() -> dict[str, str]:
    rng = random.Random(65)
    chr_t = [rng.choice("ACGT") for _ in range(300)]
    # Homopolymer used by the left-alignment case: C AAAAA G at 99..105 (1-based).
    for pos, base in zip(range(99, 106), "CAAAAAG"):
        chr_t[pos - 1] = base
    for pos, base in ((121, "G"), (140, "A"), (160, "C"), (180, "G")):
        chr_t[pos - 1] = base
    # Boundary deletion 238..241: make its last base differ from its first so
    # normalization does not move it.
    chr_t[237], chr_t[240] = "A", "C"
    chr_u = [rng.choice("ACGT") for _ in range(100)]
    return {"chrT": "".join(chr_t), "chrU": "".join(chr_u)}


HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chrT,length=300>",
    "##contig=<ID=chrU,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def base(seq: dict[str, str], contig: str, pos: int, length: int = 1) -> str:
    return seq[contig][pos - 1 : pos - 1 + length]


def write(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    seq = reference_sequences()
    fasta = "".join(
        f">{name}\n" + "".join(s[i : i + 60] + "\n" for i in range(0, len(s), 60))
        for name, s in seq.items()
    )
    (directory / "reference.fa").write_text(fasta)
    offset = 0
    fai = []
    for name, s in seq.items():
        offset += len(f">{name}\n")
        fai.append(f"{name}\t{len(s)}\t{offset}\t60\t61\n")
        offset += len(s) + (len(s) + 59) // 60
    (directory / "reference.fa.fai").write_text("".join(fai))

    def snp(contig: str, pos: int) -> tuple[str, str]:
        ref = base(seq, contig, pos)
        return ref, {"A": "C", "C": "G", "G": "T", "T": "A"}[ref]

    truth, query = [], []

    def row(target, contig, pos, ref, alt, gt):
        target.append((contig, pos, ref, alt, gt))

    def snps(contig: str, pos: int, truth_gt: str | None, query_gt: str | None) -> None:
        ref, alt = snp(contig, pos)
        if truth_gt is not None:
            row(truth, contig, pos, ref, alt, truth_gt)
        if query_gt is not None:
            row(query, contig, pos, ref, alt, query_gt)

    snps("chrT", 10, "0/1", "0/1")  # TP, GT concordant
    snps("chrT", 20, "1/1", "0/1")  # TP, GT discordant
    snps("chrT", 30, None, "0/1")  # FP (truth has no record)
    snps("chrT", 40, "0/1", None)  # FN (query has no record)
    snps("chrT", 50, "0/1", "./.")  # query no-call
    snps("chrT", 60, "./.", "0/1")  # truth no-call
    snps("chrT", 70, "0/0", "0/1")  # FP (truth hom-ref)
    snps("chrT", 80, "0/1", "0/0")  # FN (query hom-ref)
    snps("chrT", 90, "0/1", "0/.")  # query partial no-call
    row(truth, "chrT", 99, "CA", "C", "0/1")  # deletion, left-aligned
    row(
        query, "chrT", 103, "AA", "A", "0/1"
    )  # same deletion, right-shifted: TP after normalization
    ins = base(seq, "chrT", 121)
    row(truth, "chrT", 121, ins, ins + "TT", "1/1")  # insertion TP
    row(query, "chrT", 121, ins, ins + "TT", "1/1")
    row(truth, "chrT", 140, "A", "C,T", "1/2")  # multi-allelic truth ...
    row(query, "chrT", 140, "A", "C", "0/1")  # ... split in query: 2 TP
    row(query, "chrT", 140, "A", "T", "0/1")
    row(truth, "chrT", 160, "C", "A", "1/1")  # multi-allelic query: 1 TP, G not positive
    row(query, "chrT", 160, "C", "A,G", "1/1")
    row(query, "chrT", 180, "T", "C", "0/1")  # REF mismatch (reference is G): excluded
    row(truth, "chrT", 190, base(seq, "chrT", 190), "<DEL>", "0/1")  # symbolic: excluded
    # Multi-allelic records that mix a sequence ALT with a symbolic one: only the
    # symbolic allele is excluded, and an excluded unit both sides carry counts once.
    ref200, alt200 = snp("chrT", 200)
    row(truth, "chrT", 200, ref200, alt200, "0/1")  # query calls the sequence ALT: TP
    row(query, "chrT", 200, ref200, alt200 + ",*", "0/1")
    ref205, alt205 = snp("chrT", 205)
    row(truth, "chrT", 205, ref205, alt205, "0/1")  # query calls only the symbolic ALT:
    row(query, "chrT", 205, ref205, alt205 + ",*", "0/2")  # FN, plus one symbolic exclusion
    ref210, alt210 = snp("chrT", 210)
    row(truth, "chrT", 210, ref210, alt210 + ",*", "1/2")  # truth GT 1/2 over both alleles:
    row(query, "chrT", 210, ref210, alt210, "0/1")  # TP on the sequence ALT
    ref215, alt215 = snp("chrT", 215)
    row(truth, "chrT", 215, ref215, alt215 + ",*", "1/2")  # both sides carry the same `*`:
    row(query, "chrT", 215, ref215, alt215 + ",*", "1/2")  # TP, and one symbolic exclusion
    wrong = {"A": "C", "C": "G", "G": "T", "T": "A"}[base(seq, "chrT", 220)]
    other = {"A": "C", "C": "G", "G": "T", "T": "A"}[wrong]
    row(truth, "chrT", 220, wrong, other, "0/1")  # REF mismatch on both sides: one exclusion
    row(query, "chrT", 220, wrong, other, "0/1")
    row(truth, "chrZ", 5, "A", "C", "0/1")  # contig not in reference, both sides: one exclusion
    row(query, "chrZ", 5, "A", "C", "0/1")
    deleted = base(seq, "chrT", 238, 4)
    row(truth, "chrT", 238, deleted, deleted[0], "0/1")  # spans region end: excluded
    row(query, "chrT", 238, deleted, deleted[0], "0/1")
    snps("chrT", 250, "0/1", "0/1")  # outside region: excluded
    snps("chrU", 10, "0/1", "0/1")  # second contig TP

    order = {"chrT": 0, "chrU": 1, "chrZ": 2}
    for name, rows, sample in (("truth.vcf", truth, "TRUTH"), ("query.vcf", query, "QUERY")):
        lines = HEADER + [
            "\t".join(
                ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
            )
        ]
        for contig, pos, ref, alt, gt in sorted(rows, key=lambda x: (order[x[0]], x[1], x[3])):
            lines.append("\t".join([contig, str(pos), ".", ref, alt, "50", "PASS", ".", "GT", gt]))
        (directory / name).write_text("\n".join(lines) + "\n")

    (directory / "evaluation_region.bed").write_text("chrT\t0\t240\nchrU\t0\t100\n")
    (directory / "repeat.bed").write_text("chrT\t95\t110\n")
    (directory / "low_mappability.bed").write_text("chrT\t0\t45\n")
    (directory / "gc_low.bed").write_text("chrT\t0\t120\nchrU\t0\t100\n")
    (directory / "gc_high.bed").write_text("chrT\t120\t240\n")


if __name__ == "__main__":
    write(Path(sys.argv[1]))
