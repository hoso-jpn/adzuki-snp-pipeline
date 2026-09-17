"""VCF reading and allele normalization for Issue #65 comparisons.

Comparisons never match on raw POS. Every ALT allele is split into its own
biallelic record, validated against the reference, trimmed and left-aligned
(the parsimony + left-alignment rule of Tan, Abecasis & Kang 2015, which is
what `bcftools norm -f ref -m -any` implements), and only then keyed by
`(contig, pos, ref, alt)`. The unit tests pin this against records the pinned
bcftools container normalized, so the Python rule cannot quietly diverge.

Records that cannot be compared are classified, not dropped. Two of the
reasons belong to the site and one to the single ALT allele, because the
comparison unit is one ALT allele after multi-allelic splitting:

* `coordinate_mismatch` (site) -- the contig is not in the reference, or the
  REF span runs past the contig end;
* `ref_mismatch` (site) -- REF does not equal the reference bases at POS;
* `symbolic_allele` (allele) -- `<DEL>`, `*`, breakends and similar, which this
  sequence-level comparison cannot place. A record such as `ALT=G,*` keeps its
  sequence allele `G`: only the `*` allele is excluded.

Malformed input (a data row with the wrong column count, a non-integer POS, a
GT that is not a genotype) is a hard error: evaluating a broken file and
reporting numbers from it would be worse than failing.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

_BASES = re.compile(r"^[ACGTNacgtn]+$")


class MalformedVcfError(ValueError):
    """The VCF cannot be read safely."""


class Reference:
    """Random access to a FASTA through its .fai, without loading it whole."""

    def __init__(self, fasta: Path, fai: Path | None = None) -> None:
        self.fasta = Path(fasta)
        fai = Path(fai) if fai else Path(str(fasta) + ".fai")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        self.order: list[tuple[str, int]] = []
        for line in fai.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(bases), int(width))
            self.order.append((name, int(length)))
        self._handle = self.fasta.open("rb")

    def length(self, contig: str) -> int | None:
        entry = self.index.get(contig)
        return entry[0] if entry else None

    def fetch(self, contig: str, start0: int, end0: int) -> str:
        """Upper-cased bases `[start0, end0)`; the caller checks bounds first."""
        length, offset, bases, width = self.index[contig]
        if not (0 <= start0 <= end0 <= length):
            raise ValueError(f"{contig}:{start0}-{end0} is outside 0..{length}")
        out = []
        position = start0
        while position < end0:
            line, column = divmod(position, bases)
            take = min(bases - column, end0 - position)
            self._handle.seek(offset + line * width + column)
            out.append(self._handle.read(take).decode("ascii"))
            position += take
        return "".join(out).upper()

    def close(self) -> None:
        self._handle.close()


@dataclass(frozen=True)
class NormalizedAllele:
    contig: str
    pos: int
    ref: str
    alt: str
    alt_index: int  # 1-based index of this ALT in the original record
    changed_by_normalization: bool

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.contig, self.pos, self.ref, self.alt)

    @property
    def span(self) -> int:
        return len(self.ref)

    @property
    def variant_type(self) -> str:
        return "snp" if len(self.ref) == 1 and len(self.alt) == 1 else "indel"


def classify_record(reference: Reference, contig: str, pos: int, ref: str) -> str | None:
    """Site-level exclusion reason for a record, or None if its site is usable.

    This judges only what the whole record shares -- its coordinates and its REF
    -- so a record with both a sequence and a symbolic ALT is not discarded
    wholesale. Each ALT is judged separately by `is_symbolic_allele`.
    """
    length = reference.length(contig)
    if length is None or pos < 1 or pos - 1 + len(ref) > length:
        return "coordinate_mismatch"
    if reference.fetch(contig, pos - 1, pos - 1 + len(ref)) != ref.upper():
        return "ref_mismatch"
    return None


def is_symbolic_allele(alt: str) -> bool:
    """True for `<DEL>`, `*`, breakends and anything else that is not plain bases."""
    return not _BASES.match(alt)


def normalize_allele(
    reference: Reference, contig: str, pos: int, ref: str, alt: str, alt_index: int
) -> NormalizedAllele:
    original = (pos, ref.upper(), alt.upper())
    ref, alt = ref.upper(), alt.upper()
    changed = True
    while changed:
        changed = False
        if ref and alt and ref[-1] == alt[-1] and not (len(ref) == 1 and len(alt) == 1):
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not ref or not alt:
            if pos <= 1:
                break
            pos -= 1
            base = reference.fetch(contig, pos - 1, pos)
            ref, alt = base + ref, base + alt
            changed = True
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return NormalizedAllele(contig, pos, ref, alt, alt_index, (pos, ref, alt) != original)


@dataclass(frozen=True)
class VcfRecord:
    contig: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    genotypes: dict[str, tuple[int | None, ...]]
    line_number: int


def _parse_gt(value: str, where: str) -> tuple[int | None, ...]:
    if value in (".", ""):
        return (None,)
    alleles = re.split(r"[/|]", value)
    parsed: list[int | None] = []
    for allele in alleles:
        if allele == ".":
            parsed.append(None)
        elif allele.isdigit():
            parsed.append(int(allele))
        else:
            raise MalformedVcfError(f"{where}: GT {value!r} is not a genotype")
    return tuple(parsed)


def read_vcf(path: Path, samples: list[str] | None = None):
    """Yield (sample_names, VcfRecord) pairs, validating every row."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header: list[str] | None = None
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            fields = line.split("\t")
            if line.startswith("#CHROM"):
                if header is not None:
                    raise MalformedVcfError(f"{path}: line {number}: a second #CHROM line")
                if len(fields) < 10:
                    raise MalformedVcfError(f"{path}: #CHROM has no sample column")
                header = fields[9:]
                if samples is not None:
                    missing = [s for s in samples if s not in header]
                    if missing:
                        raise MalformedVcfError(f"{path}: samples {missing} are not in the VCF")
                continue
            if header is None:
                raise MalformedVcfError(f"{path}: line {number}: data before #CHROM")
            where = f"{path}: line {number}"
            if len(fields) != 9 + len(header):
                raise MalformedVcfError(
                    f"{where}: {len(fields)} columns, expected {9 + len(header)}"
                )
            if not fields[1].isdigit():
                raise MalformedVcfError(f"{where}: POS {fields[1]!r} is not an integer")
            keys = fields[8].split(":")
            if "GT" not in keys:
                raise MalformedVcfError(f"{where}: FORMAT has no GT")
            gt_index = keys.index("GT")
            wanted = samples if samples is not None else header
            genotypes = {}
            for name in wanted:
                values = fields[9 + header.index(name)].split(":")
                genotypes[name] = _parse_gt(
                    values[gt_index] if gt_index < len(values) else ".", where
                )
            yield (
                header,
                VcfRecord(
                    contig=fields[0],
                    pos=int(fields[1]),
                    ref=fields[3],
                    alts=tuple(fields[4].split(",")) if fields[4] != "." else (),
                    genotypes=genotypes,
                    line_number=number,
                ),
            )
