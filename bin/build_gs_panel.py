#!/usr/bin/env python3
"""Build a genomic selection (GS) SNP panel from the GS-eligible PASS VCF.

Reads `cohort_gs.snp.pass.vcf.gz` directly (no bcftools re-invocation, to
keep an independent `bcftools query` cross-check meaningful) and emits,
from a single pass over the file:

- a genotype dosage matrix (variant rows x sample columns, matching the
  on-disk convention of the sibling `genomic-prediction-resnet-hybrid`
  repository's own SoyNAM genotype loader, which reads a marker-rows x
  sample-columns TSV and transposes it after loading -- see
  `docs/gs_panel_data_contract.md` for the full reasoning);
- sample metadata (concern 4, separate from the matrix itself);
- variant metadata (concern 4, separate from the matrix itself);
- a genotype-encoding accounting report, broken down by *why* a cell
  became missing (concern 5's genotype-level counterpart -- never
  folding every non-standard call into one undifferentiated bucket).

Dosage encoding is -1 / 0 / +1 for homozygous-reference / heterozygous /
homozygous-alternate, matching the sibling repository's own
`GENOTYPE_ENCODING` convention (`soynam_data.py`): a real SNP/indel
dosage is `allele_count - 1`, so `0/0` -> -1, `0/1` (or `1/0`) -> 0, and
`1/1` -> +1. Missing is IEEE754 NaN, never a sentinel integer, because
a sentinel could collide with a real dosage value or require every
downstream consumer to special-case it.

Phasing (`|` vs `/`) does not affect dosage: per the VCF specification,
the separator only records whether the call is phased, not which or how
many alleles are present, so `0|1` carries exactly the same allele
count as `0/1` and must resolve to the same dosage (0). Earlier
revisions of this script treated any phased call as missing, which was
wrong -- phase is orthogonal to additive dosage. Every genotype is
still checked for how many phased calls it contained
(`phased_genotype_count` in the accounting output), but that count is
informational and never removes a cell from the dosage matrix by
itself.

A genotype is encoded as a dosage as long as it is diploid with a
biallelic-index call (`0/0`, `0/1`, `1/0`, `1/1`, in either phasing).
Every other shape -- missing (`.`, `./.`, or any allele position that
is `.`), non-diploid (haploid, triploid, ...), or an allele index
outside `{0, 1}` (defensive: the input is already biallelic-only by
construction) -- is treated as missing in the matrix, but counted
under its own specific reason so "never silently coerce" is checkable
with real numbers, not just asserted in prose.

This encoding is diploid-only by design (the genotype encoding schema
`diploid_additive_dosage_v1`, which is versioned independently of the
GS panel manifest's own `schema_version`): `--sample-ploidy`
must equal 2, checked before any other work, because a non-diploid
ploidy would make every genotype call "non-diploid-shaped" by
definition, silently producing an all-missing (but successfully
completing) panel rather than a meaningful error.

## Memory (Issue #44)

A real 20-sample cohort produced 9,252,873 variants, and this script is
expected to keep working as cohorts grow. Nothing here may therefore
retain a quantity that scales with `variant_count`.

The production entry point (`main` -> `stream_gs_panel`) makes exactly
one pass over the VCF and, for each data row, immediately classifies
that row's genotypes, writes the matrix row into a streaming gzip
compressor, writes the variant metadata row into an open TSV handle,
and folds the row's genotypes into fixed-size and per-sample counters.
Nothing about the row survives the loop iteration. Peak resident memory
is therefore bounded by:

  * the sample names (`O(sample_count)`);
  * two per-sample counter lists (`O(sample_count)`);
  * a fixed-size set of cohort-wide counters;
  * the single row being processed (`O(sample_count)`);
  * the I/O and compression buffers (fixed, `_MATRIX_FLUSH_THRESHOLD_BYTES`
    plus zlib's own window).

None of those depend on `variant_count`.

The `parse_gs_pass_vcf` / `build_*_rows` / `write_matrix` functions
below are the *reference* implementation: they materialize whole
documents in memory and are `O(variant_count x sample_count)`. They are
retained deliberately -- they express each output's shape declaratively,
and the tests use them as an independent oracle that the streaming
implementation must reproduce exactly -- but the production path must
never call them, which `tests/bin/test_build_gs_panel.py` asserts
directly. Both implementations share `classify_genotype` and the row
formatting helpers, so genotype semantics and output shape cannot drift
between them; what differs is only what is retained.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import zlib
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

DOSAGE_BY_ALT_COUNT: dict[int, str] = {0: "-1", 1: "0", 2: "1"}
MISSING_CELL_TOKEN = "nan"
NOT_APPLICABLE = "NA"

#: VCF's nine fixed columns (CHROM POS ID REF ALT QUAL FILTER INFO FORMAT)
#: before the first sample column.
FIXED_COLUMN_COUNT = 9

#: `gzip.compress(payload, mtime=0)` delegates to zlib's own gzip wrapper
#: (`wbits=31`), which writes `MTIME=0` and `OS=3`. Compressing
#: incrementally through `zlib.compressobj` with the same wbits and the
#: same compression level reproduces that byte-for-byte -- see
#: `_StreamingGzipWriter` and the regression test that pins it.
GZIP_WBITS = 31
GZIP_COMPRESS_LEVEL = 9

#: How much uncompressed matrix text to buffer before handing it to the
#: compressor. Fixed, so it does not scale with the variant count.
_MATRIX_FLUSH_THRESHOLD_BYTES = 1 << 20

GENOTYPE_ACCOUNTING_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")
SAMPLE_METADATA_HEADER: tuple[str, ...] = (
    "cohort_id",
    "sample_index",
    "sample_id",
    "missing_genotype_count",
    "missing_genotype_rate",
    "non_standard_genotype_count",
)
VARIANT_METADATA_HEADER: tuple[str, ...] = (
    "cohort_id",
    "variant_index",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "missing_genotype_count",
    "missing_genotype_rate",
)

#: Accounting metric names, in the order the accounting TSV lists them.
#: Shared by the streaming and reference implementations so the two can
#: never disagree about the document's shape.
_METRIC_BY_DOSAGE: dict[str, str] = {
    "-1": "standard_hom_ref_calls",
    "0": "standard_het_calls",
    "1": "standard_hom_alt_calls",
}
_METRIC_BY_NON_STANDARD_CATEGORY: dict[str, str] = {
    "missing": "missing_calls",
    "non_diploid": "non_diploid_calls_treated_as_missing",
    "non_biallelic_index": "non_biallelic_index_calls_treated_as_missing",
}
_ACCOUNTING_METRIC_ORDER: tuple[str, ...] = (
    "standard_hom_ref_calls",
    "standard_het_calls",
    "standard_hom_alt_calls",
    "missing_calls",
    "non_diploid_calls_treated_as_missing",
    "non_biallelic_index_calls_treated_as_missing",
)


class MalformedVcfError(Exception):
    """Raised when a GS-eligible PASS VCF cannot be turned into a panel safely."""


@dataclass(frozen=True)
class GsPassRecord:
    """The subset of one GS-eligible PASS VCF row this tool needs."""

    chrom: str
    pos: str
    ref: str
    alt: str
    qual: str
    sample_genotypes: tuple[str, ...]

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


@dataclass(frozen=True)
class GsPassVcf:
    """Every sample name and record read from a GS-eligible PASS VCF."""

    sample_names: tuple[str, ...]
    records: tuple[GsPassRecord, ...]


@dataclass(frozen=True)
class GenotypeCell:
    """One sample's classified genotype at one variant.

    ``is_phased`` is orthogonal to ``category``/``dosage``: a phased
    call that is otherwise a clean diploid biallelic-index genotype is
    ``category="standard"`` with a real dosage, exactly like its
    unphased counterpart -- phasing is tracked for informational
    accounting only, never as a reason to treat a cell as missing.
    """

    category: str
    dosage: str
    is_phased: bool


# --------------------------------------------------------------------------
# VCF reading: shared, strictly-validating primitives
#
# Issue #44 review: the previous reader accepted any data row with at
# least 10 tab-separated fields and then indexed into it positionally.
# A row carrying fewer sample columns than `#CHROM` declared was read
# without complaint -- the review demonstrated a two-sample header whose
# data row carried one genotype being parsed cleanly into a matrix row
# one cell short of the matrix header. A row carrying extra columns
# silently gained cells the header never promised. Both reach the matrix
# as a column misalignment, which is what the downstream reconciliation
# exists to catch -- but the builder should not emit it in the first
# place, and a check that only looks for "at least 10" cannot see it.
#
# Every data row must now carry exactly `9 + sample_count` fields, and
# every sample field must actually contain the subfield FORMAT's GT
# index points at: `sample_field.split(":")[gt_index]` previously raised
# a bare `IndexError` that escaped as an unhandled traceback rather than
# a diagnosable failure.
#
# Both the streaming and the reference reader go through these, so the
# two can never disagree about which inputs are well-formed.
# --------------------------------------------------------------------------


def _locate_gt_index(format_field: str, path: Path) -> int:
    keys = format_field.split(":")
    if "GT" not in keys:
        raise MalformedVcfError(f"{path}: FORMAT field has no GT subfield: {format_field}")
    return keys.index("GT")


def _parse_chrom_header(fields: list[str], path: Path) -> tuple[str, ...]:
    """Read the sample names out of a `#CHROM` header line."""
    if len(fields) <= FIXED_COLUMN_COUNT:
        raise MalformedVcfError(
            f"{path}: #CHROM header has {len(fields)} fields, expected "
            f"at least {FIXED_COLUMN_COUNT + 1} ({FIXED_COLUMN_COUNT} fixed "
            "columns plus one or more samples)"
        )
    return tuple(fields[FIXED_COLUMN_COUNT:])


def _extract_row_genotypes(
    fields: list[str],
    sample_names: tuple[str, ...],
    path: Path,
    line_number: int,
) -> tuple[str, ...]:
    """Pull one data row's GT strings out, refusing any shape mismatch.

    The column count is checked against the `#CHROM` header rather than
    against a floor: a row is wrong if it carries *any* number of sample
    columns other than the number the header declared, in either
    direction, no matter how late in the file it appears.
    """
    expected = FIXED_COLUMN_COUNT + len(sample_names)
    if len(fields) != expected:
        raise MalformedVcfError(
            f"{path}: line {line_number}: data row has {len(fields)} "
            f"tab-separated fields, expected exactly {expected} "
            f"({FIXED_COLUMN_COUNT} fixed columns plus {len(sample_names)} "
            "sample columns declared by the #CHROM header)"
        )

    gt_index = _locate_gt_index(fields[8], path)

    genotypes: list[str] = []
    for sample_position, sample_field in enumerate(fields[FIXED_COLUMN_COUNT:]):
        subfields = sample_field.split(":")
        if gt_index >= len(subfields):
            raise MalformedVcfError(
                f"{path}: line {line_number}: sample "
                f"'{sample_names[sample_position]}' has {len(subfields)} "
                f"FORMAT subfield(s), but FORMAT places GT at index {gt_index}"
            )
        genotypes.append(subfields[gt_index])

    return tuple(genotypes)


def classify_genotype(gt: str) -> GenotypeCell:
    """Classify one raw GT string into a category and its matrix dosage token.

    Phasing (``|`` vs ``/``) never changes the category or dosage: it
    only changes ``is_phased``. The VCF specification defines ``|``/``/``
    as recording phase, not allele identity or count, so a phased call
    is resolved exactly like its unphased counterpart.
    """
    is_phased = "|" in gt
    alleles = gt.split("|") if is_phased else gt.split("/")

    if any(allele in (".", "") for allele in alleles):
        return GenotypeCell(category="missing", dosage=MISSING_CELL_TOKEN, is_phased=is_phased)

    if len(alleles) != 2:
        return GenotypeCell(category="non_diploid", dosage=MISSING_CELL_TOKEN, is_phased=is_phased)

    if any(allele not in ("0", "1") for allele in alleles):
        return GenotypeCell(
            category="non_biallelic_index", dosage=MISSING_CELL_TOKEN, is_phased=is_phased
        )

    alt_count = alleles.count("1")
    return GenotypeCell(
        category="standard", dosage=DOSAGE_BY_ALT_COUNT[alt_count], is_phased=is_phased
    )


def _format_rate(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return NOT_APPLICABLE

    return f"{numerator / denominator:.6f}"


def _variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


# --------------------------------------------------------------------------
# Row formatting, shared by the streaming and reference implementations.
#
# Only the *shape* of each output is shared here. The counts fed into
# these functions are derived independently by each implementation --
# the reference one by re-scanning materialized records, the streaming
# one by folding each row in as it is read -- which is what makes the
# equivalence tests between them meaningful rather than tautological.
# --------------------------------------------------------------------------


def _sample_metadata_row(
    cohort_id: str,
    sample_index: int,
    sample_id: str,
    missing_count: int,
    non_standard_count: int,
    variant_count: int,
) -> list[str]:
    return [
        cohort_id,
        str(sample_index),
        sample_id,
        str(missing_count),
        _format_rate(missing_count, variant_count),
        str(non_standard_count),
    ]


def _variant_metadata_row(
    cohort_id: str,
    variant_index: int,
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    qual: str,
    missing_count: int,
    sample_count: int,
) -> list[str]:
    return [
        cohort_id,
        str(variant_index),
        _variant_key(chrom, pos, ref, alt),
        chrom,
        pos,
        ref,
        alt,
        qual,
        str(missing_count),
        _format_rate(missing_count, sample_count),
    ]


def _accounting_rows_from_counts(
    cohort_id: str,
    total_genotype_cells: int,
    counts: dict[str, int],
    phased_genotype_count: int,
) -> list[list[str]]:
    """Lay out the accounting TSV rows in their contractual order."""
    total_treated_as_missing = (
        counts["missing_calls"]
        + counts["non_diploid_calls_treated_as_missing"]
        + counts["non_biallelic_index_calls_treated_as_missing"]
    )

    rows = [[cohort_id, "total_genotype_cells", str(total_genotype_cells)]]
    for metric in _ACCOUNTING_METRIC_ORDER:
        rows.append([cohort_id, metric, str(counts[metric])])
    rows.append([cohort_id, "total_treated_as_missing", str(total_treated_as_missing)])
    rows.append([cohort_id, "phased_genotype_count", str(phased_genotype_count)])

    return rows


def _new_accounting_counts() -> dict[str, int]:
    return dict.fromkeys(_ACCOUNTING_METRIC_ORDER, 0)


def _summary_text_from_accounting(
    cohort_id: str,
    accounting: dict[str, str],
    variant_count: int,
    sample_count: int,
) -> str:
    """Render the human-readable genotype-encoding summary."""
    lines = [
        "GS panel genotype encoding summary",
        f"Cohort ID: {cohort_id}",
        f"Variants: {variant_count}",
        f"Samples: {sample_count}",
        f"Total genotype cells: {accounting['total_genotype_cells']}",
        f"  standard hom-ref (-1): {accounting['standard_hom_ref_calls']}",
        f"  standard het (0): {accounting['standard_het_calls']}",
        f"  standard hom-alt (+1): {accounting['standard_hom_alt_calls']}",
        f"  missing (nan): {accounting['missing_calls']}",
        (
            "  non-diploid, treated as missing (nan): "
            f"{accounting['non_diploid_calls_treated_as_missing']}"
        ),
        (
            "  non-biallelic-index, treated as missing (nan): "
            f"{accounting['non_biallelic_index_calls_treated_as_missing']}"
        ),
        f"Total cells treated as missing: {accounting['total_treated_as_missing']}",
        (
            f"Phased genotype calls: {accounting['phased_genotype_count']} "
            "(informational only -- phasing does not affect dosage or "
            "missingness; a phased call that is otherwise a clean diploid "
            "biallelic-index genotype is encoded exactly like its unphased "
            "counterpart, per the VCF specification's definition of "
            "'|' as recording phase, not allele identity)."
        ),
        (
            "Every non-standard genotype shape is counted under its own "
            "reason rather than a single undifferentiated 'missing' "
            "bucket; see the sample and variant metadata files for the "
            "same breakdown at finer granularity."
        ),
    ]

    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Reference implementation (materializing).
#
# NOT used by the production path -- see this module's docstring. These
# functions hold the whole document in memory, which is precisely what
# Issue #44 removed from the CLI; they survive because they state each
# output's content declaratively, which makes them a good oracle for the
# streaming implementation to be tested against.
# --------------------------------------------------------------------------


def parse_gs_pass_vcf(path: Path) -> GsPassVcf:
    """Parse the CHROM/POS/REF/ALT/QUAL/GT columns of a bgzipped VCF.

    Reference implementation: retains every record, so its memory grows
    with the variant count. `stream_gs_panel` is the production reader.
    """
    sample_names: tuple[str, ...] | None = None
    records: list[GsPassRecord] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#CHROM"):
                if sample_names is not None:
                    raise MalformedVcfError(
                        f"{path}: line {line_number}: a second #CHROM header line"
                    )
                sample_names = _parse_chrom_header(line.split("\t"), path)
                continue

            if line.startswith("#"):
                continue

            if sample_names is None:
                raise MalformedVcfError(f"{path}: data row seen before #CHROM header")

            fields = line.split("\t")
            genotypes = _extract_row_genotypes(fields, sample_names, path, line_number)

            chrom, pos, _id, ref, alt, qual = fields[:6]
            records.append(
                GsPassRecord(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    sample_genotypes=genotypes,
                )
            )

    if sample_names is None:
        raise MalformedVcfError(f"{path}: no #CHROM header line found")

    return GsPassVcf(sample_names=sample_names, records=tuple(records))


def build_matrix_rows(vcf: GsPassVcf) -> list[list[str]]:
    """Build the genotype matrix's data rows: one row per variant."""
    rows: list[list[str]] = []

    for record in vcf.records:
        cells = [classify_genotype(gt) for gt in record.sample_genotypes]
        rows.append([record.variant_key, *(cell.dosage for cell in cells)])

    return rows


def build_sample_metadata_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the sample metadata rows: one row per sample, in header order."""
    total_variants = len(vcf.records)
    missing_counts = [0] * len(vcf.sample_names)
    non_standard_counts = [0] * len(vcf.sample_names)

    for record in vcf.records:
        for sample_position, gt in enumerate(record.sample_genotypes):
            cell = classify_genotype(gt)
            if cell.category == "missing":
                missing_counts[sample_position] += 1
            elif cell.category != "standard":
                missing_counts[sample_position] += 1
                non_standard_counts[sample_position] += 1

    return [
        _sample_metadata_row(
            cohort_id,
            sample_index,
            sample_id,
            missing_counts[sample_index],
            non_standard_counts[sample_index],
            total_variants,
        )
        for sample_index, sample_id in enumerate(vcf.sample_names)
    ]


def build_variant_metadata_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the variant metadata rows: one row per variant, in file order."""
    total_samples = len(vcf.sample_names)
    rows: list[list[str]] = []

    for variant_index, record in enumerate(vcf.records):
        missing_count = sum(
            1
            for gt in record.sample_genotypes
            if classify_genotype(gt).category != "standard"
        )
        rows.append(
            _variant_metadata_row(
                cohort_id,
                variant_index,
                record.chrom,
                record.pos,
                record.ref,
                record.alt,
                record.qual,
                missing_count,
                total_samples,
            )
        )

    return rows


def build_genotype_accounting_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the cohort-wide genotype-encoding accounting rows.

    ``phased_genotype_count`` is reported separately from every other
    metric here: it counts calls that were phased, regardless of
    whether they were standard, missing, or otherwise non-standard, and
    is never added into ``total_treated_as_missing`` -- a phased call
    that resolves to a real dosage is not missing.
    """
    counts = _new_accounting_counts()

    total_genotype_cells = 0
    phased_genotype_count = 0
    for record in vcf.records:
        for gt in record.sample_genotypes:
            total_genotype_cells += 1
            cell = classify_genotype(gt)
            if cell.is_phased:
                phased_genotype_count += 1
            if cell.category == "standard":
                counts[_METRIC_BY_DOSAGE[cell.dosage]] += 1
            else:
                counts[_METRIC_BY_NON_STANDARD_CATEGORY[cell.category]] += 1

    return _accounting_rows_from_counts(
        cohort_id, total_genotype_cells, counts, phased_genotype_count
    )


def build_genotype_accounting_summary_text(cohort_id: str, vcf: GsPassVcf) -> str:
    """Build the human-readable genotype-encoding accounting summary."""
    accounting = {row[1]: row[2] for row in build_genotype_accounting_rows(cohort_id, vcf)}

    return _summary_text_from_accounting(
        cohort_id, accounting, len(vcf.records), len(vcf.sample_names)
    )


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_matrix(path: Path, vcf: GsPassVcf) -> None:
    """Write the gzipped genotype matrix: variant rows x sample columns.

    Reference implementation (see this module's docstring): it builds
    the entire matrix text, encodes it, and compresses it in one shot,
    so its memory grows with the variant count. `stream_gs_panel` writes
    the same bytes incrementally.

    The header is always written, even with zero variants, so an empty
    panel never loses the sample list -- only the data rows are absent.

    Compressed with ``mtime=0`` so that identical logical content always
    produces byte-identical compressed output: ``gzip.open`` embeds the
    current wall-clock time in the gzip header by default, which would
    otherwise make two runs over the same input produce different
    checksums in the manifest even though nothing about the data
    changed -- undermining the "same input reproduces the same panel"
    guarantee this contract requires.
    """
    header = ["variant_key", *vcf.sample_names]
    rows = build_matrix_rows(vcf)

    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    text = "\n".join(lines) + "\n"

    path.write_bytes(gzip.compress(text.encode("utf-8"), mtime=0))


# --------------------------------------------------------------------------
# Streaming implementation: the production path.
# --------------------------------------------------------------------------


class _StreamingGzipWriter:
    """Write gzip bytes incrementally, identically to `gzip.compress(..., mtime=0)`.

    Issue #44 review asked, correctly, that swapping the one-shot
    `gzip.compress` for a streaming writer not be *assumed* to preserve
    the compressed bytes. It does not, for the obvious candidate:
    `gzip.GzipFile` writes its own header with `OS=255` ("unknown"),
    whereas `gzip.compress(payload, mtime=0)` delegates to zlib, whose
    gzip wrapper writes `OS=3`. Same decompressed bytes, different file.

    Going through `zlib.compressobj` with `wbits=31` uses that same zlib
    wrapper, at the same compression level, and deflate's output does
    not depend on how the input was chunked as long as no intermediate
    flush is forced. The result is byte-identical to the old one-shot
    call -- which keeps every previously published matrix checksum
    valid. `test_build_gs_panel.py` pins this against a live
    `gzip.compress` oracle rather than a hardcoded digest, so a future
    Python or zlib that broke the equivalence would fail the suite
    loudly instead of silently changing published checksums.

    Text is buffered up to a fixed threshold before being handed to the
    compressor, so neither the buffer nor the compressor scales with the
    number of variants written.
    """

    def __init__(self, path: Path) -> None:
        self._handle = path.open("wb")
        self._compressor = zlib.compressobj(GZIP_COMPRESS_LEVEL, zlib.DEFLATED, GZIP_WBITS)
        self._pending: list[str] = []
        self._pending_length = 0

    def write(self, text: str) -> None:
        self._pending.append(text)
        self._pending_length += len(text)
        if self._pending_length >= _MATRIX_FLUSH_THRESHOLD_BYTES:
            self._drain()

    def _drain(self) -> None:
        if not self._pending:
            return
        chunk = "".join(self._pending).encode("utf-8")
        self._pending.clear()
        self._pending_length = 0
        compressed = self._compressor.compress(chunk)
        if compressed:
            self._handle.write(compressed)

    def close(self) -> None:
        try:
            self._drain()
            self._handle.write(self._compressor.flush())
            self._handle.flush()
            os.fsync(self._handle.fileno())
        finally:
            self._handle.close()

    def __enter__(self) -> _StreamingGzipWriter:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()


@dataclass
class _PanelTotals:
    """Everything the post-scan outputs need, in bounded space.

    Every field here is either fixed-size or `O(sample_count)`. Nothing
    accumulates per variant.
    """

    sample_names: tuple[str, ...]
    variant_count: int
    total_genotype_cells: int
    phased_genotype_count: int
    counts: dict[str, int]
    sample_missing_counts: list[int]
    sample_non_standard_counts: list[int]


def _write_tsv_header(handle, header: tuple[str, ...]) -> None:
    handle.write("\t".join(header) + "\n")


def _write_tsv_row(handle, row: list[str]) -> None:
    handle.write("\t".join(row) + "\n")


def stream_gs_panel(
    *,
    gs_pass_vcf: Path,
    cohort_id: str,
    matrix_path: Path,
    sample_metadata_path: Path,
    variant_metadata_path: Path,
    genotype_accounting_path: Path,
    genotype_accounting_summary_path: Path,
) -> _PanelTotals:
    """Build all five GS panel outputs in one bounded-memory pass.

    The five paths are written directly; callers that need the outputs
    to appear only on success (`main` does) pass staging paths here and
    publish them afterwards.

    Each genotype is classified exactly once. The resulting cell feeds
    the matrix token, the variant's missing count, that sample's
    counters, and the cohort-wide accounting -- rather than being
    re-derived per output, as the reference implementation's five
    separate scans do.
    """
    sample_names: tuple[str, ...] | None = None
    variant_count = 0
    total_genotype_cells = 0
    phased_genotype_count = 0
    counts = _new_accounting_counts()
    sample_missing_counts: list[int] = []
    sample_non_standard_counts: list[int] = []

    with ExitStack() as stack:
        handle = stack.enter_context(gzip.open(gs_pass_vcf, "rt", encoding="utf-8"))
        matrix = stack.enter_context(_StreamingGzipWriter(matrix_path))
        variant_metadata = stack.enter_context(
            variant_metadata_path.open("w", encoding="utf-8")
        )
        _write_tsv_header(variant_metadata, VARIANT_METADATA_HEADER)

        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#CHROM"):
                if sample_names is not None:
                    raise MalformedVcfError(
                        f"{gs_pass_vcf}: line {line_number}: a second #CHROM header line"
                    )
                sample_names = _parse_chrom_header(line.split("\t"), gs_pass_vcf)
                sample_missing_counts = [0] * len(sample_names)
                sample_non_standard_counts = [0] * len(sample_names)
                # The matrix header is written as soon as the sample list
                # is known, so a zero-variant panel still publishes it.
                matrix.write("\t".join(["variant_key", *sample_names]) + "\n")
                continue

            if line.startswith("#"):
                continue

            if sample_names is None:
                raise MalformedVcfError(
                    f"{gs_pass_vcf}: data row seen before #CHROM header"
                )

            fields = line.split("\t")
            genotypes = _extract_row_genotypes(
                fields, sample_names, gs_pass_vcf, line_number
            )
            chrom, pos, _id, ref, alt, qual = fields[:6]

            dosages: list[str] = []
            row_missing_count = 0
            for sample_position, gt in enumerate(genotypes):
                cell = classify_genotype(gt)
                dosages.append(cell.dosage)

                total_genotype_cells += 1
                if cell.is_phased:
                    phased_genotype_count += 1

                if cell.category == "standard":
                    counts[_METRIC_BY_DOSAGE[cell.dosage]] += 1
                    continue

                counts[_METRIC_BY_NON_STANDARD_CATEGORY[cell.category]] += 1
                row_missing_count += 1
                sample_missing_counts[sample_position] += 1
                if cell.category != "missing":
                    sample_non_standard_counts[sample_position] += 1

            matrix.write(
                "\t".join([_variant_key(chrom, pos, ref, alt), *dosages]) + "\n"
            )
            _write_tsv_row(
                variant_metadata,
                _variant_metadata_row(
                    cohort_id,
                    variant_count,
                    chrom,
                    pos,
                    ref,
                    alt,
                    qual,
                    row_missing_count,
                    len(sample_names),
                ),
            )
            variant_count += 1

    if sample_names is None:
        raise MalformedVcfError(f"{gs_pass_vcf}: no #CHROM header line found")

    totals = _PanelTotals(
        sample_names=sample_names,
        variant_count=variant_count,
        total_genotype_cells=total_genotype_cells,
        phased_genotype_count=phased_genotype_count,
        counts=counts,
        sample_missing_counts=sample_missing_counts,
        sample_non_standard_counts=sample_non_standard_counts,
    )

    _write_post_scan_outputs(
        cohort_id=cohort_id,
        totals=totals,
        sample_metadata_path=sample_metadata_path,
        genotype_accounting_path=genotype_accounting_path,
        genotype_accounting_summary_path=genotype_accounting_summary_path,
    )

    return totals


def _write_post_scan_outputs(
    *,
    cohort_id: str,
    totals: _PanelTotals,
    sample_metadata_path: Path,
    genotype_accounting_path: Path,
    genotype_accounting_summary_path: Path,
) -> None:
    """Write the three outputs that only need the accumulated counters.

    All three are `O(sample_count)` or fixed size, so materializing them
    is bounded regardless of how many variants were read.
    """
    sample_rows = [
        _sample_metadata_row(
            cohort_id,
            sample_index,
            sample_id,
            totals.sample_missing_counts[sample_index],
            totals.sample_non_standard_counts[sample_index],
            totals.variant_count,
        )
        for sample_index, sample_id in enumerate(totals.sample_names)
    ]
    write_tsv(sample_metadata_path, SAMPLE_METADATA_HEADER, sample_rows)

    accounting_rows = _accounting_rows_from_counts(
        cohort_id,
        totals.total_genotype_cells,
        totals.counts,
        totals.phased_genotype_count,
    )
    write_tsv(genotype_accounting_path, GENOTYPE_ACCOUNTING_HEADER, accounting_rows)

    accounting = {row[1]: row[2] for row in accounting_rows}
    genotype_accounting_summary_path.write_text(
        _summary_text_from_accounting(
            cohort_id, accounting, totals.variant_count, len(totals.sample_names)
        ),
        encoding="utf-8",
    )


# --------------------------------------------------------------------------
# Publication.
#
# What this does and does not guarantee (Issue #44 review asked that the
# difference be stated rather than glossed):
#
#   * Every output is built into a staging file beside its final path,
#     so nothing appears at a final path until the whole VCF has been
#     read and all five documents have been produced. A malformed row in
#     the last line of the input therefore leaves no final output, and
#     an existing panel from a previous successful run is untouched.
#   * Each individual `os.replace` is atomic, and a failure partway
#     through the sequence rolls the already-replaced files back from
#     the copies moved aside a moment earlier.
#   * It is NOT a filesystem transaction over five files. `os.replace`
#     is atomic per file; there is no primitive that commits five of
#     them together. A crash (SIGKILL, power loss) in the middle of the
#     replace sequence, or a failure during the rollback itself, can
#     leave a mix of new and old files. The rollback narrows that window
#     to the sequence itself rather than the whole build, which is the
#     honest description of the guarantee.
# --------------------------------------------------------------------------


def _staging_path(final_path: Path) -> Path:
    return final_path.with_name(f".{final_path.name}.partial")


def _rollback_path(final_path: Path) -> Path:
    return final_path.with_name(f".{final_path.name}.previous")


def publish_outputs(final_paths: list[Path]) -> None:
    """Move each staged output onto its final path, rolling back on failure."""
    moved_aside: list[tuple[Path, Path]] = []
    replaced: list[Path] = []

    try:
        for final_path in final_paths:
            rollback_path = _rollback_path(final_path)
            if final_path.exists():
                os.replace(final_path, rollback_path)
                moved_aside.append((final_path, rollback_path))
            os.replace(_staging_path(final_path), final_path)
            replaced.append(final_path)
    except OSError:
        for final_path in replaced:
            final_path.unlink(missing_ok=True)
        for final_path, rollback_path in moved_aside:
            try:
                os.replace(rollback_path, final_path)
            except OSError:
                pass
        raise
    else:
        for _final_path, rollback_path in moved_aside:
            rollback_path.unlink(missing_ok=True)


def discard_staged_outputs(final_paths: list[Path]) -> None:
    """Remove any staging files left behind by a failed build."""
    for final_path in final_paths:
        _staging_path(final_path).unlink(missing_ok=True)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the GS panel builder CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a GS genotype matrix, sample/variant metadata, and a "
            "genotype-encoding accounting report from the GS-eligible PASS VCF."
        )
    )
    parser.add_argument(
        "--gs-pass-vcf",
        required=True,
        type=Path,
        help="Path to cohort_gs.snp.pass.vcf.gz.",
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--sample-ploidy",
        required=True,
        type=int,
        help="The pipeline's configured sample ploidy (params.sample_ploidy); "
        "this schema is diploid-only, so any value other than 2 is a hard error.",
    )
    parser.add_argument(
        "--matrix-output",
        required=True,
        type=Path,
        help="Output path for the gzipped genotype matrix TSV.",
    )
    parser.add_argument(
        "--sample-metadata-output",
        required=True,
        type=Path,
        help="Output path for the sample metadata TSV.",
    )
    parser.add_argument(
        "--variant-metadata-output",
        required=True,
        type=Path,
        help="Output path for the variant metadata TSV.",
    )
    parser.add_argument(
        "--genotype-accounting-output",
        required=True,
        type=Path,
        help="Output path for the genotype-encoding accounting TSV.",
    )
    parser.add_argument(
        "--genotype-accounting-summary-output",
        required=True,
        type=Path,
        help="Output path for the human-readable genotype-encoding summary.",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code.

    Builds every output into a staging file first and publishes the set
    only after the whole VCF has been read successfully, so a malformed
    row anywhere in the input -- including its last line -- leaves no
    final output and does not disturb a previous run's panel.
    """
    args = parse_args(argv)

    if args.sample_ploidy != 2:
        print(
            "build_gs_panel.py: error: this GS panel genotype encoding "
            "(diploid_additive_dosage_v1) is diploid-only, "
            f"but --sample-ploidy was {args.sample_ploidy}. Every genotype call "
            "would be classified as non-diploid-shaped and encoded as missing, "
            "which would silently produce an all-missing panel rather than a "
            "meaningful error. A generalized encoding is tracked as future work; "
            "until then, this schema cannot be used for non-diploid cohorts.",
            file=sys.stderr,
        )
        return 1

    # Order matters only in that it is the order the outputs are
    # published in; every one of them is fully built before any is moved.
    final_paths = [
        args.matrix_output,
        args.sample_metadata_output,
        args.variant_metadata_output,
        args.genotype_accounting_output,
        args.genotype_accounting_summary_output,
    ]

    # `finally`, not an exception handler: the malformed-input paths below
    # leave by `return 1`, which no `except` clause would see, and those
    # are exactly the runs that leave half-written staging files behind.
    # On success the staging files have already been consumed by the
    # publish, so the sweep is a no-op.
    try:
        try:
            stream_gs_panel(
                gs_pass_vcf=args.gs_pass_vcf,
                cohort_id=args.cohort_id,
                matrix_path=_staging_path(args.matrix_output),
                sample_metadata_path=_staging_path(args.sample_metadata_output),
                variant_metadata_path=_staging_path(args.variant_metadata_output),
                genotype_accounting_path=_staging_path(args.genotype_accounting_output),
                genotype_accounting_summary_path=_staging_path(
                    args.genotype_accounting_summary_output
                ),
            )
        except OSError as error:
            print(
                f"build_gs_panel.py: error: cannot read {args.gs_pass_vcf}: {error}",
                file=sys.stderr,
            )
            return 1
        except MalformedVcfError as error:
            print(f"build_gs_panel.py: error: {error}", file=sys.stderr)
            return 1

        publish_outputs(final_paths)
    finally:
        discard_staged_outputs(final_paths)

    return 0


if __name__ == "__main__":
    sys.exit(main())
