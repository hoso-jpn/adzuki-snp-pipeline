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
import json
import os
import struct
import sys
import zlib
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path

# Issue #64: the one genotype quality evaluator the matrix, the masked VCF
# and the verifier all share. See bin/gs_genotype_quality.py.
import gs_genotype_quality as quality

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

#: Issue #64: metadata headers when a genotype quality policy is enabled.
#: The existing columns keep their names, positions and meaning -- a masked
#: call is a `nan` cell, so it is part of `missing_genotype_count` -- and one
#: column is appended that says how many of those cells the policy masked.
#: With the policy disabled the historical headers above are written as-is.
QUALITY_MASKED_COLUMN = "quality_masked_genotype_count"
SAMPLE_METADATA_HEADER_WITH_QUALITY: tuple[str, ...] = (
    *SAMPLE_METADATA_HEADER,
    QUALITY_MASKED_COLUMN,
)
VARIANT_METADATA_HEADER_WITH_QUALITY: tuple[str, ...] = (
    *VARIANT_METADATA_HEADER,
    QUALITY_MASKED_COLUMN,
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
    #: Issue #64: the row's FORMAT and whole sample fields, so the reference
    #: implementation can evaluate a genotype quality policy too. Defaults
    #: keep every pre-existing construction of this record valid.
    format_field: str = "GT"
    sample_fields: tuple[str, ...] = ()

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


@dataclass
class _QualityCounts:
    """Issue #64: every genotype quality count, in bounded space.

    Fixed size for a given policy -- at most two fields, six statuses each,
    and fifteen reason pairs -- plus nothing per variant. Shared by the
    streaming and reference implementations so the accounting shape cannot
    differ between them, while each still derives its counts on its own.
    """

    policy: quality.GenotypeQualityPolicy
    evaluated: int = 0
    passed: int = 0
    kept_unevaluated: int = 0
    masked: int = 0
    masked_by_dosage: dict[str, int] = field(
        default_factory=lambda: dict.fromkeys(("-1", "0", "1"), 0)
    )
    status_counts: list[dict[str, int]] = field(init=False)
    reason_counts: dict[tuple[str, ...], int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.status_counts = [
            dict.fromkeys(quality.FIELD_STATUSES, 0) for _ in self.policy.thresholds
        ]

    def record(self, dosage: str, verdict: quality.CallQuality) -> None:
        self.evaluated += 1
        for counts, status in zip(self.status_counts, verdict.statuses):
            counts[status] += 1
        if verdict.masked:
            self.masked += 1
            self.masked_by_dosage[dosage] += 1
            self.reason_counts[verdict.reasons] = self.reason_counts.get(verdict.reasons, 0) + 1
        elif verdict.unevaluated:
            self.kept_unevaluated += 1
        else:
            self.passed += 1

    def rows(self, cohort_id: str) -> list[list[str]]:
        policy = self.policy
        rows = [
            [cohort_id, "genotype_quality_policy_hash", policy.policy_hash()],
            [cohort_id, "quality_evaluated_calls", str(self.evaluated)],
            [cohort_id, "quality_passed_calls", str(self.passed)],
            [cohort_id, "quality_kept_unevaluated_calls", str(self.kept_unevaluated)],
            [cohort_id, "quality_masked_calls", str(self.masked)],
            [cohort_id, "quality_masked_hom_ref_calls", str(self.masked_by_dosage["-1"])],
            [cohort_id, "quality_masked_het_calls", str(self.masked_by_dosage["0"])],
            [cohort_id, "quality_masked_hom_alt_calls", str(self.masked_by_dosage["1"])],
        ]
        by_metric = {
            quality.reason_metric(policy, reasons): count
            for reasons, count in self.reason_counts.items()
        }
        for metric in quality.reason_metrics(policy):
            rows.append([cohort_id, metric, str(by_metric.get(metric, 0))])
        for (key, _), counts in zip(policy.thresholds, self.status_counts):
            for status in quality.FIELD_STATUSES:
                rows.append([cohort_id, quality.status_metric(key, status), str(counts[status])])
        return rows


def _accounting_rows_from_counts(
    cohort_id: str,
    total_genotype_cells: int,
    counts: dict[str, int],
    phased_genotype_count: int,
    quality_counts: _QualityCounts | None = None,
) -> list[list[str]]:
    """Lay out the accounting TSV rows in their contractual order.

    With a quality policy, the historical rows keep their names and order:
    the three `standard_*_calls` rows count the calls that are still encoded
    as a dosage, and `total_treated_as_missing` gains the masked calls, so
    `total_genotype_cells = standard + total_treated_as_missing` holds in
    both modes. The quality rows are appended after them.
    """
    total_treated_as_missing = (
        counts["missing_calls"]
        + counts["non_diploid_calls_treated_as_missing"]
        + counts["non_biallelic_index_calls_treated_as_missing"]
    )
    if quality_counts is not None:
        total_treated_as_missing += quality_counts.masked

    rows = [[cohort_id, "total_genotype_cells", str(total_genotype_cells)]]
    for metric in _ACCOUNTING_METRIC_ORDER:
        rows.append([cohort_id, metric, str(counts[metric])])
    rows.append([cohort_id, "total_treated_as_missing", str(total_treated_as_missing)])
    rows.append([cohort_id, "phased_genotype_count", str(phased_genotype_count)])
    if quality_counts is not None:
        rows.extend(quality_counts.rows(cohort_id))

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

    if "genotype_quality_policy_hash" in accounting:
        lines.extend(_quality_summary_lines(accounting))

    return "\n".join(lines) + "\n"


def _quality_summary_lines(accounting: dict[str, str]) -> list[str]:
    """Issue #64: the genotype quality section, present only when a policy ran."""
    lines = [
        "",
        "Genotype quality policy (DP/GQ masking, separate from the site hard filter)",
        f"Policy hash: {accounting['genotype_quality_policy_hash']}",
        (
            "Evaluated calls (diploid biallelic-index calls with no missing allele): "
            f"{accounting['quality_evaluated_calls']}"
        ),
        f"  passed: {accounting['quality_passed_calls']}",
        f"  kept, a field could not be evaluated: {accounting['quality_kept_unevaluated_calls']}",
        f"  masked to nan: {accounting['quality_masked_calls']}",
        (
            f"    by original dosage: hom-ref {accounting['quality_masked_hom_ref_calls']}, "
            f"het {accounting['quality_masked_het_calls']}, "
            f"hom-alt {accounting['quality_masked_hom_alt_calls']}"
        ),
        "    by reason pair (each masked call in exactly one):",
    ]
    lines.extend(
        f"      {metric.removeprefix('quality_masked_reason.')}: {value}"
        for metric, value in accounting.items()
        if metric.startswith("quality_masked_reason.")
    )
    lines.append("  per-field status (each evaluated call in exactly one, per field):")
    lines.extend(
        f"    {metric}: {value}"
        for metric, value in accounting.items()
        if metric.startswith(("dp_status.", "gq_status."))
    )
    lines.append(
        "The standard hom-ref/het/hom-alt counts above are the calls still encoded "
        "after masking, and 'Total cells treated as missing' includes the masked "
        "calls. This is a sensitivity to the chosen thresholds, not a measure of "
        "genotype accuracy: no truth set was used."
    )
    return lines


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
                    format_field=fields[8],
                    sample_fields=tuple(fields[FIXED_COLUMN_COUNT:]),
                )
            )

    if sample_names is None:
        raise MalformedVcfError(f"{path}: no #CHROM header line found")

    return GsPassVcf(sample_names=sample_names, records=tuple(records))


@dataclass(frozen=True)
class _ResolvedCell:
    """Reference implementation only: one cell after any quality policy.

    ``category`` is the genotype-shape category, or ``quality_masked`` for a
    standard call the policy masked; ``original_dosage`` keeps what the call
    would have encoded to, so masked calls can be accounted by dosage.
    """

    category: str
    dosage: str
    is_phased: bool
    original_dosage: str
    verdict: quality.CallQuality | None


def _resolve_record_cells(
    record: GsPassRecord, policy: quality.GenotypeQualityPolicy | None
) -> list[_ResolvedCell]:
    evaluator = (
        quality.RowQualityEvaluator(policy, record.format_field) if policy is not None else None
    )
    cells: list[_ResolvedCell] = []
    for position, gt in enumerate(record.sample_genotypes):
        cell = classify_genotype(gt)
        if evaluator is None or cell.category != "standard":
            cells.append(
                _ResolvedCell(cell.category, cell.dosage, cell.is_phased, cell.dosage, None)
            )
            continue
        verdict = evaluator.evaluate(record.sample_fields[position].split(":"))
        if verdict.masked:
            cells.append(
                _ResolvedCell(
                    "quality_masked", MISSING_CELL_TOKEN, cell.is_phased, cell.dosage, verdict
                )
            )
        else:
            cells.append(
                _ResolvedCell("standard", cell.dosage, cell.is_phased, cell.dosage, verdict)
            )
    return cells


def build_matrix_rows(
    vcf: GsPassVcf, quality_policy: quality.GenotypeQualityPolicy | None = None
) -> list[list[str]]:
    """Build the genotype matrix's data rows: one row per variant."""
    rows: list[list[str]] = []

    for record in vcf.records:
        cells = _resolve_record_cells(record, quality_policy)
        rows.append([record.variant_key, *(cell.dosage for cell in cells)])

    return rows


def build_sample_metadata_rows(
    cohort_id: str, vcf: GsPassVcf, quality_policy: quality.GenotypeQualityPolicy | None = None
) -> list[list[str]]:
    """Build the sample metadata rows: one row per sample, in header order."""
    total_variants = len(vcf.records)
    missing_counts = [0] * len(vcf.sample_names)
    non_standard_counts = [0] * len(vcf.sample_names)
    masked_counts = [0] * len(vcf.sample_names)

    for record in vcf.records:
        for sample_position, cell in enumerate(_resolve_record_cells(record, quality_policy)):
            if cell.category == "standard":
                continue
            missing_counts[sample_position] += 1
            if cell.category == "quality_masked":
                masked_counts[sample_position] += 1
            elif cell.category != "missing":
                non_standard_counts[sample_position] += 1

    rows = []
    for sample_index, sample_id in enumerate(vcf.sample_names):
        row = _sample_metadata_row(
            cohort_id,
            sample_index,
            sample_id,
            missing_counts[sample_index],
            non_standard_counts[sample_index],
            total_variants,
        )
        if quality_policy is not None:
            row.append(str(masked_counts[sample_index]))
        rows.append(row)
    return rows


def build_variant_metadata_rows(
    cohort_id: str, vcf: GsPassVcf, quality_policy: quality.GenotypeQualityPolicy | None = None
) -> list[list[str]]:
    """Build the variant metadata rows: one row per variant, in file order."""
    total_samples = len(vcf.sample_names)
    rows: list[list[str]] = []

    for variant_index, record in enumerate(vcf.records):
        cells = _resolve_record_cells(record, quality_policy)
        missing_count = sum(1 for cell in cells if cell.category != "standard")
        row = _variant_metadata_row(
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
        if quality_policy is not None:
            row.append(str(sum(1 for cell in cells if cell.category == "quality_masked")))
        rows.append(row)

    return rows


def build_genotype_accounting_rows(
    cohort_id: str, vcf: GsPassVcf, quality_policy: quality.GenotypeQualityPolicy | None = None
) -> list[list[str]]:
    """Build the cohort-wide genotype-encoding accounting rows.

    ``phased_genotype_count`` is reported separately from every other
    metric here: it counts calls that were phased, regardless of
    whether they were standard, missing, or otherwise non-standard, and
    is never added into ``total_treated_as_missing`` -- a phased call
    that resolves to a real dosage is not missing.
    """
    counts = _new_accounting_counts()
    quality_counts = _QualityCounts(quality_policy) if quality_policy is not None else None

    total_genotype_cells = 0
    phased_genotype_count = 0
    for record in vcf.records:
        for cell in _resolve_record_cells(record, quality_policy):
            total_genotype_cells += 1
            if cell.is_phased:
                phased_genotype_count += 1
            if cell.verdict is not None and quality_counts is not None:
                quality_counts.record(cell.original_dosage, cell.verdict)
            if cell.category == "standard":
                counts[_METRIC_BY_DOSAGE[cell.dosage]] += 1
            elif cell.category != "quality_masked":
                counts[_METRIC_BY_NON_STANDARD_CATEGORY[cell.category]] += 1

    return _accounting_rows_from_counts(
        cohort_id, total_genotype_cells, counts, phased_genotype_count, quality_counts
    )


def build_genotype_accounting_summary_text(
    cohort_id: str, vcf: GsPassVcf, quality_policy: quality.GenotypeQualityPolicy | None = None
) -> str:
    """Build the human-readable genotype-encoding accounting summary."""
    accounting = {
        row[1]: row[2] for row in build_genotype_accounting_rows(cohort_id, vcf, quality_policy)
    }

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


#: Issue #64: BGZF, the blocked gzip htslib reads and indexes. The quality-
#: masked VCF is written here, by the same pass that writes the matrix, so a
#: second tool never re-serializes it: bcftools only builds its index.
_BGZF_MAX_BLOCK_INPUT_BYTES = 0xFF00
_BGZF_COMPRESS_LEVEL = 6
_BGZF_HEADER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


class _StreamingBgzfWriter:
    """Write BGZF incrementally: gzip members of at most 64 KiB, then the EOF block.

    Each block is one gzip member carrying the `BC` extra subfield with the
    block's own size, exactly as the SAM/BAM specification defines BGZF, so
    the result is readable by any gzip reader and indexable by `bcftools
    index`. Buffered input never exceeds one block, so memory does not grow
    with the number of rows. Deterministic: no timestamp or file name is
    written, so identical input reproduces identical bytes.
    """

    def __init__(self, path: Path) -> None:
        self._handle = path.open("wb")
        self._pending = bytearray()

    def write(self, text: str) -> None:
        self._pending += text.encode("utf-8")
        while len(self._pending) >= _BGZF_MAX_BLOCK_INPUT_BYTES:
            self._write_block(bytes(self._pending[:_BGZF_MAX_BLOCK_INPUT_BYTES]))
            del self._pending[:_BGZF_MAX_BLOCK_INPUT_BYTES]

    def _write_block(self, data: bytes) -> None:
        compressor = zlib.compressobj(_BGZF_COMPRESS_LEVEL, zlib.DEFLATED, -15)
        payload = compressor.compress(data) + compressor.flush()
        if len(payload) + 26 > 0x10000:
            # Incompressible input: a stored deflate block always fits.
            compressor = zlib.compressobj(0, zlib.DEFLATED, -15)
            payload = compressor.compress(data) + compressor.flush()
        self._handle.write(_BGZF_HEADER)
        self._handle.write(struct.pack("<H", len(payload) + 25))
        self._handle.write(payload)
        self._handle.write(struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data)))

    def close(self) -> None:
        try:
            if self._pending:
                self._write_block(bytes(self._pending))
                self._pending.clear()
            self._handle.write(_BGZF_EOF)
            self._handle.flush()
            os.fsync(self._handle.fileno())
        finally:
            self._handle.close()

    def __enter__(self) -> _StreamingBgzfWriter:
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
    #: Issue #64: present only when a genotype quality policy ran.
    quality_counts: _QualityCounts | None = None
    sample_masked_counts: list[int] | None = None


def _write_tsv_header(handle, header: tuple[str, ...]) -> None:
    handle.write("\t".join(header) + "\n")


def _write_tsv_row(handle, row: list[str]) -> None:
    handle.write("\t".join(row) + "\n")


def _split_row_calls(
    fields: list[str],
    sample_names: tuple[str, ...],
    path: Path,
    line_number: int,
) -> tuple[int, list[list[str]]]:
    """Issue #64: like `_extract_row_genotypes`, but keep every sample's subfields.

    The same shape checks and the same messages, so enabling a quality
    policy never changes which inputs are accepted as well-formed.
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
    calls: list[list[str]] = []
    for sample_position, sample_field in enumerate(fields[FIXED_COLUMN_COUNT:]):
        subfields = sample_field.split(":")
        if gt_index >= len(subfields):
            raise MalformedVcfError(
                f"{path}: line {line_number}: sample "
                f"'{sample_names[sample_position]}' has {len(subfields)} "
                f"FORMAT subfield(s), but FORMAT places GT at index {gt_index}"
            )
        calls.append(subfields)
    return gt_index, calls


_REQUIRED_INFO_DEFINITIONS: tuple[str, ...] = tuple(
    f"##INFO=<ID={key}," for key in quality.RECOMPUTED_INFO_KEYS
)


def stream_gs_panel(
    *,
    gs_pass_vcf: Path,
    cohort_id: str,
    matrix_path: Path,
    sample_metadata_path: Path,
    variant_metadata_path: Path,
    genotype_accounting_path: Path,
    genotype_accounting_summary_path: Path,
    quality_policy: quality.GenotypeQualityPolicy | None = None,
    quality_masked_vcf_path: Path | None = None,
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

    Issue #64: with `quality_policy`, the same pass also evaluates each
    standard call against the policy and writes the quality-masked VCF to
    `quality_masked_vcf_path`. Without it, this is the historical loop,
    unchanged, so a disabled policy cannot alter a single output byte.
    """
    if (quality_policy is None) != (quality_masked_vcf_path is None):
        raise ValueError("a quality policy and a quality-masked VCF path go together")
    if quality_policy is not None:
        assert quality_masked_vcf_path is not None
        return _stream_gs_panel_with_quality(
            gs_pass_vcf=gs_pass_vcf,
            cohort_id=cohort_id,
            matrix_path=matrix_path,
            sample_metadata_path=sample_metadata_path,
            variant_metadata_path=variant_metadata_path,
            genotype_accounting_path=genotype_accounting_path,
            genotype_accounting_summary_path=genotype_accounting_summary_path,
            policy=quality_policy,
            quality_masked_vcf_path=quality_masked_vcf_path,
        )

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
        variant_metadata = stack.enter_context(variant_metadata_path.open("w", encoding="utf-8"))
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
                raise MalformedVcfError(f"{gs_pass_vcf}: data row seen before #CHROM header")

            fields = line.split("\t")
            genotypes = _extract_row_genotypes(fields, sample_names, gs_pass_vcf, line_number)
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

            matrix.write("\t".join([_variant_key(chrom, pos, ref, alt), *dosages]) + "\n")
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


def _stream_gs_panel_with_quality(
    *,
    gs_pass_vcf: Path,
    cohort_id: str,
    matrix_path: Path,
    sample_metadata_path: Path,
    variant_metadata_path: Path,
    genotype_accounting_path: Path,
    genotype_accounting_summary_path: Path,
    policy: quality.GenotypeQualityPolicy,
    quality_masked_vcf_path: Path,
) -> _PanelTotals:
    """Issue #64: the same single pass, with a genotype quality policy.

    Kept as its own loop rather than as branches inside the historical one,
    so the disabled path stays exactly the code whose outputs are already
    published, and so nothing here costs the disabled path a per-cell test.

    For each row: every call is classified once; a standard call is judged
    once by the shared evaluator; a masked call becomes `nan` in the matrix
    and a missing GT in the VCF, from that one verdict. A row with no
    masked call is copied to the VCF verbatim. A row with any masked call
    gets its GTs replaced and AC/AN/AF recomputed from the masked GTs, and
    every other column and FORMAT subfield is kept as it was. Memory is
    still the single row plus per-sample and fixed-size counters.
    """
    sample_names: tuple[str, ...] | None = None
    variant_count = 0
    total_genotype_cells = 0
    phased_genotype_count = 0
    counts = _new_accounting_counts()
    quality_counts = _QualityCounts(policy)
    sample_missing_counts: list[int] = []
    sample_non_standard_counts: list[int] = []
    sample_masked_counts: list[int] = []
    declared_info: set[str] = set()
    header_line = quality.mask_header_line(policy)

    with ExitStack() as stack:
        handle = stack.enter_context(gzip.open(gs_pass_vcf, "rt", encoding="utf-8"))
        matrix = stack.enter_context(_StreamingGzipWriter(matrix_path))
        masked_vcf = stack.enter_context(_StreamingBgzfWriter(quality_masked_vcf_path))
        variant_metadata = stack.enter_context(variant_metadata_path.open("w", encoding="utf-8"))
        _write_tsv_header(variant_metadata, VARIANT_METADATA_HEADER_WITH_QUALITY)

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
                missing = [
                    definition
                    for definition in _REQUIRED_INFO_DEFINITIONS
                    if definition not in declared_info
                ]
                if missing:
                    raise MalformedVcfError(
                        f"{gs_pass_vcf}: the header does not define "
                        + ", ".join(d.removeprefix("##INFO=<ID=").rstrip(",") for d in missing)
                        + " as INFO; a genotype quality policy recomputes AC, AN and AF on "
                        "rows with masked calls, so they must be declared"
                    )
                sample_missing_counts = [0] * len(sample_names)
                sample_non_standard_counts = [0] * len(sample_names)
                sample_masked_counts = [0] * len(sample_names)
                matrix.write("\t".join(["variant_key", *sample_names]) + "\n")
                masked_vcf.write(header_line + "\n")
                masked_vcf.write(line + "\n")
                continue

            if line.startswith("#"):
                if sample_names is not None:
                    # The historical reader ignores a stray header line after
                    # #CHROM; enabling a policy must not change which inputs
                    # are accepted, and it must not write one into the body.
                    continue
                for definition in _REQUIRED_INFO_DEFINITIONS:
                    if line.startswith(definition):
                        declared_info.add(definition)
                masked_vcf.write(line + "\n")
                continue

            if sample_names is None:
                raise MalformedVcfError(f"{gs_pass_vcf}: data row seen before #CHROM header")

            fields = line.split("\t")
            gt_index, calls = _split_row_calls(fields, sample_names, gs_pass_vcf, line_number)
            chrom, pos, _id, ref, alt, qual = fields[:6]
            evaluator = quality.RowQualityEvaluator(policy, fields[8])

            dosages: list[str] = []
            row_missing_count = 0
            row_masked_count = 0
            for sample_position, subfields in enumerate(calls):
                gt = subfields[gt_index]
                cell = classify_genotype(gt)

                total_genotype_cells += 1
                if cell.is_phased:
                    phased_genotype_count += 1

                if cell.category != "standard":
                    dosages.append(cell.dosage)
                    counts[_METRIC_BY_NON_STANDARD_CATEGORY[cell.category]] += 1
                    row_missing_count += 1
                    sample_missing_counts[sample_position] += 1
                    if cell.category != "missing":
                        sample_non_standard_counts[sample_position] += 1
                    continue

                try:
                    verdict = evaluator.evaluate(subfields)
                except quality.GenotypeQualityPolicyError as error:
                    raise quality.GenotypeQualityPolicyError(
                        f"{gs_pass_vcf}: line {line_number}: {chrom}:{pos} sample "
                        f"'{sample_names[sample_position]}': {error}"
                    ) from error
                quality_counts.record(cell.dosage, verdict)

                if not verdict.masked:
                    dosages.append(cell.dosage)
                    counts[_METRIC_BY_DOSAGE[cell.dosage]] += 1
                    continue

                dosages.append(MISSING_CELL_TOKEN)
                subfields[gt_index] = quality.masked_genotype(gt)
                row_masked_count += 1
                row_missing_count += 1
                sample_missing_counts[sample_position] += 1
                sample_masked_counts[sample_position] += 1

            if row_masked_count:
                ac, an = quality.allele_counts([subfields[gt_index] for subfields in calls])
                try:
                    info = quality.rewrite_allele_info(fields[7], ac, an)
                except ValueError as error:
                    raise MalformedVcfError(
                        f"{gs_pass_vcf}: line {line_number}: {chrom}:{pos}: {error}"
                    ) from error
                masked_vcf.write(
                    "\t".join(
                        [
                            *fields[:7],
                            info,
                            fields[8],
                            *(":".join(subfields) for subfields in calls),
                        ]
                    )
                    + "\n"
                )
            else:
                masked_vcf.write(line + "\n")

            matrix.write("\t".join([_variant_key(chrom, pos, ref, alt), *dosages]) + "\n")
            row = _variant_metadata_row(
                cohort_id,
                variant_count,
                chrom,
                pos,
                ref,
                alt,
                qual,
                row_missing_count,
                len(sample_names),
            )
            row.append(str(row_masked_count))
            _write_tsv_row(variant_metadata, row)
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
        quality_counts=quality_counts,
        sample_masked_counts=sample_masked_counts,
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
    sample_rows = []
    for sample_index, sample_id in enumerate(totals.sample_names):
        row = _sample_metadata_row(
            cohort_id,
            sample_index,
            sample_id,
            totals.sample_missing_counts[sample_index],
            totals.sample_non_standard_counts[sample_index],
            totals.variant_count,
        )
        if totals.sample_masked_counts is not None:
            row.append(str(totals.sample_masked_counts[sample_index]))
        sample_rows.append(row)
    header = (
        SAMPLE_METADATA_HEADER
        if totals.quality_counts is None
        else SAMPLE_METADATA_HEADER_WITH_QUALITY
    )
    write_tsv(sample_metadata_path, header, sample_rows)

    accounting_rows = _accounting_rows_from_counts(
        cohort_id,
        totals.total_genotype_cells,
        totals.counts,
        totals.phased_genotype_count,
        totals.quality_counts,
    )
    write_tsv(genotype_accounting_path, GENOTYPE_ACCOUNTING_HEADER, accounting_rows)

    accounting = {row[1]: row[2] for row in accounting_rows}
    genotype_accounting_summary_path.write_text(
        _summary_text_from_accounting(
            cohort_id, accounting, totals.variant_count, len(totals.sample_names)
        ),
        encoding="utf-8",
    )


def write_policy_document(path: Path, policy: quality.GenotypeQualityPolicy) -> None:
    """Issue #64: the machine-readable policy, byte-deterministic for one policy."""
    path.write_text(
        json.dumps(policy.document(), indent=2, sort_keys=True, ensure_ascii=True) + "\n",
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

    # Issue #64: an optional genotype quality policy. Off unless
    # --genotype-quality-mask is given; no threshold has a default, and every
    # other option below is refused when the mask is off rather than ignored.
    group = parser.add_argument_group(
        "genotype quality policy (Issue #64; off by default, no threshold is recommended)"
    )
    group.add_argument(
        "--genotype-quality-mask",
        action="store_true",
        help="Enable genotype-level DP/GQ masking. Requires --genotype-min-dp and/or "
        "--genotype-min-gq, --quality-masked-vcf-output and --genotype-quality-policy-output.",
    )
    group.add_argument(
        "--genotype-min-dp",
        type=int,
        default=None,
        help="Mask a call whose FORMAT/DP is below this value; the value itself passes.",
    )
    group.add_argument(
        "--genotype-min-gq",
        type=int,
        default=None,
        help="Mask a call whose FORMAT/GQ is below this value; the value itself passes.",
    )
    for option, what in (
        ("missing-format-field", "a configured field is absent from the row's FORMAT"),
        ("missing-value", "the sample's value is '.' or its subfields stop before it"),
        ("malformed-value", "the value is not a non-negative decimal integer"),
    ):
        group.add_argument(
            f"--genotype-{option}-policy",
            choices=quality.POLICY_ACTIONS,
            default=None,
            help=f"What to do when {what}: reject (fail the build), unevaluated (keep "
            "the call and count it), or mask. Defaults to reject when the mask is enabled.",
        )
    group.add_argument(
        "--quality-masked-vcf-output",
        type=Path,
        default=None,
        help="Output path for the BGZF quality-masked VCF.",
    )
    group.add_argument(
        "--genotype-quality-policy-output",
        type=Path,
        default=None,
        help="Output path for the machine-readable policy JSON.",
    )

    return parser.parse_args(argv)


def _quality_policy_from_args(args: argparse.Namespace) -> quality.GenotypeQualityPolicy | None:
    """Resolve the CLI's quality options, refusing any silently ignored one."""
    policy_options = {
        "--genotype-min-dp": args.genotype_min_dp,
        "--genotype-min-gq": args.genotype_min_gq,
        "--genotype-missing-format-field-policy": args.genotype_missing_format_field_policy,
        "--genotype-missing-value-policy": args.genotype_missing_value_policy,
        "--genotype-malformed-value-policy": args.genotype_malformed_value_policy,
        "--quality-masked-vcf-output": args.quality_masked_vcf_output,
        "--genotype-quality-policy-output": args.genotype_quality_policy_output,
    }
    if not args.genotype_quality_mask:
        given = [name for name, value in policy_options.items() if value is not None]
        if given:
            raise quality.InvalidPolicyError(
                f"{', '.join(given)} given without --genotype-quality-mask; the mask is off, "
                "so these would be silently ignored"
            )
        return None

    if args.quality_masked_vcf_output is None or args.genotype_quality_policy_output is None:
        raise quality.InvalidPolicyError(
            "--genotype-quality-mask needs --quality-masked-vcf-output and "
            "--genotype-quality-policy-output"
        )
    return quality.GenotypeQualityPolicy(
        min_dp=args.genotype_min_dp,
        min_gq=args.genotype_min_gq,
        missing_format_field=args.genotype_missing_format_field_policy or quality.REJECT,
        missing_value=args.genotype_missing_value_policy or quality.REJECT,
        malformed_value=args.genotype_malformed_value_policy or quality.REJECT,
    )


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

    try:
        quality_policy = _quality_policy_from_args(args)
    except quality.InvalidPolicyError as error:
        print(f"build_gs_panel.py: error: {error}", file=sys.stderr)
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
    if quality_policy is not None:
        # Published with the panel, or not at all: a masked VCF next to an
        # unmasked matrix from an earlier run would be two lineages at once.
        final_paths.extend([args.quality_masked_vcf_output, args.genotype_quality_policy_output])

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
                quality_policy=quality_policy,
                quality_masked_vcf_path=(
                    _staging_path(args.quality_masked_vcf_output)
                    if quality_policy is not None
                    else None
                ),
            )
            if quality_policy is not None:
                write_policy_document(
                    _staging_path(args.genotype_quality_policy_output), quality_policy
                )
        except OSError as error:
            print(
                f"build_gs_panel.py: error: cannot read {args.gs_pass_vcf}: {error}",
                file=sys.stderr,
            )
            return 1
        except (MalformedVcfError, quality.GenotypeQualityPolicyError) as error:
            print(f"build_gs_panel.py: error: {error}", file=sys.stderr)
            return 1

        publish_outputs(final_paths)
    finally:
        discard_staged_outputs(final_paths)

    return 0


if __name__ == "__main__":
    sys.exit(main())
