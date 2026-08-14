#!/usr/bin/env python3
"""Reconcile record counts across the full GS panel lineage.

No single Nextflow task in the GS lineage (raw/all -> normalize ->
classify -> hard-filter -> matrix) ever sees every stage's record count
at once, so this tool re-counts raw/all, the normalized VCF, and the
GS-eligible PASS VCF directly (independent of any other script's
numbers, for the same reason `bin/reconcile_variant_type_counts.py`
re-counts rather than trusting another script's totals), cross-
references the classification stage's own accounting (re-deriving its
shape-classification and duplicate-key logic here would duplicate that
logic rather than verify it), and -- critically -- reads the actual
genotype matrix and metadata files rather than assuming they agree with
the VCF they were built from. An earlier revision of this tool derived
`final_matrix_variant_records` from the variant metadata file's row
count without ever opening the matrix itself, which meant a bug that
corrupted only the matrix (a dropped row, a misaligned column, a
truncated write) would go completely undetected as long as the metadata
files still looked correct.

A later revision closed that gap for *counts*, but a subtler bug
remained possible: the matrix, the PASS VCF, and the metadata files
could all have the *same number* of rows while disagreeing on *which*
variant or sample each row actually describes -- a shuffled sample
column, a variant_key swapped for a different one, or the matrix and
variant metadata simply listing their rows in a different order.
Counting alone cannot catch any of this, so this tool now compares the
actual, ordered sequence of sample IDs and variant keys across every
artifact, not just their lengths, and also checks that every matrix
data row has exactly `1 + sample_count` columns (a column-count check
`summarize_matrix()` did not previously perform at all).

The full reconciliation is:

    raw_all_records
      -> normalized_records (post bcftools norm -m- splitting)
      -> classified_biallelic_snp_records (post shape reclassification
         and duplicate-key exclusion; from classify_normalized_variants.py)
      -> gs_pass_records (post GS-specific hard filtering)
      -> matrix_variant_records / variant_metadata_records (both read
         directly, and must agree with gs_pass_records and each other
         in count, identity, AND order)

Sample counts are reconciled the same way: the GS-eligible PASS VCF's
own sample IDs, in header order, must exactly match the matrix header
and the sample metadata file's own `sample_id` column, in order.

Every one of these cross-file checks is a **hard error**
(`InconsistentGsPanelError`, exit 1, no output written), not a warning:
there is no legitimate scientific scenario where the matrix and its own
source VCF disagree on record or sample counts, identity, or order --
unlike, say, a negative `records_not_selected` in the primary lineage (a
real consequence of GATK's multiallelic handling), any mismatch here can
only mean a bug in this pipeline's own wiring or scripts, and should
stop the run rather than be reported alongside a set of numbers that
cannot all be trusted at once.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

OUTPUT_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")
CLASSIFIED_RECORDS_METRIC = "output_records"
SAMPLE_ID_COLUMN = "sample_id"
VARIANT_KEY_COLUMN = "variant_key"


class MalformedVcfError(Exception):
    """Raised when a VCF cannot be summarized safely."""


class MalformedMatrixError(Exception):
    """Raised when the genotype matrix cannot be summarized safely."""


class MalformedAccountingError(Exception):
    """Raised when an upstream accounting or metadata file is missing data."""


class InconsistentGsPanelError(Exception):
    """Raised when the GS panel's own artifacts disagree with each other.

    Always a bug, never a legitimate scientific outcome -- see the
    module docstring for why this is fatal rather than a warning.
    """


@dataclass(frozen=True)
class VcfSummary:
    """Record/sample counts, identity, and row-shape read directly from one VCF."""

    record_count: int
    sample_count: int
    sample_ids: tuple[str, ...]
    variant_keys: tuple[str, ...]
    all_rows_have_expected_sample_count: bool


@dataclass(frozen=True)
class MatrixSummary:
    """Variant/sample counts, identity, and row-shape read directly from the matrix."""

    variant_count: int
    sample_count: int
    sample_ids: tuple[str, ...]
    variant_keys: tuple[str, ...]
    all_rows_have_expected_width: bool


@dataclass(frozen=True)
class ReconciliationResult:
    """Every metric this tool computes, once all consistency checks pass."""

    raw_all_records: int
    normalized_records: int
    classified_biallelic_snp_records: int
    gs_hard_filter_excluded_records: int
    gs_pass_records: int
    gs_pass_sample_count: int
    matrix_variant_records: int
    matrix_sample_count: int
    variant_metadata_records: int
    sample_metadata_records: int
    panel_status: str


def summarize_vcf(path: Path) -> VcfSummary:
    """Read sample IDs and per-row variant keys from a bgzipped VCF.

    ``all_rows_have_expected_sample_count`` defends against a
    truncated or otherwise malformed data row that has the right
    number of fixed columns but the wrong number of sample columns --
    a shape error `MalformedVcfError`'s plain field-count check alone
    would not catch, since it only requires *at least* 10 fields. A
    row's variant key is still recorded even when its sample-field
    count is wrong: the CHROM/POS/REF/ALT columns are unaffected by a
    malformed sample-field count, and withholding the key would only
    make the resulting hard error harder to diagnose.
    """
    sample_ids: tuple[str, ...] | None = None
    sample_count: int | None = None
    variant_keys: list[str] = []
    all_rows_consistent = True

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) <= 9:
                    raise MalformedVcfError(
                        f"{path}: #CHROM header has {len(fields)} fields, expected "
                        "at least 10 (9 fixed columns plus one or more samples)"
                    )
                sample_ids = tuple(fields[9:])
                sample_count = len(sample_ids)
                continue

            if line.startswith("#"):
                continue

            if sample_count is None:
                raise MalformedVcfError(f"{path}: data row seen before #CHROM header")

            fields = line.split("\t")
            if len(fields) < 10:
                raise MalformedVcfError(
                    f"{path}: data row has {len(fields)} tab-separated fields, "
                    "expected at least 10"
                )

            variant_keys.append(f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}")
            if len(fields) - 9 != sample_count:
                all_rows_consistent = False

    if sample_ids is None or sample_count is None:
        raise MalformedVcfError(f"{path}: no #CHROM header line found")

    return VcfSummary(
        record_count=len(variant_keys),
        sample_count=sample_count,
        sample_ids=sample_ids,
        variant_keys=tuple(variant_keys),
        all_rows_have_expected_sample_count=all_rows_consistent,
    )


def summarize_matrix(path: Path) -> MatrixSummary:
    """Read the genotype matrix's own header, sample IDs, and variant keys.

    ``all_rows_have_expected_width`` defends against a data row with
    too few or too many dosage columns -- for example one dropped
    column from a truncated write -- which the previous revision of
    this function could not detect at all, since it only counted
    non-empty lines without ever checking a single row's field count
    against the header's own declared sample count.
    """
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n")

        if not header_line:
            raise MalformedMatrixError(f"{path}: file is empty")

        header_fields = header_line.split("\t")
        if header_fields[0] != "variant_key":
            raise MalformedMatrixError(
                f"{path}: expected 'variant_key' as the first header column, "
                f"got {header_fields[0]!r}"
            )

        sample_ids = tuple(header_fields[1:])
        sample_count = len(sample_ids)

        variant_keys: list[str] = []
        all_rows_expected_width = True

        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")
            variant_keys.append(fields[0])
            if len(fields) != 1 + sample_count:
                all_rows_expected_width = False

    return MatrixSummary(
        variant_count=len(variant_keys),
        sample_count=sample_count,
        sample_ids=sample_ids,
        variant_keys=tuple(variant_keys),
        all_rows_have_expected_width=all_rows_expected_width,
    )


def read_accounting_metric(path: Path, metric: str) -> str:
    """Read one metric's value from a (cohort_id, metric, value) TSV."""
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            raise MalformedAccountingError(
                f"{path}: row has {len(fields)} tab-separated fields, expected at least 3"
            )

        if fields[1] == metric:
            return fields[2]

    raise MalformedAccountingError(f"{path}: missing required metric: {metric}")


def read_metadata_column(path: Path, column: str) -> tuple[str, ...]:
    """Read one named column's values, in row order, from a metadata TSV.

    Locating the column by name (rather than a hardcoded index) means
    this keeps working even if `build_gs_panel.py` ever reorders its
    own metadata columns.
    """
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    header = lines[0].split("\t")
    if column not in header:
        raise MalformedAccountingError(f"{path}: missing required column: {column}")
    index = header.index(column)

    values: list[str] = []
    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) <= index:
            raise MalformedAccountingError(
                f"{path}: row has {len(fields)} tab-separated fields, expected "
                f"at least {index + 1} to read column {column!r}"
            )

        values.append(fields[index])

    return tuple(values)


def _describe_sequence_mismatch(labeled: list[tuple[str, tuple[str, ...]]]) -> str:
    """Build a bounded-size diagnostic for 2+ named sequences that disagree.

    Never dumps every element: a GS panel can hold thousands of
    variants, so the message instead reports either a length
    disagreement or the first index at which the sequences diverge --
    always enough to start debugging, never proportional to panel size.
    """
    lengths = {name: len(seq) for name, seq in labeled}
    if len(set(lengths.values())) > 1:
        length_desc = ", ".join(f"{name}={length}" for name, length in lengths.items())
        return f"lengths differ ({length_desc})"

    common_length = next(iter(lengths.values()))
    for index in range(common_length):
        values_at_index = {name: seq[index] for name, seq in labeled}
        if len(set(values_at_index.values())) > 1:
            detail = ", ".join(
                f"{name}[{index}]={value!r}" for name, value in values_at_index.items()
            )
            return f"same length ({common_length}), first disagreement at index {index}: {detail}"

    return "sequences are equal"


def reconcile(
    *,
    raw_all: VcfSummary,
    normalized: VcfSummary,
    classified_biallelic_snp_records: int,
    gs_pass: VcfSummary,
    matrix: MatrixSummary,
    variant_metadata_keys: tuple[str, ...],
    sample_metadata_ids: tuple[str, ...],
) -> ReconciliationResult:
    """Cross-check every GS panel artifact and compute the full lineage.

    Raises `InconsistentGsPanelError` immediately on the first
    disagreement found, rather than collecting every mismatch: once one
    artifact disagrees with its own source, none of the remaining
    numbers can be trusted either. Identity/order checks (tuple
    equality) subsume the plain count checks they replace -- two
    sequences cannot be equal without also being the same length -- so
    there is no separate, redundant count-only comparison left here.
    """
    if not raw_all.all_rows_have_expected_sample_count:
        raise InconsistentGsPanelError(
            "raw/all VCF has at least one data row whose sample-field count "
            "does not match its own #CHROM header"
        )
    if not normalized.all_rows_have_expected_sample_count:
        raise InconsistentGsPanelError(
            "normalized VCF has at least one data row whose sample-field count "
            "does not match its own #CHROM header"
        )
    if not gs_pass.all_rows_have_expected_sample_count:
        raise InconsistentGsPanelError(
            "GS-eligible PASS VCF has at least one data row whose sample-field "
            "count does not match its own #CHROM header"
        )
    if not matrix.all_rows_have_expected_width:
        raise InconsistentGsPanelError(
            "genotype matrix has at least one data row whose column count does "
            "not equal '1 + sample_count' implied by its own header"
        )

    gs_hard_filter_excluded_records = classified_biallelic_snp_records - gs_pass.record_count
    if gs_hard_filter_excluded_records < 0:
        raise InconsistentGsPanelError(
            "gs_hard_filter_excluded_records is negative "
            f"(classified_biallelic_snp_records={classified_biallelic_snp_records}, "
            f"gs_pass_records={gs_pass.record_count}): hard filtering can only "
            "remove records, never add them, so this can only be a bug"
        )

    if not (gs_pass.sample_ids == matrix.sample_ids == sample_metadata_ids):
        raise InconsistentGsPanelError(
            "sample identity and/or order disagree across the GS panel's own "
            "artifacts: "
            + _describe_sequence_mismatch(
                [
                    ("gs_pass_vcf", gs_pass.sample_ids),
                    ("matrix", matrix.sample_ids),
                    ("sample_metadata", sample_metadata_ids),
                ]
            )
        )

    if not (gs_pass.variant_keys == matrix.variant_keys == variant_metadata_keys):
        raise InconsistentGsPanelError(
            "variant identity and/or order disagree across the GS panel's own "
            "artifacts: "
            + _describe_sequence_mismatch(
                [
                    ("gs_pass_vcf", gs_pass.variant_keys),
                    ("matrix", matrix.variant_keys),
                    ("variant_metadata", variant_metadata_keys),
                ]
            )
        )

    return ReconciliationResult(
        raw_all_records=raw_all.record_count,
        normalized_records=normalized.record_count,
        classified_biallelic_snp_records=classified_biallelic_snp_records,
        gs_hard_filter_excluded_records=gs_hard_filter_excluded_records,
        gs_pass_records=gs_pass.record_count,
        gs_pass_sample_count=gs_pass.sample_count,
        matrix_variant_records=matrix.variant_count,
        matrix_sample_count=matrix.sample_count,
        variant_metadata_records=len(variant_metadata_keys),
        sample_metadata_records=len(sample_metadata_ids),
        panel_status="empty" if matrix.variant_count == 0 else "populated",
    )


def build_output_rows(cohort_id: str, result: ReconciliationResult) -> list[list[str]]:
    """Build the reconciliation TSV data rows.

    Every count here is independently sourced (VCF, matrix, or
    metadata file); by the time this function is called, `reconcile()`
    has already confirmed they all agree in count, identity, and
    order, so the values shown are not merely one script's opinion of
    the panel's shape.
    """
    return [
        [cohort_id, "raw_all_records", str(result.raw_all_records)],
        [cohort_id, "normalized_records", str(result.normalized_records)],
        [
            cohort_id,
            "classified_biallelic_snp_records",
            str(result.classified_biallelic_snp_records),
        ],
        [
            cohort_id,
            "gs_hard_filter_excluded_records",
            str(result.gs_hard_filter_excluded_records),
        ],
        [cohort_id, "gs_pass_records", str(result.gs_pass_records)],
        [cohort_id, "gs_pass_sample_count", str(result.gs_pass_sample_count)],
        [cohort_id, "matrix_variant_records", str(result.matrix_variant_records)],
        [cohort_id, "matrix_sample_count", str(result.matrix_sample_count)],
        [cohort_id, "variant_metadata_records", str(result.variant_metadata_records)],
        [cohort_id, "sample_metadata_records", str(result.sample_metadata_records)],
        [cohort_id, "panel_status", result.panel_status],
    ]


def build_summary_text(cohort_id: str, result: ReconciliationResult) -> str:
    """Build the human-readable GS panel record-accounting summary."""
    lines = [
        "GS panel record accounting summary",
        f"Cohort ID: {cohort_id}",
        f"raw/all records: {result.raw_all_records}",
        f"normalized records (post bcftools norm -m- splitting): {result.normalized_records}",
        (
            "classified biallelic-SNP records (post shape reclassification "
            f"and duplicate-key exclusion): {result.classified_biallelic_snp_records}"
        ),
        (
            "GS hard-filter excluded records = classified_biallelic_snp_records "
            f"({result.classified_biallelic_snp_records}) - gs_pass_records "
            f"({result.gs_pass_records}) = {result.gs_hard_filter_excluded_records}"
        ),
        f"GS-eligible PASS records: {result.gs_pass_records} ({result.gs_pass_sample_count} sample(s))",
        (
            "Cross-checked directly against the panel's own artifacts "
            f"(never assumed to agree): matrix has {result.matrix_variant_records} "
            f"variant row(s) and {result.matrix_sample_count} sample column(s); "
            f"variant metadata has {result.variant_metadata_records} row(s); "
            f"sample metadata has {result.sample_metadata_records} row(s). "
            "The counts, the identity (variant_key / sample ID), and the row "
            "order all agree across every artifact -- this file would not "
            "exist otherwise, since any disagreement on any of the three is "
            "a hard error."
        ),
        f"panel_status: {result.panel_status}",
    ]

    if result.panel_status == "empty":
        lines.append(
            "This is a normal outcome, not an error: zero GS-eligible PASS "
            "records means every record was excluded by normalization, "
            "classification, or hard filtering, not that the run failed. "
            f"The sample list ({result.sample_metadata_records} sample(s)) "
            "is preserved in the matrix and metadata regardless."
        )

    return "\n".join(lines) + "\n"


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the reconciliation CLI."""
    parser = argparse.ArgumentParser(
        description="Reconcile record counts across the full GS panel lineage."
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--raw-all-vcf", required=True, type=Path, help="Path to cohort.raw.vcf.gz."
    )
    parser.add_argument(
        "--normalized-vcf",
        required=True,
        type=Path,
        help="Path to the GS_NORMALIZE_VARIANTS output VCF.",
    )
    parser.add_argument(
        "--normalization-accounting",
        required=True,
        type=Path,
        help="Path to classify_normalized_variants.py's accounting TSV.",
    )
    parser.add_argument(
        "--gs-pass-vcf",
        required=True,
        type=Path,
        help="Path to the GS-eligible PASS VCF (cohort_gs.snp.pass.vcf.gz).",
    )
    parser.add_argument(
        "--matrix",
        required=True,
        type=Path,
        help="Path to build_gs_panel.py's gzipped genotype matrix.",
    )
    parser.add_argument(
        "--variant-metadata",
        required=True,
        type=Path,
        help="Path to build_gs_panel.py's variant metadata TSV.",
    )
    parser.add_argument(
        "--sample-metadata",
        required=True,
        type=Path,
        help="Path to build_gs_panel.py's sample metadata TSV.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Output path for the reconciliation TSV."
    )
    parser.add_argument(
        "--summary-output",
        required=True,
        type=Path,
        help="Output path for the human-readable summary text file.",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        raw_all = summarize_vcf(args.raw_all_vcf)
        normalized = summarize_vcf(args.normalized_vcf)
        classified_biallelic_snp_records = int(
            read_accounting_metric(args.normalization_accounting, CLASSIFIED_RECORDS_METRIC)
        )
        gs_pass = summarize_vcf(args.gs_pass_vcf)
        matrix = summarize_matrix(args.matrix)
        variant_metadata_keys = read_metadata_column(args.variant_metadata, VARIANT_KEY_COLUMN)
        sample_metadata_ids = read_metadata_column(args.sample_metadata, SAMPLE_ID_COLUMN)
    except OSError as error:
        print(f"reconcile_gs_panel_accounting.py: error: {error}", file=sys.stderr)
        return 1
    except (MalformedVcfError, MalformedMatrixError, MalformedAccountingError) as error:
        print(f"reconcile_gs_panel_accounting.py: error: {error}", file=sys.stderr)
        return 1

    try:
        result = reconcile(
            raw_all=raw_all,
            normalized=normalized,
            classified_biallelic_snp_records=classified_biallelic_snp_records,
            gs_pass=gs_pass,
            matrix=matrix,
            variant_metadata_keys=variant_metadata_keys,
            sample_metadata_ids=sample_metadata_ids,
        )
    except InconsistentGsPanelError as error:
        print(f"reconcile_gs_panel_accounting.py: error: {error}", file=sys.stderr)
        return 1

    write_tsv(args.output, OUTPUT_HEADER, build_output_rows(args.cohort_id, result))
    args.summary_output.write_text(
        build_summary_text(args.cohort_id, result),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
