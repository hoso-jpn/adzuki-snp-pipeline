#!/usr/bin/env python3
"""Reconcile record counts across the full GS panel lineage.

No single Nextflow task in the GS lineage (raw/all -> normalize ->
classify -> hard-filter -> matrix) ever sees every stage's record count
at once, so this tool re-counts raw/all and the normalized VCF directly
(independent of any other script's numbers, for the same reason
`bin/reconcile_variant_type_counts.py` re-counts rather than trusting
another script's totals), cross-references the classification stage's
own accounting (re-deriving its shape-classification and duplicate-key
logic here would duplicate that logic rather than verify it), and
counts the final GS-pass VCF and the matrix's own metadata files, to
report one end-to-end reconciliation:

    raw_all_records
      -> normalized_records (post bcftools norm -m- splitting)
      -> classified_biallelic_snp_records (post shape reclassification
         and duplicate-key exclusion; from classify_normalized_variants.py)
      -> gs_pass_records (post GS-specific hard filtering)
      -> final_matrix_variant_records (should equal gs_pass_records)

`gs_hard_filter_excluded_records` is computed as
`classified_biallelic_snp_records - gs_pass_records` and is never
expected to be negative (hard filtering can only remove records, not
add them); this tool still reports it as computed and flags it rather
than hiding a surprise, matching this repository's existing convention
for record-accounting tools.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

OUTPUT_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")
CLASSIFIED_RECORDS_METRIC = "output_records"


class MalformedVcfError(Exception):
    """Raised when a VCF cannot be counted safely."""


class MalformedAccountingError(Exception):
    """Raised when an upstream accounting or metadata file is missing data."""


@dataclass(frozen=True)
class ReconciliationResult:
    """Every metric this tool computes."""

    raw_all_records: int
    normalized_records: int
    classified_biallelic_snp_records: int
    gs_hard_filter_excluded_records: int
    gs_hard_filter_excluded_is_negative: bool
    gs_pass_records: int
    final_matrix_variant_records: int
    final_matrix_variant_count_matches_pass: bool
    final_matrix_sample_count: int
    panel_status: str


def count_vcf_records(path: Path) -> int:
    """Count data rows (non-header, non-blank lines) in a bgzipped VCF."""
    count = 0

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            count += 1

    return count


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


def count_metadata_rows(path: Path) -> int:
    """Count data rows (excluding the header) in a metadata TSV."""
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    return sum(1 for line in lines[1:] if line)


def reconcile(
    *,
    raw_all_records: int,
    normalized_records: int,
    classified_biallelic_snp_records: int,
    gs_pass_records: int,
    final_matrix_variant_records: int,
    final_matrix_sample_count: int,
) -> ReconciliationResult:
    """Compute the full GS panel lineage reconciliation."""
    gs_hard_filter_excluded_records = (
        classified_biallelic_snp_records - gs_pass_records
    )

    return ReconciliationResult(
        raw_all_records=raw_all_records,
        normalized_records=normalized_records,
        classified_biallelic_snp_records=classified_biallelic_snp_records,
        gs_hard_filter_excluded_records=gs_hard_filter_excluded_records,
        gs_hard_filter_excluded_is_negative=gs_hard_filter_excluded_records < 0,
        gs_pass_records=gs_pass_records,
        final_matrix_variant_records=final_matrix_variant_records,
        final_matrix_variant_count_matches_pass=(
            final_matrix_variant_records == gs_pass_records
        ),
        final_matrix_sample_count=final_matrix_sample_count,
        panel_status="empty" if final_matrix_variant_records == 0 else "populated",
    )


def build_output_rows(cohort_id: str, result: ReconciliationResult) -> list[list[str]]:
    """Build the reconciliation TSV data rows."""
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
        [
            cohort_id,
            "gs_hard_filter_excluded_is_negative",
            "true" if result.gs_hard_filter_excluded_is_negative else "false",
        ],
        [cohort_id, "gs_pass_records", str(result.gs_pass_records)],
        [
            cohort_id,
            "final_matrix_variant_records",
            str(result.final_matrix_variant_records),
        ],
        [
            cohort_id,
            "final_matrix_variant_count_matches_pass",
            "true" if result.final_matrix_variant_count_matches_pass else "false",
        ],
        [cohort_id, "final_matrix_sample_count", str(result.final_matrix_sample_count)],
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
    ]

    if result.gs_hard_filter_excluded_is_negative:
        lines.append(
            "WARNING: gs_hard_filter_excluded_records is negative. Hard "
            "filtering can only remove records, never add them, so this "
            "means gs_pass_records exceeds classified_biallelic_snp_records "
            "-- investigate a possible mismatch between the two inputs "
            "(for example, files from different runs)."
        )

    lines.append(f"GS-eligible PASS records: {result.gs_pass_records}")
    lines.append(
        f"Final matrix variant records: {result.final_matrix_variant_records} "
        f"(samples: {result.final_matrix_sample_count})"
    )

    if not result.final_matrix_variant_count_matches_pass:
        lines.append(
            "WARNING: final_matrix_variant_records does not match "
            "gs_pass_records. The matrix should describe exactly the "
            "GS-eligible PASS records; investigate a possible mismatch "
            "between the matrix-building input and the PASS VCF."
        )

    lines.append(f"panel_status: {result.panel_status}")

    if result.panel_status == "empty":
        lines.append(
            "This is a normal outcome, not an error: zero GS-eligible PASS "
            "records means every record was excluded by normalization, "
            "classification, or hard filtering, not that the run failed. "
            f"The sample list ({result.final_matrix_sample_count} sample(s)) "
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
        raw_all_records = count_vcf_records(args.raw_all_vcf)
        normalized_records = count_vcf_records(args.normalized_vcf)
        classified_biallelic_snp_records = int(
            read_accounting_metric(args.normalization_accounting, CLASSIFIED_RECORDS_METRIC)
        )
        gs_pass_records = count_vcf_records(args.gs_pass_vcf)
        final_matrix_variant_records = count_metadata_rows(args.variant_metadata)
        final_matrix_sample_count = count_metadata_rows(args.sample_metadata)
    except OSError as error:
        print(f"reconcile_gs_panel_accounting.py: error: {error}", file=sys.stderr)
        return 1
    except (MalformedVcfError, MalformedAccountingError) as error:
        print(f"reconcile_gs_panel_accounting.py: error: {error}", file=sys.stderr)
        return 1

    result = reconcile(
        raw_all_records=raw_all_records,
        normalized_records=normalized_records,
        classified_biallelic_snp_records=classified_biallelic_snp_records,
        gs_pass_records=gs_pass_records,
        final_matrix_variant_records=final_matrix_variant_records,
        final_matrix_sample_count=final_matrix_sample_count,
    )

    write_tsv(args.output, OUTPUT_HEADER, build_output_rows(args.cohort_id, result))
    args.summary_output.write_text(
        build_summary_text(args.cohort_id, result),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
