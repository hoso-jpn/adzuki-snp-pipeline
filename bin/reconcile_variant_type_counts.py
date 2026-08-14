#!/usr/bin/env python3
"""Reconcile raw/all record counts against the raw/snp and raw/indel splits.

`GATK_SELECTVARIANTS` produces `cohort.snp.vcf.gz` and
`cohort.indel.vcf.gz` independently from `cohort.raw.vcf.gz`. Under
GATK's `--select-type-to-include` contract, each VariantContext is
classified as a whole into exactly one overall type (SNP, INDEL,
MIXED, MNP, ... via `vc.getType()`), and is selected only on an exact
match against that single type -- confirmed empirically against the
pinned GATK 4.6.2.0 container: a MIXED-type record (one site carrying
both a SNP-type and an indel-type ALT allele) is excluded from *both*
the SNP and the INDEL selections. It is never selected into both.

A single Nextflow task in the existing seven-channel variant QC (one
`BCFTOOLS_STATS` call per raw/all, raw/snp, or raw/indel selection)
cannot answer whether the three record counts are mutually consistent,
because each of those tasks only ever sees one of the three VCFs. This
tool reads all three files together and reports:

- how many raw/all records are excluded from both type-specific
  selections (`records_not_selected`) -- expected to be >= 0 in normal
  operation, since it counts MIXED/MNP/symbolic/other-typed records,
  not double-selected ones, and
- how many records with an identical (CHROM, POS, REF, ALT) appear in
  *both* the raw/snp and raw/indel selections
  (`snp_indel_duplicate_records`) -- expected to be 0 in normal
  operation; a non-zero value is direct evidence of a duplicated
  output record, which is the only way `records_not_selected` could
  legitimately go negative.

`records_not_selected` can still be negative if that invariant is
violated for any reason; this tool always reports the value as
computed and separately flags it as a warning, rather than omitting or
clamping it, since silently hiding an unexpected negative value would
be worse than surfacing one whose cause needs investigating.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


OUTPUT_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")

# Metrics cross-referenced from the existing raw/all variant_qc.tsv
# (see bin/summarize_variant_qc.py), reported here for context on how
# MNP/other/multiallelic records relate to the snp/indel split.
CROSS_REFERENCED_METRICS: tuple[str, ...] = (
    "number_of_mnps",
    "number_of_others",
    "number_of_multiallelic_sites",
)


class MalformedVcfError(Exception):
    """Raised when a VCF cannot be reconciled safely."""


class MalformedVariantQcError(Exception):
    """Raised when the cross-referenced variant_qc.tsv is missing a metric."""


@dataclass(frozen=True)
class VcfSites:
    """The record count and per-record identity multiset read from one VCF."""

    record_count: int
    variant_keys: Counter[tuple[str, str, str, str]]


@dataclass(frozen=True)
class ReconciliationResult:
    """Every metric this tool computes."""

    raw_all_records: int
    raw_snp_records: int
    raw_indel_records: int
    records_not_selected: int
    records_not_selected_is_negative: bool
    snp_indel_duplicate_records: int
    cross_referenced: dict[str, str]


def parse_vcf_sites(path: Path) -> VcfSites:
    """Count data rows and collect a (CHROM, POS, REF, ALT) multiset from a VCF.

    Identity is keyed on all four fields, not just (CHROM, POS): a real
    SNP and a real indel can legitimately share the same coordinate
    without being the same record (or any kind of double-count), so
    position alone would misclassify that as a duplicate.
    """
    record_count = 0
    variant_keys: Counter[tuple[str, str, str, str]] = Counter()

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 5:
                raise MalformedVcfError(
                    f"{path}: line {line_number}: data row has {len(fields)} "
                    "tab-separated fields, expected at least 5"
                )

            record_count += 1
            variant_keys[(fields[0], fields[1], fields[3], fields[4])] += 1

    return VcfSites(record_count=record_count, variant_keys=variant_keys)


def read_cross_referenced_metrics(path: Path) -> dict[str, str]:
    """Read the MNP/other/multiallelic metrics from a raw/all variant_qc.tsv.

    The file is the standard (cohort_id, stage, variant_type, metric,
    value) long-format table `bin/summarize_variant_qc.py` writes; only
    the metric/value columns (the last two) are used here.
    """
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedVariantQcError(f"{path}: file is empty")

    values: dict[str, str] = {}

    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 5:
            raise MalformedVariantQcError(
                f"{path}: row has {len(fields)} tab-separated fields, expected at least 5"
            )

        values[fields[3]] = fields[4]

    missing = [metric for metric in CROSS_REFERENCED_METRICS if metric not in values]

    if missing:
        raise MalformedVariantQcError(
            f"{path}: missing required metric(s): {', '.join(missing)}"
        )

    return {metric: values[metric] for metric in CROSS_REFERENCED_METRICS}


def reconcile(
    raw_all: VcfSites,
    raw_snp: VcfSites,
    raw_indel: VcfSites,
    cross_referenced: dict[str, str],
) -> ReconciliationResult:
    """Compute the raw/all vs. raw/snp/raw/indel reconciliation."""
    records_not_selected = raw_all.record_count - raw_snp.record_count - raw_indel.record_count
    duplicate_records = sum((raw_snp.variant_keys & raw_indel.variant_keys).values())

    return ReconciliationResult(
        raw_all_records=raw_all.record_count,
        raw_snp_records=raw_snp.record_count,
        raw_indel_records=raw_indel.record_count,
        records_not_selected=records_not_selected,
        records_not_selected_is_negative=records_not_selected < 0,
        snp_indel_duplicate_records=duplicate_records,
        cross_referenced=cross_referenced,
    )


def build_output_rows(cohort_id: str, result: ReconciliationResult) -> list[list[str]]:
    """Build the reconciliation TSV data rows."""
    rows = [
        [cohort_id, "raw_all_records", str(result.raw_all_records)],
        [cohort_id, "raw_snp_records", str(result.raw_snp_records)],
        [cohort_id, "raw_indel_records", str(result.raw_indel_records)],
        [cohort_id, "records_not_selected", str(result.records_not_selected)],
        [
            cohort_id,
            "records_not_selected_is_negative",
            "true" if result.records_not_selected_is_negative else "false",
        ],
        [cohort_id, "snp_indel_duplicate_records", str(result.snp_indel_duplicate_records)],
    ]

    for metric in CROSS_REFERENCED_METRICS:
        rows.append([cohort_id, f"raw_all_{metric}", result.cross_referenced[metric]])

    return rows


def build_summary_text(cohort_id: str, result: ReconciliationResult) -> str:
    """Build the human-readable variant_type_accounting.summary.txt content."""
    lines = [
        "Variant type accounting summary",
        f"Cohort ID: {cohort_id}",
        f"raw/all records: {result.raw_all_records}",
        f"raw/snp records: {result.raw_snp_records}",
        f"raw/indel records: {result.raw_indel_records}",
        (
            f"records_not_selected = raw_all ({result.raw_all_records}) - raw_snp "
            f"({result.raw_snp_records}) - raw_indel ({result.raw_indel_records}) = "
            f"{result.records_not_selected}"
        ),
    ]

    lines.append(
        "Under GATK SelectVariants' contract, --select-type-to-include "
        "classifies each VariantContext as a whole (a single overall "
        "type: SNP, INDEL, MIXED, MNP, ...) and selects it only on an "
        "exact match. A MIXED-type record (for example one site with "
        "both a SNP and an indel ALT allele) does not match SNP and "
        "does not match INDEL, so it is excluded from *both* "
        "cohort.snp.vcf.gz and cohort.indel.vcf.gz -- it is not "
        "selected into both. records_not_selected is therefore expected "
        "to be >= 0 in normal operation, counting records excluded from "
        "both selections (MIXED, MNP, symbolic, or other non-SNP/"
        "non-INDEL types)."
    )

    if result.records_not_selected_is_negative:
        cause = (
            "which is direct evidence of duplicated output records."
            if result.snp_indel_duplicate_records > 0
            else "so duplicated records are not the cause here -- check "
            "for a wiring or input-file mismatch instead."
        )
        lines.append(
            "WARNING: records_not_selected is negative. This violates "
            "the expected raw_snp/raw_indel disjointness under GATK "
            "SelectVariants' per-VariantContext single-type "
            "classification and needs investigation. It is NOT "
            "explained by MIXED-type records being selected into both "
            "outputs (SelectVariants does not do that; see above). "
            f"snp_indel_duplicate_records = "
            f"{result.snp_indel_duplicate_records} record(s) with an "
            "identical (CHROM, POS, REF, ALT) appear in both the "
            f"raw/snp and raw/indel selections, {cause}"
        )
    else:
        lines.append(
            f"snp_indel_duplicate_records: {result.snp_indel_duplicate_records} "
            "record(s) with an identical (CHROM, POS, REF, ALT) appear in "
            "both the raw/snp and raw/indel selections (expected to be 0 "
            "in normal operation; a non-zero value indicates duplicated "
            "output records)."
        )

    lines.append(
        "SNP + indel is not guaranteed to equal raw/all: MIXED, MNP, "
        "other, and symbolic-typed records (see raw/all's own "
        "variant_qc.tsv for related counts) are excluded from both "
        "type-specific selections by GATK's classification contract."
    )
    lines.append(
        f"raw/all number_of_mnps: {result.cross_referenced['number_of_mnps']}"
    )
    lines.append(
        f"raw/all number_of_others: {result.cross_referenced['number_of_others']}"
    )
    lines.append(
        "raw/all number_of_multiallelic_sites: "
        f"{result.cross_referenced['number_of_multiallelic_sites']} (a "
        "multiallelic site whose ALT alleles are all the same "
        "elementary type, e.g. two SNP alts, is still classified as "
        "pure SNP or INDEL and is not part of records_not_selected; "
        "only MIXED-typed multiallelic sites are excluded from both "
        "selections, so this count and records_not_selected are not "
        "expected to match exactly)"
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
        description=(
            "Reconcile raw/all record counts against the raw/snp and "
            "raw/indel type-specific selections."
        )
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--raw-all-vcf", required=True, type=Path, help="Path to cohort.raw.vcf.gz."
    )
    parser.add_argument(
        "--raw-snp-vcf", required=True, type=Path, help="Path to cohort.snp.vcf.gz."
    )
    parser.add_argument(
        "--raw-indel-vcf", required=True, type=Path, help="Path to cohort.indel.vcf.gz."
    )
    parser.add_argument(
        "--raw-all-variant-qc",
        required=True,
        type=Path,
        help="Path to the existing raw/all cohort.raw.all.variant_qc.tsv.",
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
        raw_all = parse_vcf_sites(args.raw_all_vcf)
        raw_snp = parse_vcf_sites(args.raw_snp_vcf)
        raw_indel = parse_vcf_sites(args.raw_indel_vcf)
        cross_referenced = read_cross_referenced_metrics(args.raw_all_variant_qc)
    except OSError as error:
        print(f"reconcile_variant_type_counts.py: error: {error}", file=sys.stderr)
        return 1
    except (MalformedVcfError, MalformedVariantQcError) as error:
        print(f"reconcile_variant_type_counts.py: error: {error}", file=sys.stderr)
        return 1

    result = reconcile(raw_all, raw_snp, raw_indel, cross_referenced)

    write_tsv(args.output, OUTPUT_HEADER, build_output_rows(args.cohort_id, result))
    args.summary_output.write_text(
        build_summary_text(args.cohort_id, result),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
