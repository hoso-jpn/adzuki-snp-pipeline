#!/usr/bin/env python3
"""Reconcile raw/all record counts against the raw/snp and raw/indel splits.

`GATK_SELECTVARIANTS` produces `cohort.snp.vcf.gz` and
`cohort.indel.vcf.gz` independently from `cohort.raw.vcf.gz`, each by
asking whether *any* allele at a site matches the requested type. A
single Nextflow task in the existing seven-channel variant QC (one
`BCFTOOLS_STATS` call per raw/all, raw/snp, or raw/indel selection)
cannot answer whether the three record counts are mutually consistent,
because each of those tasks only ever sees one of the three VCFs. This
tool reads all three files together and reports:

- how many raw/all records are not present in either type-specific
  selection (`records_not_selected`), and
- how many (contig, position) sites appear in *both* the snp and indel
  selections (`snp_indel_overlap_records`), the direct evidence for a
  multiallelic site with both a SNP and an indel allele being counted
  in both outputs.

`records_not_selected` can be negative when overlap inflates
`raw_snp_records + raw_indel_records` past `raw_all_records`; this tool
always reports the value as computed and separately flags it, rather
than omitting or clamping it.
"""

from __future__ import annotations

import argparse
import gzip
import sys
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
    """The record count and (contig, position) set read from one VCF."""

    record_count: int
    positions: frozenset[tuple[str, str]]


@dataclass(frozen=True)
class ReconciliationResult:
    """Every metric this tool computes."""

    raw_all_records: int
    raw_snp_records: int
    raw_indel_records: int
    records_not_selected: int
    records_not_selected_is_negative: bool
    snp_indel_overlap_records: int
    cross_referenced: dict[str, str]


def parse_vcf_sites(path: Path) -> VcfSites:
    """Count data rows and collect (contig, position) pairs from a VCF."""
    record_count = 0
    positions: set[tuple[str, str]] = set()

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                raise MalformedVcfError(
                    f"{path}: line {line_number}: data row has {len(fields)} "
                    "tab-separated fields, expected at least 2"
                )

            record_count += 1
            positions.add((fields[0], fields[1]))

    return VcfSites(record_count=record_count, positions=frozenset(positions))


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
    overlap = len(raw_snp.positions & raw_indel.positions)

    return ReconciliationResult(
        raw_all_records=raw_all.record_count,
        raw_snp_records=raw_snp.record_count,
        raw_indel_records=raw_indel.record_count,
        records_not_selected=records_not_selected,
        records_not_selected_is_negative=records_not_selected < 0,
        snp_indel_overlap_records=overlap,
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
        [cohort_id, "snp_indel_overlap_records", str(result.snp_indel_overlap_records)],
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

    if result.records_not_selected_is_negative:
        lines.append(
            "WARNING: records_not_selected is negative. GATK SelectVariants "
            "selects a record into cohort.snp.vcf.gz or cohort.indel.vcf.gz "
            "whenever *any* allele at that site matches the requested type, "
            "so a multiallelic site carrying both a SNP and an indel allele "
            "is counted in both outputs and inflates raw_snp + raw_indel "
            f"past raw_all. snp_indel_overlap_records = "
            f"{result.snp_indel_overlap_records} site(s) appear in both "
            "the raw/snp and raw/indel selections, which is direct "
            "evidence for this double-counting."
        )
    else:
        lines.append(
            f"snp_indel_overlap_records: {result.snp_indel_overlap_records} "
            "site(s) appear in both the raw/snp and raw/indel selections."
        )

    lines.append(
        "SNP + indel is not guaranteed to equal raw/all: MNP, other, and "
        "multiallelic records (from raw/all's own variant_qc.tsv) may "
        "fall outside both type-specific selections, or (see above) be "
        "counted in both."
    )
    lines.append(
        f"raw/all number_of_mnps: {result.cross_referenced['number_of_mnps']}"
    )
    lines.append(
        f"raw/all number_of_others: {result.cross_referenced['number_of_others']}"
    )
    lines.append(
        "raw/all number_of_multiallelic_sites: "
        f"{result.cross_referenced['number_of_multiallelic_sites']}"
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
