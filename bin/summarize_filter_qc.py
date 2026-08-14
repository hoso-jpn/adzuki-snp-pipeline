#!/usr/bin/env python3
"""Summarize FILTER-value breakdown and annotation coverage for a filtered VCF.

Reads a GATK VariantFiltration output directly (the FILTER column and the
INFO annotations it was computed from) and reports two things `bcftools
stats` cannot: how FILTER values break down by exact combination and by
individual tag, and -- separately from filter *outcome* -- whether each
hard-filter annotation was even present to evaluate. A record with a
missing annotation (for example `MQRankSum` or `ReadPosRankSum`, which
GATK omits when there is no comparable ref/alt read population) is not
evaluated by the filter expression that reads it; reporting it as
"0 records failed this filter" without also reporting "N records could
not be evaluated for this filter" would misrepresent an unevaluated
record as one that passed a threshold.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path


ANNOTATION_NAMES: tuple[str, ...] = (
    "QD",
    "QUAL",
    "SOR",
    "FS",
    "MQ",
    "MQRankSum",
    "ReadPosRankSum",
)

# Which hard-filter tag (if any) reads each annotation, by variant type.
# None means this variant type has no hard filter that uses the
# annotation at all (out of scope for that type's filtering, not a
# missing-data case). Keep in sync with the filter definitions in
# workflows/adzuki_snp_pipeline.nf -- only the tag *names* are
# duplicated here, not the numeric thresholds.
FILTER_TAG_BY_VARIANT_TYPE: dict[str, dict[str, str | None]] = {
    "snp": {
        "QD": "SNP_QD_LOW",
        "QUAL": "SNP_QUAL_LOW",
        "SOR": "SNP_SOR_HIGH",
        "FS": "SNP_FS_HIGH",
        "MQ": "SNP_MQ_LOW",
        "MQRankSum": "SNP_MQRANKSUM_LOW",
        "ReadPosRankSum": "SNP_READPOSRANKSUM_LOW",
    },
    "indel": {
        "QD": "INDEL_QD_LOW",
        "QUAL": "INDEL_QUAL_LOW",
        "SOR": None,
        "FS": "INDEL_FS_HIGH",
        "MQ": None,
        "MQRankSum": None,
        "ReadPosRankSum": "INDEL_READPOSRANKSUM_LOW",
    },
}

# Values that mean "this annotation was not evaluable for this record",
# whether the VCF omits the key entirely, writes the fixed-column QUAL
# as ".", or (defensively) some tool writes a NaN-like placeholder.
MISSING_VALUE_TOKENS = frozenset({"", ".", "nan", "NaN", "NAN"})

FILTER_BREAKDOWN_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "category",
    "key",
    "record_count",
)

ANNOTATION_QC_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "annotation",
    "total_records",
    "present_records",
    "missing_records",
    "evaluable_rate",
    "filter_tag",
    "filter_tagged_records",
    "filter_hit_rate",
)

NOT_APPLICABLE = "NA"


class MalformedVcfError(Exception):
    """Raised when a VCF cannot be summarized safely."""


@dataclass(frozen=True)
class VcfRecord:
    """The subset of one VCF data row this tool needs."""

    chrom: str
    pos: str
    qual: str
    filter_value: str
    info: dict[str, str]


@dataclass(frozen=True)
class FilterBreakdown:
    """FILTER-column accounting for a filtered VCF."""

    total_records: int
    pass_records: int
    non_pass_records: int
    multi_tag_records: int
    combination_counts: dict[str, int]
    tag_counts: dict[str, int]


@dataclass(frozen=True)
class AnnotationCoverage:
    """Presence and filter-hit accounting for one annotation."""

    annotation: str
    total_records: int
    present_records: int
    missing_records: int
    evaluable_rate: str
    filter_tag: str
    filter_tagged_records: str
    filter_hit_rate: str


def parse_filtered_vcf(path: Path) -> list[VcfRecord]:
    """Parse the FILTER, QUAL, and INFO columns of a bgzipped VCF.

    Every other column (ID, REF, ALT, FORMAT, sample genotypes) is
    ignored; this tool only needs the filter outcome and the annotations
    it was computed from.
    """
    records: list[VcfRecord] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 8:
                raise MalformedVcfError(
                    f"{path}: line {line_number}: data row has {len(fields)} "
                    "tab-separated fields, expected at least 8"
                )

            records.append(
                VcfRecord(
                    chrom=fields[0],
                    pos=fields[1],
                    qual=fields[5],
                    filter_value=fields[6],
                    info=_parse_info(fields[7]),
                )
            )

    return records


def _parse_info(info_field: str) -> dict[str, str]:
    """Parse a VCF INFO field into a key/value mapping.

    Bare flags (no `=`) are kept with an empty value; none of the
    annotations this tool tracks are ever encoded as flags, so a flag
    can never be mistaken for one of them.
    """
    info: dict[str, str] = {}

    if info_field in (".", ""):
        return info

    for entry in info_field.split(";"):
        key, separator, value = entry.partition("=")
        info[key] = value if separator else ""

    return info


def _annotation_value(record: VcfRecord, annotation: str) -> str | None:
    """Look up one annotation's raw string value for a record.

    QUAL is a fixed VCF column, not an INFO key, and is handled
    separately from the other six annotations.
    """
    if annotation == "QUAL":
        return record.qual

    return record.info.get(annotation)


def _is_present(value: str | None) -> bool:
    """Return whether an annotation's raw value counts as present."""
    if value is None:
        return False

    return value not in MISSING_VALUE_TOKENS


def _filter_tags(record: VcfRecord) -> list[str]:
    """Split a record's FILTER column into its individual tag names."""
    if record.filter_value in ("PASS", ".", ""):
        return []

    return record.filter_value.split(";")


def _format_rate(numerator: int, denominator: int) -> str:
    """Format a rate to six decimal places, or "NA" if undefined."""
    if denominator <= 0:
        return NOT_APPLICABLE

    return f"{numerator / denominator:.6f}"


def compute_filter_breakdown(records: list[VcfRecord]) -> FilterBreakdown:
    """Compute FILTER-value accounting: total, PASS, tags, and combinations."""
    combination_counts: dict[str, int] = {}
    tag_counts: dict[str, int] = {}
    pass_records = 0
    multi_tag_records = 0

    for record in records:
        combination_counts[record.filter_value] = (
            combination_counts.get(record.filter_value, 0) + 1
        )

        if record.filter_value == "PASS":
            pass_records += 1
            continue

        tags = _filter_tags(record)

        if len(tags) > 1:
            multi_tag_records += 1

        for tag in tags:
            tag_counts[tag] = tag_counts.get(tag, 0) + 1

    total_records = len(records)

    return FilterBreakdown(
        total_records=total_records,
        pass_records=pass_records,
        non_pass_records=total_records - pass_records,
        multi_tag_records=multi_tag_records,
        combination_counts=combination_counts,
        tag_counts=tag_counts,
    )


def compute_annotation_coverage(
    records: list[VcfRecord],
    variant_type: str,
) -> list[AnnotationCoverage]:
    """Compute per-annotation presence and filter-hit accounting.

    `filter_hit_rate` is relative to `present_records`, not
    `total_records`: a missing annotation cannot have been evaluated by
    its filter, so it must not dilute the rate for records that could be.
    """
    tag_by_annotation = FILTER_TAG_BY_VARIANT_TYPE[variant_type]
    total_records = len(records)
    coverages: list[AnnotationCoverage] = []

    for annotation in ANNOTATION_NAMES:
        present_records = sum(
            1
            for record in records
            if _is_present(_annotation_value(record, annotation))
        )
        missing_records = total_records - present_records
        evaluable_rate = _format_rate(present_records, total_records)
        filter_tag = tag_by_annotation[annotation]

        if filter_tag is None:
            coverages.append(
                AnnotationCoverage(
                    annotation=annotation,
                    total_records=total_records,
                    present_records=present_records,
                    missing_records=missing_records,
                    evaluable_rate=evaluable_rate,
                    filter_tag=NOT_APPLICABLE,
                    filter_tagged_records=NOT_APPLICABLE,
                    filter_hit_rate=NOT_APPLICABLE,
                )
            )
            continue

        filter_tagged_records = sum(
            1 for record in records if filter_tag in _filter_tags(record)
        )

        coverages.append(
            AnnotationCoverage(
                annotation=annotation,
                total_records=total_records,
                present_records=present_records,
                missing_records=missing_records,
                evaluable_rate=evaluable_rate,
                filter_tag=filter_tag,
                filter_tagged_records=str(filter_tagged_records),
                filter_hit_rate=_format_rate(filter_tagged_records, present_records),
            )
        )

    return coverages


def build_filter_breakdown_rows(
    cohort_id: str,
    stage: str,
    variant_type: str,
    breakdown: FilterBreakdown,
) -> list[list[str]]:
    """Build filter_breakdown.tsv data rows: summary, combination, and tag rows."""
    rows: list[list[str]] = [
        [cohort_id, stage, variant_type, "summary", "total_records", str(breakdown.total_records)],
        [cohort_id, stage, variant_type, "summary", "pass_records", str(breakdown.pass_records)],
        [cohort_id, stage, variant_type, "summary", "non_pass_records", str(breakdown.non_pass_records)],
        [cohort_id, stage, variant_type, "summary", "multi_tag_records", str(breakdown.multi_tag_records)],
        [
            cohort_id,
            stage,
            variant_type,
            "summary",
            "distinct_filter_combinations",
            str(len(breakdown.combination_counts)),
        ],
        [
            cohort_id,
            stage,
            variant_type,
            "summary",
            "distinct_filter_tags",
            str(len(breakdown.tag_counts)),
        ],
    ]

    for combination in sorted(breakdown.combination_counts):
        rows.append(
            [
                cohort_id,
                stage,
                variant_type,
                "combination",
                combination,
                str(breakdown.combination_counts[combination]),
            ]
        )

    for tag in sorted(breakdown.tag_counts):
        rows.append(
            [cohort_id, stage, variant_type, "tag", tag, str(breakdown.tag_counts[tag])]
        )

    return rows


def build_annotation_qc_rows(
    cohort_id: str,
    stage: str,
    variant_type: str,
    coverages: list[AnnotationCoverage],
) -> list[list[str]]:
    """Build annotation_qc.tsv data rows, one per tracked annotation."""
    return [
        [
            cohort_id,
            stage,
            variant_type,
            coverage.annotation,
            str(coverage.total_records),
            str(coverage.present_records),
            str(coverage.missing_records),
            coverage.evaluable_rate,
            coverage.filter_tag,
            coverage.filter_tagged_records,
            coverage.filter_hit_rate,
        ]
        for coverage in coverages
    ]


def build_filter_qc_summary_text(
    cohort_id: str,
    stage: str,
    variant_type: str,
    breakdown: FilterBreakdown,
    coverages: list[AnnotationCoverage],
) -> str:
    """Build the human-readable filter_qc.summary.txt content."""
    lines = [
        "Filter QC summary",
        f"Cohort ID: {cohort_id}",
        f"Stage: {stage}",
        f"Variant type: {variant_type}",
        f"Total records: {breakdown.total_records}",
        f"PASS records: {breakdown.pass_records}",
        f"Non-PASS records: {breakdown.non_pass_records}",
        (
            f"Reconciliation: total ({breakdown.total_records}) = PASS "
            f"({breakdown.pass_records}) + non-PASS ({breakdown.non_pass_records})"
        ),
        f"Multi-tag records: {breakdown.multi_tag_records}",
        (
            "Note: the sum of per-tag record counts can exceed non-PASS "
            "records whenever records carry more than one filter tag."
        ),
        "",
        (
            "Annotation coverage (present/missing does not by itself mean "
            "a filter passed or failed: a missing annotation means GATK "
            "could not evaluate that filter expression for the record, "
            "not that the record satisfied the threshold):"
        ),
    ]

    for coverage in coverages:
        if coverage.filter_tag == NOT_APPLICABLE:
            applicability = "no hard filter uses this annotation for this variant type"
        else:
            applicability = (
                f"filter={coverage.filter_tag}, tagged={coverage.filter_tagged_records}, "
                f"hit_rate_among_present={coverage.filter_hit_rate}"
            )

        lines.append(
            f"  {coverage.annotation}: present {coverage.present_records}/"
            f"{coverage.total_records} (evaluable_rate={coverage.evaluable_rate}), "
            f"{applicability}"
        )

    return "\n".join(lines) + "\n"


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the filter QC summarizer CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Summarize FILTER-value breakdown and annotation coverage for "
            "a GATK VariantFiltration output."
        )
    )
    parser.add_argument(
        "--filtered-vcf",
        required=True,
        type=Path,
        help="Path to a bgzipped VariantFiltration output VCF.",
    )
    parser.add_argument(
        "--cohort-id",
        required=True,
        help="Cohort identifier recorded in every output row.",
    )
    parser.add_argument(
        "--stage",
        required=True,
        help="QC stage label (expected: filtered).",
    )
    parser.add_argument(
        "--variant-type",
        required=True,
        choices=sorted(FILTER_TAG_BY_VARIANT_TYPE),
        help="Variant type label; selects which hard filters apply.",
    )
    parser.add_argument(
        "--filter-breakdown-output",
        required=True,
        type=Path,
        help="Output path for the FILTER-value breakdown TSV.",
    )
    parser.add_argument(
        "--annotation-qc-output",
        required=True,
        type=Path,
        help="Output path for the annotation coverage TSV.",
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
        records = parse_filtered_vcf(args.filtered_vcf)
    except OSError as error:
        print(
            f"summarize_filter_qc.py: error: cannot read {args.filtered_vcf}: {error}",
            file=sys.stderr,
        )
        return 1
    except MalformedVcfError as error:
        print(f"summarize_filter_qc.py: error: {error}", file=sys.stderr)
        return 1

    breakdown = compute_filter_breakdown(records)
    coverages = compute_annotation_coverage(records, args.variant_type)

    write_tsv(
        args.filter_breakdown_output,
        FILTER_BREAKDOWN_HEADER,
        build_filter_breakdown_rows(args.cohort_id, args.stage, args.variant_type, breakdown),
    )
    write_tsv(
        args.annotation_qc_output,
        ANNOTATION_QC_HEADER,
        build_annotation_qc_rows(args.cohort_id, args.stage, args.variant_type, coverages),
    )
    args.summary_output.write_text(
        build_filter_qc_summary_text(
            args.cohort_id,
            args.stage,
            args.variant_type,
            breakdown,
            coverages,
        ),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
