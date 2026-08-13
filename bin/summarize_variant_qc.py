#!/usr/bin/env python3
"""Summarize a `bcftools stats` report into the pipeline's variant QC contract.

Reads the plain-text output of `bcftools stats` and produces the three
downstream QC artifacts the workflow publishes for every cohort/stage/
variant-type combination: a machine-readable cohort-level QC table, a
machine-readable per-sample QC table, and a human-readable summary.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path


REQUIRED_SUMMARY_KEYS: tuple[str, ...] = (
    "number of samples:",
    "number of records:",
    "number of SNPs:",
    "number of MNPs:",
    "number of indels:",
    "number of others:",
    "number of multiallelic sites:",
)

VARIANT_QC_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "metric",
    "value",
)

SAMPLE_QC_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "sample",
    "reference_homozygous",
    "non_reference_homozygous",
    "heterozygous",
    "missing",
    "missingness_rate",
    "average_depth",
    "singletons",
)


class MalformedBcftoolsStatsError(Exception):
    """Raised when a `bcftools stats` report cannot be summarized safely."""


@dataclass(frozen=True)
class SummaryNumbers:
    """Cohort-wide counts parsed from `bcftools stats` SN rows."""

    number_of_samples: int
    number_of_records: int
    number_of_snps: int
    number_of_mnps: int
    number_of_indels: int
    number_of_others: int
    number_of_multiallelic_sites: int


@dataclass(frozen=True)
class TransitionTransversion:
    """Transition/transversion counts parsed from a `bcftools stats` TSTV row."""

    transitions: int
    transversions: int
    ratio: str


@dataclass(frozen=True)
class PerSampleCounts:
    """Per-sample genotype counts parsed from one `bcftools stats` PSC row."""

    sample: str
    reference_homozygous: int
    non_reference_homozygous: int
    heterozygous: int
    average_depth: str
    singletons: int
    missing: int


@dataclass(frozen=True)
class ParsedStats:
    """Everything this tool needs from one `bcftools stats` report."""

    summary: SummaryNumbers
    transition_transversion: TransitionTransversion | None
    per_sample: tuple[PerSampleCounts, ...]


def parse_bcftools_stats(text: str, source: Path) -> ParsedStats:
    """Parse the SN, TSTV, and PSC sections of a `bcftools stats` report.

    Every other record type (ID, SiS, AF, QUAL, ST, DP, PSI, HWE, VAF, and
    `#`-prefixed comment/header lines) is ignored: `bcftools stats` always
    emits sections this tool does not need, and their presence is not an
    error. Raises MalformedBcftoolsStatsError, naming `source` and the
    offending line, when a row this tool does rely on is truncated or
    holds a non-integer value where an integer is required.
    """
    summary_values: dict[str, str] = {}
    transition_transversion: TransitionTransversion | None = None
    per_sample: list[PerSampleCounts] = []

    for line_number, line in enumerate(text.splitlines(), start=1):
        fields = line.split("\t")
        record_type = fields[0]

        if record_type == "SN":
            _record_summary_value(fields, summary_values, source, line_number)
        elif record_type == "TSTV":
            transition_transversion = _parse_transition_transversion(
                fields,
                source,
                line_number,
            )
        elif record_type == "PSC":
            per_sample.append(_parse_per_sample_counts(fields, source, line_number))

    summary = _build_summary_numbers(summary_values, source)

    return ParsedStats(
        summary=summary,
        transition_transversion=transition_transversion,
        per_sample=tuple(per_sample),
    )


def _record_summary_value(
    fields: list[str],
    summary_values: dict[str, str],
    source: Path,
    line_number: int,
) -> None:
    """Record one SN row's key/value pair for later validation."""
    if len(fields) < 4:
        raise MalformedBcftoolsStatsError(
            f"{source}: line {line_number}: SN row has {len(fields)} "
            "tab-separated fields, expected at least 4"
        )

    key, value = fields[2], fields[3]
    summary_values[key] = value


def _build_summary_numbers(
    summary_values: dict[str, str],
    source: Path,
) -> SummaryNumbers:
    """Validate and convert the SN section into a SummaryNumbers.

    Every key in REQUIRED_SUMMARY_KEYS is present in real `bcftools
    stats` output regardless of record count, so a missing key means the
    input is truncated or otherwise malformed rather than legitimately
    empty; this is reported as an error instead of defaulting to zero.
    """
    missing_keys = [key for key in REQUIRED_SUMMARY_KEYS if key not in summary_values]

    if missing_keys:
        raise MalformedBcftoolsStatsError(
            f"{source}: missing required SN metric(s): {', '.join(missing_keys)}"
        )

    parsed_values: dict[str, int] = {}

    for key in REQUIRED_SUMMARY_KEYS:
        raw_value = summary_values[key]

        try:
            parsed_values[key] = int(raw_value)
        except ValueError as error:
            raise MalformedBcftoolsStatsError(
                f"{source}: SN metric {key!r} has non-integer value {raw_value!r}"
            ) from error

    return SummaryNumbers(
        number_of_samples=parsed_values["number of samples:"],
        number_of_records=parsed_values["number of records:"],
        number_of_snps=parsed_values["number of SNPs:"],
        number_of_mnps=parsed_values["number of MNPs:"],
        number_of_indels=parsed_values["number of indels:"],
        number_of_others=parsed_values["number of others:"],
        number_of_multiallelic_sites=parsed_values["number of multiallelic sites:"],
    )


def _parse_transition_transversion(
    fields: list[str],
    source: Path,
    line_number: int,
) -> TransitionTransversion:
    """Parse one TSTV row. The ratio field is kept as-is (for example "0.00")."""
    if len(fields) < 5:
        raise MalformedBcftoolsStatsError(
            f"{source}: line {line_number}: TSTV row has {len(fields)} "
            "tab-separated fields, expected at least 5"
        )

    try:
        transitions = int(fields[2])
        transversions = int(fields[3])
    except ValueError as error:
        raise MalformedBcftoolsStatsError(
            f"{source}: line {line_number}: TSTV row has non-integer "
            f"transitions/transversions: {fields[2]!r}, {fields[3]!r}"
        ) from error

    return TransitionTransversion(
        transitions=transitions,
        transversions=transversions,
        ratio=fields[4],
    )


def _parse_per_sample_counts(
    fields: list[str],
    source: Path,
    line_number: int,
) -> PerSampleCounts:
    """Parse one PSC row into its per-sample genotype counts."""
    if len(fields) < 14:
        raise MalformedBcftoolsStatsError(
            f"{source}: line {line_number}: PSC row has {len(fields)} "
            "tab-separated fields, expected at least 14"
        )

    sample = fields[2]

    try:
        reference_homozygous = int(fields[3])
        non_reference_homozygous = int(fields[4])
        heterozygous = int(fields[5])
        singletons = int(fields[10])
        missing = int(fields[13])
    except ValueError as error:
        raise MalformedBcftoolsStatsError(
            f"{source}: line {line_number}: PSC row for sample {sample!r} "
            "has a non-integer genotype count"
        ) from error

    return PerSampleCounts(
        sample=sample,
        reference_homozygous=reference_homozygous,
        non_reference_homozygous=non_reference_homozygous,
        heterozygous=heterozygous,
        average_depth=fields[9],
        singletons=singletons,
        missing=missing,
    )


def _format_missingness_rate(missing: int, denominator: int) -> str:
    """Format a missingness rate to six decimal places, or "NA" if undefined."""
    if denominator <= 0:
        return "NA"

    return f"{missing / denominator:.6f}"


def build_variant_qc_rows(
    cohort_id: str,
    stage: str,
    variant_type: str,
    parsed: ParsedStats,
) -> list[list[str]]:
    """Build the cohort-level variant_qc.tsv data rows (excluding the header)."""
    summary = parsed.summary
    tstv = parsed.transition_transversion
    transitions = tstv.transitions if tstv is not None else 0
    transversions = tstv.transversions if tstv is not None else 0
    ratio = tstv.ratio if tstv is not None else "NA"

    cohort_missing_genotypes = sum(sample.missing for sample in parsed.per_sample)
    cohort_total_genotypes = summary.number_of_records * summary.number_of_samples
    cohort_missingness_rate = _format_missingness_rate(
        cohort_missing_genotypes,
        cohort_total_genotypes,
    )
    sample_names = ",".join(sample.sample for sample in parsed.per_sample)

    metrics: list[tuple[str, str]] = [
        ("number_of_samples", str(summary.number_of_samples)),
        ("number_of_records", str(summary.number_of_records)),
        ("number_of_snps", str(summary.number_of_snps)),
        ("number_of_mnps", str(summary.number_of_mnps)),
        ("number_of_indels", str(summary.number_of_indels)),
        ("number_of_others", str(summary.number_of_others)),
        ("number_of_multiallelic_sites", str(summary.number_of_multiallelic_sites)),
        ("transitions", str(transitions)),
        ("transversions", str(transversions)),
        ("transition_transversion_ratio", ratio),
        ("cohort_missing_genotypes", str(cohort_missing_genotypes)),
        ("cohort_total_genotypes", str(cohort_total_genotypes)),
        ("cohort_missingness_rate", cohort_missingness_rate),
        ("sample_names", sample_names),
    ]

    return [[cohort_id, stage, variant_type, metric, value] for metric, value in metrics]


def build_sample_qc_rows(
    cohort_id: str,
    stage: str,
    variant_type: str,
    parsed: ParsedStats,
) -> list[list[str]]:
    """Build the per-sample sample_qc.tsv data rows, in PSC row order."""
    rows: list[list[str]] = []

    for sample in parsed.per_sample:
        missingness_rate = _format_missingness_rate(
            sample.missing,
            parsed.summary.number_of_records,
        )
        rows.append(
            [
                cohort_id,
                stage,
                variant_type,
                sample.sample,
                str(sample.reference_homozygous),
                str(sample.non_reference_homozygous),
                str(sample.heterozygous),
                str(sample.missing),
                missingness_rate,
                sample.average_depth,
                str(sample.singletons),
            ]
        )

    return rows


def build_summary_text(
    cohort_id: str,
    stage: str,
    variant_type: str,
    parsed: ParsedStats,
) -> str:
    """Build the human-readable summary.txt content."""
    summary = parsed.summary
    tstv = parsed.transition_transversion
    transitions = tstv.transitions if tstv is not None else 0
    transversions = tstv.transversions if tstv is not None else 0
    ratio = tstv.ratio if tstv is not None else "NA"

    cohort_missing_genotypes = sum(sample.missing for sample in parsed.per_sample)
    cohort_total_genotypes = summary.number_of_records * summary.number_of_samples
    cohort_missingness_rate = _format_missingness_rate(
        cohort_missing_genotypes,
        cohort_total_genotypes,
    )
    sample_names = ",".join(sample.sample for sample in parsed.per_sample)

    lines = [
        "Variant QC summary",
        f"Cohort ID: {cohort_id}",
        f"Stage: {stage}",
        f"Variant type: {variant_type}",
        f"Samples ({summary.number_of_samples}): {sample_names}",
        f"Records: {summary.number_of_records}",
        f"SNPs: {summary.number_of_snps}",
        f"MNPs: {summary.number_of_mnps}",
        f"Indels: {summary.number_of_indels}",
        f"Other variants: {summary.number_of_others}",
        f"Multiallelic sites: {summary.number_of_multiallelic_sites}",
        f"Transitions: {transitions}",
        f"Transversions: {transversions}",
        f"Ti/Tv ratio: {ratio}",
        (
            f"Missing genotypes: {cohort_missing_genotypes}/{cohort_total_genotypes} "
            f"({cohort_missingness_rate})"
        ),
    ]

    return "\n".join(lines) + "\n"


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the summarizer CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Summarize a bcftools stats report into the pipeline's "
            "variant_qc.tsv, sample_qc.tsv, and summary.txt contract."
        )
    )
    parser.add_argument(
        "--bcftools-stats",
        required=True,
        type=Path,
        help="Path to a bcftools stats report (bcftools stats output).",
    )
    parser.add_argument(
        "--cohort-id",
        required=True,
        help="Cohort identifier recorded in every output row.",
    )
    parser.add_argument(
        "--stage",
        required=True,
        help="QC stage label (for example raw, filtered, or pass).",
    )
    parser.add_argument(
        "--variant-type",
        required=True,
        help="Variant type label (for example all, snp, or indel).",
    )
    parser.add_argument(
        "--variant-qc-output",
        required=True,
        type=Path,
        help="Output path for the cohort-level variant QC TSV.",
    )
    parser.add_argument(
        "--sample-qc-output",
        required=True,
        type=Path,
        help="Output path for the per-sample QC TSV.",
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
        text = args.bcftools_stats.read_text(encoding="utf-8")
    except OSError as error:
        print(
            f"summarize_variant_qc.py: error: cannot read {args.bcftools_stats}: {error}",
            file=sys.stderr,
        )
        return 1

    try:
        parsed = parse_bcftools_stats(text, args.bcftools_stats)
    except MalformedBcftoolsStatsError as error:
        print(f"summarize_variant_qc.py: error: {error}", file=sys.stderr)
        return 1

    variant_qc_rows = build_variant_qc_rows(
        args.cohort_id,
        args.stage,
        args.variant_type,
        parsed,
    )
    sample_qc_rows = build_sample_qc_rows(
        args.cohort_id,
        args.stage,
        args.variant_type,
        parsed,
    )
    summary_text = build_summary_text(
        args.cohort_id,
        args.stage,
        args.variant_type,
        parsed,
    )

    write_tsv(args.variant_qc_output, VARIANT_QC_HEADER, variant_qc_rows)
    write_tsv(args.sample_qc_output, SAMPLE_QC_HEADER, sample_qc_rows)
    args.summary_output.write_text(summary_text, encoding="utf-8")

    return 0


if __name__ == "__main__":
    sys.exit(main())
