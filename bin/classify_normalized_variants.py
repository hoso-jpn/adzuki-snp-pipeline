#!/usr/bin/env python3
"""Reclassify bcftools-norm output into a biallelic-SNP-only GS panel input.

`bcftools norm -m-` splits every multiallelic record (including a
GATK-MIXED record, one site with both a SNP-type and an indel-type ALT
allele -- see `bin/reconcile_variant_type_counts.py` for how GATK's own
type label works pre-split) into one biallelic row per original ALT
allele. Once split, GATK's pre-split type label no longer applies to
any individual output row, so this tool reclassifies every row from
scratch by its own REF/ALT shape:

- `snp`: REF and ALT are each a single base.
- `mnp`: REF and ALT are equal length, longer than one base.
- `indel`: REF and ALT differ in length.
- `symbolic_or_star`: ALT contains symbolic/breakend syntax (`<`, `[`,
  `]`) or is the spanning-deletion marker `*`.

Only `snp`-classified rows are eligible for the GS panel. A row whose
(CHROM, POS, REF, ALT) key is not unique after classification is
excluded entirely -- every occurrence, not just the extras -- because
there is no automatic way to decide which of two colliding records is
correct; keeping one at random would be a silent, unreviewable choice.

This tool also resets every kept record's FILTER to "." (not "PASS"):
bcftools norm propagates the *original* record's FILTER value to every
split child regardless of that child's new shape, which can leave a
semantically mismatched tag (e.g. a SNP-named filter on a row that is
now shaped like an indel). Resetting to "." marks these records as not
yet evaluated, matching this repository's own `raw` stage convention,
so the GS-specific hard-filter step downstream makes its own PASS/FAIL
decision from a clean slate.

The output is a plain-text (uncompressed) VCF, not `.vcf.gz`: this
script only ever needs to *read* bgzip-compressed input (BGZF is valid
gzip and decompresses correctly with the standard library, the same as
every other script in this repository), but the standard `gzip` module
writes plain DEFLATE, which `tabix` cannot index. Compressing and
indexing the output is left to a dedicated bcftools-container step.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

VARIANT_CLASSES: tuple[str, ...] = ("snp", "mnp", "indel", "symbolic_or_star")
ELIGIBLE_CLASS = "snp"
RESET_FILTER_VALUE = "."
SYMBOLIC_MARKERS = ("<", "[", "]")
SPANNING_DELETION_ALLELE = "*"

ACCOUNTING_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")


class MalformedVcfError(Exception):
    """Raised when a normalized VCF cannot be classified safely."""


@dataclass(frozen=True)
class VcfHeader:
    """The header lines and sample columns read from a VCF."""

    meta_lines: tuple[str, ...]
    chrom_line: str
    sample_names: tuple[str, ...]


@dataclass(frozen=True)
class ClassifiedRecord:
    """One post-split data row, classified by its own REF/ALT shape."""

    chrom: str
    pos: str
    id_: str
    ref: str
    alt: str
    qual: str
    info: str
    format_: str
    sample_fields: tuple[str, ...]
    variant_class: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.chrom, self.pos, self.ref, self.alt)


@dataclass(frozen=True)
class ClassificationResult:
    """Every record and count this tool computes."""

    header: VcfHeader
    total_input_records: int
    class_counts: dict[str, int]
    duplicate_key_records: int
    distinct_duplicate_keys: tuple[tuple[str, str, str, str], ...]
    output_records: tuple[ClassifiedRecord, ...]


def classify_variant(ref: str, alt: str) -> str:
    """Classify one post-split (single-ALT) record by REF/ALT shape."""
    if any(marker in alt for marker in SYMBOLIC_MARKERS) or alt == SPANNING_DELETION_ALLELE:
        return "symbolic_or_star"

    if "," in alt:
        raise MalformedVcfError(
            f"record {ref}>{alt} still has multiple ALT alleles after "
            "bcftools norm -m- splitting; expected exactly one"
        )

    if len(ref) == len(alt):
        return "snp" if len(ref) == 1 else "mnp"

    return "indel"


def _parse_header(lines: list[str]) -> VcfHeader:
    meta_lines: list[str] = []
    chrom_line: str | None = None

    for line in lines:
        if line.startswith("##"):
            if line.startswith("##FILTER="):
                continue
            meta_lines.append(line)
        elif line.startswith("#CHROM"):
            chrom_line = line
            break

    if chrom_line is None:
        raise MalformedVcfError("no #CHROM header line found")

    fields = chrom_line.split("\t")
    if len(fields) <= 9:
        raise MalformedVcfError(
            f"#CHROM header has {len(fields)} fields, expected at least 10 "
            "(9 fixed columns plus one or more samples)"
        )

    return VcfHeader(
        meta_lines=tuple(meta_lines),
        chrom_line=chrom_line,
        sample_names=tuple(fields[9:]),
    )


def parse_normalized_vcf(path: Path) -> ClassificationResult:
    """Read a bcftools-norm output VCF and classify every data row."""
    header: VcfHeader | None = None
    total_input_records = 0
    class_counts: dict[str, int] = {name: 0 for name in VARIANT_CLASSES}
    key_counts: Counter[tuple[str, str, str, str]] = Counter()
    candidates: list[ClassifiedRecord] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_lines: list[str] = []

        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("#CHROM"):
                    header = _parse_header(header_lines)
                continue

            if header is None:
                raise MalformedVcfError(f"{path}: data row seen before #CHROM header")

            fields = line.split("\t")
            if len(fields) < 10:
                raise MalformedVcfError(
                    f"{path}: data row has {len(fields)} tab-separated fields, "
                    "expected at least 10"
                )

            total_input_records += 1
            chrom, pos, id_, ref, alt, qual, _filter, info, format_ = fields[:9]
            sample_fields = tuple(fields[9:])

            variant_class = classify_variant(ref, alt)
            class_counts[variant_class] += 1

            if variant_class != ELIGIBLE_CLASS:
                continue

            record = ClassifiedRecord(
                chrom=chrom,
                pos=pos,
                id_=id_,
                ref=ref,
                alt=alt,
                qual=qual,
                info=info,
                format_=format_,
                sample_fields=sample_fields,
                variant_class=variant_class,
            )
            key_counts[record.key] += 1
            candidates.append(record)

    if header is None:
        raise MalformedVcfError(f"{path}: no #CHROM header line found")

    distinct_duplicate_keys = tuple(
        sorted(key for key, count in key_counts.items() if count > 1)
    )
    duplicate_key_records = sum(
        key_counts[key] for key in distinct_duplicate_keys
    )
    output_records = tuple(
        record for record in candidates if key_counts[record.key] == 1
    )

    return ClassificationResult(
        header=header,
        total_input_records=total_input_records,
        class_counts=class_counts,
        duplicate_key_records=duplicate_key_records,
        distinct_duplicate_keys=distinct_duplicate_keys,
        output_records=output_records,
    )


def build_accounting_rows(cohort_id: str, result: ClassificationResult) -> list[list[str]]:
    """Build the normalization/classification accounting TSV data rows."""
    rows = [[cohort_id, "total_input_records", str(result.total_input_records)]]

    for variant_class in VARIANT_CLASSES:
        rows.append(
            [cohort_id, f"{variant_class}_records", str(result.class_counts[variant_class])]
        )

    rows.append([cohort_id, "duplicate_key_records", str(result.duplicate_key_records)])
    rows.append(
        [cohort_id, "distinct_duplicate_keys", str(len(result.distinct_duplicate_keys))]
    )
    rows.append([cohort_id, "output_records", str(len(result.output_records))])

    return rows


def build_summary_text(cohort_id: str, result: ClassificationResult) -> str:
    """Build the human-readable normalization/classification summary."""
    lines = [
        "Variant normalization and classification summary",
        f"Cohort ID: {cohort_id}",
        f"Total input records (post bcftools norm -m- splitting): {result.total_input_records}",
    ]

    for variant_class in VARIANT_CLASSES:
        lines.append(f"  {variant_class}: {result.class_counts[variant_class]}")

    lines.append(
        "Only 'snp'-classified records are eligible for the GS panel; "
        "mnp/indel/symbolic_or_star records are excluded here, not "
        "silently dropped further downstream."
    )
    lines.append(
        f"Duplicate (CHROM, POS, REF, ALT) keys: {len(result.distinct_duplicate_keys)} "
        f"distinct key(s), {result.duplicate_key_records} record(s) total -- "
        "every occurrence of a colliding key is excluded, not just the extras, "
        "since there is no automatic way to decide which one is correct."
    )

    if result.distinct_duplicate_keys:
        lines.append("Colliding keys:")
        for chrom, pos, ref, alt in result.distinct_duplicate_keys:
            lines.append(f"  {chrom}:{pos} {ref}>{alt}")

    lines.append(
        f"Output records (eligible for GS hard-filtering): {len(result.output_records)}"
    )
    lines.append(
        "FILTER has been reset to '.' on every output record: bcftools norm "
        "propagates the original record's FILTER value to every split child "
        "regardless of that child's new shape, which can leave a "
        "semantically mismatched tag; the GS-specific hard-filter step "
        "downstream makes its own PASS/FAIL decision from a clean slate."
    )

    return "\n".join(lines) + "\n"


def write_vcf(path: Path, result: ClassificationResult) -> None:
    """Write the eligible, FILTER-reset records as a plain-text VCF."""
    lines = list(result.header.meta_lines)
    lines.append(result.header.chrom_line)

    for record in result.output_records:
        fields = [
            record.chrom,
            record.pos,
            record.id_,
            record.ref,
            record.alt,
            record.qual,
            RESET_FILTER_VALUE,
            record.info,
            record.format_,
            *record.sample_fields,
        ]
        lines.append("\t".join(fields))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the classification CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Reclassify bcftools-norm output by post-split REF/ALT shape "
            "and select the biallelic-SNP-only subset eligible for a GS panel."
        )
    )
    parser.add_argument(
        "--normalized-vcf",
        required=True,
        type=Path,
        help="Path to a bgzipped bcftools norm -m- output VCF.",
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the plain-text, biallelic-SNP-only VCF.",
    )
    parser.add_argument(
        "--accounting-output",
        required=True,
        type=Path,
        help="Output path for the classification accounting TSV.",
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
        result = parse_normalized_vcf(args.normalized_vcf)
    except OSError as error:
        print(
            f"classify_normalized_variants.py: error: cannot read {args.normalized_vcf}: {error}",
            file=sys.stderr,
        )
        return 1
    except MalformedVcfError as error:
        print(f"classify_normalized_variants.py: error: {error}", file=sys.stderr)
        return 1

    write_vcf(args.output, result)
    write_tsv(
        args.accounting_output,
        ACCOUNTING_HEADER,
        build_accounting_rows(args.cohort_id, result),
    )
    args.summary_output.write_text(
        build_summary_text(args.cohort_id, result),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
