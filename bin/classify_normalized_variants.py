#!/usr/bin/env python3
"""Locus-stream normalized VCF records into the biallelic-SNP GS input."""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

VARIANT_CLASSES = ("snp", "mnp", "indel", "symbolic_or_star", "no_alt")
ELIGIBLE_CLASS = "snp"
RESET_FILTER_VALUE = "."
ACCOUNTING_HEADER = ("cohort_id", "metric", "value")


class MalformedVcfError(Exception):
    """Raised when a normalized VCF cannot be classified safely."""


@dataclass
class ClassificationStats:
    total_input_records: int
    class_counts: dict[str, int]
    duplicate_key_records: int = 0
    distinct_duplicate_keys: int = 0
    output_records: int = 0


@dataclass(frozen=True)
class LocusRecord:
    ref: str
    alt: str
    output_line: str

    @property
    def allele_key(self) -> tuple[str, str]:
        return self.ref, self.alt


def classify_variant(ref: str, alt: str) -> str:
    """Classify one post-split (single-ALT) record by REF/ALT shape."""
    if alt == ".":
        return "no_alt"
    if alt == "*" or any(marker in alt for marker in ("<", "[", "]")):
        return "symbolic_or_star"
    if "," in alt:
        raise MalformedVcfError(
            f"record {ref}>{alt} still has multiple ALT alleles after "
            "bcftools norm -m- splitting; expected exactly one"
        )
    if len(ref) == len(alt):
        return "snp" if len(ref) == 1 else "mnp"
    return "indel"


def _temporary_path(final_path: Path) -> Path:
    final_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{final_path.name}.", suffix=".tmp", dir=final_path.parent
    )
    os.close(descriptor)
    # Nextflow may run the container as a different uid from the publishing
    # process. Match ordinary pipeline artifacts (umask 022) after mkstemp's
    # deliberately private creation mode so publishDir can copy finalized files.
    os.chmod(name, 0o644)
    return Path(name)


def _flush_locus(chrom: str, pos: int, records: list[LocusRecord], output: TextIO,
                 duplicates: TextIO, stats: ClassificationStats) -> None:
    counts = Counter(record.allele_key for record in records)
    for record in records:
        if counts[record.allele_key] == 1:
            output.write(record.output_line + "\n")
            stats.output_records += 1
    for (ref, alt), count in counts.items():
        if count > 1:
            stats.distinct_duplicate_keys += 1
            stats.duplicate_key_records += count
            duplicates.write(f"  {chrom}:{pos} {ref}>{alt}\n")


def _write_accounting(path: Path, cohort_id: str, stats: ClassificationStats) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(ACCOUNTING_HEADER) + "\n")
        handle.write(f"{cohort_id}\ttotal_input_records\t{stats.total_input_records}\n")
        for name in VARIANT_CLASSES:
            handle.write(f"{cohort_id}\t{name}_records\t{stats.class_counts[name]}\n")
        handle.write(f"{cohort_id}\tduplicate_key_records\t{stats.duplicate_key_records}\n")
        handle.write(f"{cohort_id}\tdistinct_duplicate_keys\t{stats.distinct_duplicate_keys}\n")
        handle.write(f"{cohort_id}\toutput_records\t{stats.output_records}\n")


def _write_summary(path: Path, cohort_id: str, stats: ClassificationStats,
                   duplicate_path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Variant normalization and classification summary\n")
        handle.write(f"Cohort ID: {cohort_id}\n")
        handle.write("Total input records (post bcftools norm -m- splitting): "
                     f"{stats.total_input_records}\n")
        for name in VARIANT_CLASSES:
            handle.write(f"  {name}: {stats.class_counts[name]}\n")
        handle.write("Only 'snp'-classified records are eligible for the GS panel; "
                     "mnp/indel/symbolic_or_star records are excluded here, not "
                     "silently dropped further downstream.\n")
        handle.write("Duplicate (CHROM, POS, REF, ALT) keys: "
                     f"{stats.distinct_duplicate_keys} distinct key(s), "
                     f"{stats.duplicate_key_records} record(s) total -- every occurrence "
                     "of a colliding key is excluded, not just the extras, since there is "
                     "no automatic way to decide which one is correct.\n")
        if stats.distinct_duplicate_keys:
            handle.write("Colliding keys:\n")
            with duplicate_path.open(encoding="utf-8") as evidence:
                shutil.copyfileobj(evidence, handle)
        handle.write(f"Output records (eligible for GS hard-filtering): {stats.output_records}\n")
        handle.write("FILTER has been reset to '.' on every output record: bcftools norm "
                     "propagates the original record's FILTER value to every split child "
                     "regardless of that child's new shape, which can leave a semantically "
                     "mismatched tag; the GS-specific hard-filter step downstream makes its "
                     "own PASS/FAIL decision from a clean slate.\n")


def stream_classify_vcf(normalized_vcf: Path, cohort_id: str, output_path: Path,
                        accounting_path: Path, summary_path: Path) -> ClassificationStats:
    """Process a contig-grouped, position-sorted VCF with locus-bounded memory."""
    finals = (output_path, accounting_path, summary_path)
    temps = tuple(_temporary_path(path) for path in finals)
    duplicate_path = _temporary_path(summary_path.with_name(summary_path.name + ".duplicates"))
    stats = ClassificationStats(0, {name: 0 for name in VARIANT_CLASSES})
    try:
        with (gzip.open(normalized_vcf, "rt", encoding="utf-8") as source,
              temps[0].open("w", encoding="utf-8") as output,
              duplicate_path.open("w", encoding="utf-8") as duplicates):
            saw_header = False
            current_locus: tuple[str, int] | None = None
            locus_records: list[LocusRecord] = []
            current_contig: str | None = None
            closed_contigs: set[str] = set()
            previous_pos: int | None = None
            for line_number, raw_line in enumerate(source, 1):
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##"):
                    if saw_header:
                        raise MalformedVcfError(f"{normalized_vcf}:{line_number}: metadata after #CHROM header")
                    if not line.startswith("##FILTER="):
                        output.write(line + "\n")
                    continue
                if line.startswith("#CHROM"):
                    if saw_header:
                        raise MalformedVcfError(f"{normalized_vcf}:{line_number}: duplicate #CHROM header")
                    header_fields = line.split("\t")
                    if len(header_fields) <= 9:
                        raise MalformedVcfError(
                            f"{normalized_vcf}: #CHROM header has {len(header_fields)} fields, "
                            "expected at least 10 (9 fixed columns plus one or more samples)"
                        )
                    saw_header = True
                    output.write(line + "\n")
                    continue
                if line.startswith("#"):
                    continue
                if not saw_header:
                    raise MalformedVcfError(f"{normalized_vcf}: data row seen before #CHROM header")

                fields = line.split("\t", 9)
                if len(fields) < 10:
                    raise MalformedVcfError(
                        f"{normalized_vcf}:{line_number}: data row has {len(fields)} "
                        "tab-separated fields, expected at least 10"
                    )
                chrom, pos_text, id_, ref, alt, qual, _filter, info, format_, samples = fields
                try:
                    pos = int(pos_text)
                except ValueError as error:
                    raise MalformedVcfError(
                        f"{normalized_vcf}:{line_number}: POS is not an integer: {pos_text!r}"
                    ) from error
                if pos < 1:
                    raise MalformedVcfError(f"{normalized_vcf}:{line_number}: POS must be positive: {pos}")

                if chrom != current_contig:
                    if chrom in closed_contigs:
                        raise MalformedVcfError(
                            f"{normalized_vcf}:{line_number}: contig {chrom!r} re-entered after it was closed"
                        )
                    if current_contig is not None:
                        closed_contigs.add(current_contig)
                    current_contig, previous_pos = chrom, None
                elif previous_pos is not None and pos < previous_pos:
                    raise MalformedVcfError(
                        f"{normalized_vcf}:{line_number}: POS decreased within contig {chrom!r}: "
                        f"{pos} after {previous_pos}"
                    )
                previous_pos = pos

                locus = chrom, pos
                if current_locus is not None and locus != current_locus:
                    _flush_locus(*current_locus, locus_records, output, duplicates, stats)
                    locus_records.clear()
                current_locus = locus
                stats.total_input_records += 1
                variant_class = classify_variant(ref, alt)
                stats.class_counts[variant_class] += 1
                if variant_class == ELIGIBLE_CLASS:
                    locus_records.append(LocusRecord(
                        ref, alt, "\t".join((chrom, pos_text, id_, ref, alt, qual,
                                             RESET_FILTER_VALUE, info, format_, samples))))
            if not saw_header:
                raise MalformedVcfError(f"{normalized_vcf}: no #CHROM header line found")
            if current_locus is not None:
                _flush_locus(*current_locus, locus_records, output, duplicates, stats)

        _write_accounting(temps[1], cohort_id, stats)
        _write_summary(temps[2], cohort_id, stats, duplicate_path)
        for temporary, final in zip(temps, finals, strict=True):
            os.replace(temporary, final)
        return stats
    finally:
        for temporary in (*temps, duplicate_path):
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Locus-stream normalized variants for GS.")
    parser.add_argument("--normalized-vcf", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--accounting-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        stream_classify_vcf(args.normalized_vcf, args.cohort_id, args.output,
                            args.accounting_output, args.summary_output)
    except (OSError, MalformedVcfError) as error:
        print(f"classify_normalized_variants.py: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
