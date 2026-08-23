#!/usr/bin/env python3
"""Validate that a FASTA .fai index and a GATK sequence dictionary (.dict)
describe the exact same, identically ordered, identically sized contig list.

`GATK_HAPLOTYPECALLER`, `GATK_GENOMICSDBIMPORT`, and `GATK_GENOTYPEGVCFS` all
take the FASTA, its `.fai`, and its `.dict` as three separate inputs (Issue
#4's reference-bundle contract allows each to be either pipeline-generated or
independently supplied via `--reference_fai`/`--reference_dict`). Nothing
upstream of those GATK invocations currently confirms the `.fai` and `.dict`
actually describe the same reference in the same contig order -- a caller who
supplies a mismatched pair (e.g. a `.dict` regenerated from a different
FASTA, or a `.fai`/`.dict` with contigs listed in a different order) would
not find out until whichever GATK process happens to notice first, with
whatever error message that specific tool chooses to print, at whatever point
in a run that may be.

This tool is deliberately narrow: it does not itself decide whether a
mismatch is fatal (that is the caller's job -- see
`modules/local/validate_reference_contigs.nf`, which runs this once, before
any GATK process starts, for both the pipeline-generated and the
user-supplied prebuilt reference-bundle path). It only answers one question
--exactly, and with an actionable message identifying the first point of
disagreement -- whether the two contig lists are identical, in the same
order, matched by name, and by declared length.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ContigRecord:
    """One contig's name and declared length, in file-appearance order."""

    name: str
    length: int


class MalformedFaiError(Exception):
    """A `.fai` line could not be parsed as `NAME\tLENGTH\t...`."""


class MalformedDictError(Exception):
    """A `.dict` `@SQ` line was missing a required `SN:`/`LN:` field."""


class DuplicateContigError(Exception):
    """The same contig name appears more than once in a single file."""


def parse_fai(path: Path) -> list[ContigRecord]:
    """Parse a `samtools faidx`-format `.fai` index into ordered contig records.

    Each non-blank line is `NAME\tLENGTH\tOFFSET\tLINEBASES\tLINEWIDTH`; only
    the first two fields are used here.
    """
    records: list[ContigRecord] = []
    seen: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise MalformedFaiError(
                    f"{path}: line {line_number} has {len(fields)} field(s), "
                    "expected at least NAME and LENGTH"
                )
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as error:
                raise MalformedFaiError(
                    f"{path}: line {line_number}: length {fields[1]!r} is not an integer"
                ) from error
            seen[name] += 1
            records.append(ContigRecord(name=name, length=length))

    duplicates = sorted(name for name, count in seen.items() if count > 1)
    if duplicates:
        raise DuplicateContigError(f"{path}: duplicate contig name(s): {duplicates}")

    return records


def parse_dict(path: Path) -> list[ContigRecord]:
    """Parse a GATK/Picard sequence dictionary's `@SQ` lines into ordered contig records.

    Each `@SQ` line carries tab-separated `KEY:VALUE` fields in no fixed
    order beyond `@SQ` itself being first; only `SN:` (name) and `LN:`
    (length) are required here.
    """
    records: list[ContigRecord] = []
    seen: Counter[str] = Counter()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line.startswith("@SQ"):
                continue
            name: str | None = None
            length: int | None = None
            for field in line.split("\t")[1:]:
                if field.startswith("SN:"):
                    name = field[len("SN:") :]
                elif field.startswith("LN:"):
                    length_field = field[len("LN:") :]
                    try:
                        length = int(length_field)
                    except ValueError as error:
                        raise MalformedDictError(
                            f"{path}: line {line_number}: LN value {length_field!r} "
                            "is not an integer"
                        ) from error
            if name is None or length is None:
                raise MalformedDictError(
                    f"{path}: line {line_number}: @SQ record is missing SN: and/or LN:"
                )
            seen[name] += 1
            records.append(ContigRecord(name=name, length=length))

    duplicates = sorted(name for name, count in seen.items() if count > 1)
    if duplicates:
        raise DuplicateContigError(f"{path}: duplicate contig name(s): {duplicates}")

    return records


def find_first_mismatch(
    fai_records: list[ContigRecord], dict_records: list[ContigRecord]
) -> str | None:
    """Return a human-readable description of the first disagreement, or `None`.

    Compares position-by-position (not as sets), so two files that list the
    exact same contig names and lengths but in a different order are
    correctly reported as a mismatch, not silently accepted.
    """
    if len(fai_records) != len(dict_records):
        return (
            f"contig count differs: FAI has {len(fai_records)} contig(s), "
            f"sequence dictionary has {len(dict_records)} contig(s)"
        )

    for index, (fai_record, dict_record) in enumerate(zip(fai_records, dict_records, strict=True)):
        if fai_record.name != dict_record.name or fai_record.length != dict_record.length:
            return (
                f"first mismatch at index {index}: "
                f"FAI: {fai_record.name} length={fai_record.length}, "
                f"DICT: {dict_record.name} length={dict_record.length}"
            )

    return None


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the reference-contig validation CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Validate that a FASTA .fai index and a GATK sequence dictionary "
            "(.dict) describe the same contigs, in the same order, with the "
            "same declared lengths."
        )
    )
    parser.add_argument("--fai", required=True, type=Path, help="Path to the FASTA .fai index.")
    parser.add_argument(
        "--dict", required=True, type=Path, dest="dict_path", help="Path to the .dict file."
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        fai_records = parse_fai(args.fai)
        dict_records = parse_dict(args.dict_path)
    except OSError as error:
        print(f"validate_reference_contigs.py: error: {error}", file=sys.stderr)
        return 1
    except (MalformedFaiError, MalformedDictError, DuplicateContigError) as error:
        print(f"validate_reference_contigs.py: error: {error}", file=sys.stderr)
        return 1

    if not fai_records:
        print(
            f"validate_reference_contigs.py: error: {args.fai}: contains no contigs",
            file=sys.stderr,
        )
        return 1
    if not dict_records:
        print(
            f"validate_reference_contigs.py: error: {args.dict_path}: contains no @SQ records",
            file=sys.stderr,
        )
        return 1

    mismatch = find_first_mismatch(fai_records, dict_records)
    if mismatch is not None:
        print(
            "validate_reference_contigs.py: error: reference FAI and sequence "
            f"dictionary are inconsistent: {mismatch}",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
