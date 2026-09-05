#!/usr/bin/env python3
"""Record one read group's input FASTQ provenance as a single TSV row.

Issue #42: the run-level manifest has to name every raw input file the
run consumed, with a checksum, so a reader can tell whether two runs saw
the same data. The obvious way to do that -- hand every FASTQ in the
cohort to the manifest-building process and hash them all there -- does
not survive contact with a real cohort: this pipeline is meant to scale
to 327 samples, and that design would re-stage every raw FASTQ in the
run into one task and hash them serially, purely to write a JSON file.

So the hashing happens here instead, once per read group, in tasks that
fan out exactly like the rest of the pipeline and that Nextflow can
cache and resume individually. Each task writes one small row; the
manifest builder reads the rows, not the reads.

Output contract (headerless, tab-separated, exactly one row):

    rank  sample_id  read_group_id  library_id  platform  platform_unit
          fastq_1_filename  fastq_1_checksum  fastq_2_filename
          fastq_2_checksum

`rank` is the row's zero-based position in the samplesheet, zero-padded
to 8 digits. It exists so the manifest can restore samplesheet order
after Nextflow's parallel, order-free execution -- Nextflow gives no
ordering guarantee across tasks, and "the order the samplesheet listed
them in" is the only ordering a reader of the manifest can check
against their own input. It is a transport detail: `build_run_manifest.py`
sorts on it and then drops it, so it never appears in the manifest.

Headerless because these rows are concatenated across tasks; a header
per task would interleave with the data. The row is deliberately
narrow: filenames, never paths (a path would record the host's
directory layout in a shareable artifact), and checksums, never
sequence content.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from manifest_utils import sha256_file

# The transport format between this process and build_run_manifest.py.
# Both sides name the columns from this one tuple so a change here
# cannot desynchronize the writer from the reader.
INPUT_PROVENANCE_COLUMNS: tuple[str, ...] = (
    "rank",
    "sample_id",
    "read_group_id",
    "library_id",
    "platform",
    "platform_unit",
    "fastq_1_filename",
    "fastq_1_checksum",
    "fastq_2_filename",
    "fastq_2_checksum",
)

RANK_WIDTH = 8


class MalformedProvenanceFieldError(Exception):
    """Raised when a metadata value cannot be represented in the TSV row."""


def validate_field(name: str, value: str) -> str:
    """Reject any value that would corrupt the row's own tab/newline framing."""
    if "\t" in value or "\n" in value or "\r" in value:
        raise MalformedProvenanceFieldError(
            f"{name} contains a tab or newline, which the TSV row format cannot "
            f"represent: {value!r}"
        )
    return value


def build_row(
    *,
    rank: int,
    sample_id: str,
    read_group_id: str,
    library_id: str,
    platform: str,
    platform_unit: str,
    fastq_1: Path,
    fastq_2: Path,
) -> str:
    """Build the single TSV row for one read group, hashing both FASTQs."""
    if rank < 0:
        raise MalformedProvenanceFieldError(f"rank must not be negative: {rank}")

    fields = [
        str(rank).zfill(RANK_WIDTH),
        validate_field("sample_id", sample_id),
        validate_field("read_group_id", read_group_id),
        validate_field("library_id", library_id),
        validate_field("platform", platform),
        validate_field("platform_unit", platform_unit),
        validate_field("fastq_1", Path(fastq_1).name),
        f"sha256:{sha256_file(Path(fastq_1))}",
        validate_field("fastq_2", Path(fastq_2).name),
        f"sha256:{sha256_file(Path(fastq_2))}",
    ]
    assert len(fields) == len(INPUT_PROVENANCE_COLUMNS)
    return "\t".join(fields)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the input-provenance CLI."""
    parser = argparse.ArgumentParser(
        description="Record one read group's input FASTQ provenance as a TSV row."
    )
    parser.add_argument(
        "--rank",
        required=True,
        type=int,
        help="Zero-based samplesheet row position, used only to restore "
        "samplesheet order in the run manifest.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-group-id", required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument(
        "--platform-unit",
        required=True,
        help="May be empty: the samplesheet contract allows an empty platform_unit.",
    )
    parser.add_argument("--fastq-1", required=True, type=Path)
    parser.add_argument("--fastq-2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        row = build_row(
            rank=args.rank,
            sample_id=args.sample_id,
            read_group_id=args.read_group_id,
            library_id=args.library_id,
            platform=args.platform,
            platform_unit=args.platform_unit,
            fastq_1=args.fastq_1,
            fastq_2=args.fastq_2,
        )
    except OSError as error:
        print(f"hash_input_fastqs.py: error: {error}", file=sys.stderr)
        return 1
    except MalformedProvenanceFieldError as error:
        print(f"hash_input_fastqs.py: error: {error}", file=sys.stderr)
        return 1

    args.output.write_text(row + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
