#!/usr/bin/env python3
"""Record the reference bundle's file identity as a small TSV.

Issue #42: a run manifest that names the reference only by
`reference_id` cannot tell a reader whether two runs actually used the
same bytes. This records a checksum for every reference file the run
actually consumed -- the FASTA, its FAI, its sequence dictionary, and
the five BWA-MEM2 index files -- so "same reference" is checkable
rather than asserted.

The BWA-MEM2 index is included because it is a real mapping input, not
a derived detail: the pipeline accepts a prebuilt index
(`--bwa_index_prefix`) as readily as one it builds itself, and a
mismatched or stale prebuilt index changes alignment results while
leaving `reference_id` and the FASTA checksum identical.

Hashing happens here rather than in the manifest-building process for
the same reason as `hash_input_fastqs.py`: a real reference bundle is
several gigabytes (the index alone dominates), and staging it into the
final task purely to hash it would make the manifest the most expensive
step of the run.

Output contract (headerless, tab-separated):

    role  filename  checksum

`role` is one of `fasta`, `fai`, `dict`, `bwa_index` -- assigned from
the CLI flag the file arrived on, never guessed from its extension, so
a renamed file cannot be silently reclassified. Filenames only, never
paths: this feeds a shareable artifact.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from manifest_utils import sha256_file

REFERENCE_PROVENANCE_COLUMNS: tuple[str, ...] = ("role", "filename", "checksum")

# Roles carrying exactly one file, in the order the manifest presents
# them. `bwa_index` is handled separately: it is a set of five files.
SINGLE_FILE_ROLES: tuple[str, ...] = ("fasta", "fai", "dict")

# The BWA-MEM2 index is exactly these five files, and this is the
# pipeline's existing contract, not a new rule invented here: `main.nf`
# validates the same five suffixes before accepting a prebuilt
# `--bwa_index_prefix`, and `BWA_MEM2_INDEX` produces exactly this set.
#
# The count is checked because a short index is a real, silent failure
# mode: a partially staged or partially copied index bundle would
# otherwise be recorded as if it were the whole mapping input, and the
# manifest would understate what would be needed to reproduce the run.
# The suffixes are checked here, rather than in the manifest reader,
# because this is the side that sees the actual files.
BWA_INDEX_SUFFIXES: tuple[str, ...] = (
    ".0123",
    ".amb",
    ".ann",
    ".bwt.2bit.64",
    ".pac",
)


class MalformedReferenceBundleError(Exception):
    """Raised when the reference bundle cannot be recorded unambiguously."""


def _index_suffix(filename: str) -> str:
    """Return the BWA-MEM2 index suffix a filename carries, or the empty string.

    Matched against the known suffixes rather than parsed off the end of
    the name: `.bwt.2bit.64` is three dot-separated pieces, so anything
    that split on the last dot would classify it as `.64`.
    """
    for suffix in BWA_INDEX_SUFFIXES:
        if filename.endswith(suffix):
            return suffix
    return ""


def build_rows(
    *, fasta: Path, fai: Path, dict_file: Path, bwa_indexes: list[Path]
) -> list[str]:
    """Build one `role<TAB>filename<TAB>checksum` row per reference file."""
    if len(bwa_indexes) != len(BWA_INDEX_SUFFIXES):
        raise MalformedReferenceBundleError(
            f"expected exactly {len(BWA_INDEX_SUFFIXES)} BWA-MEM2 index files "
            f"({', '.join(BWA_INDEX_SUFFIXES)}), got {len(bwa_indexes)}. The run "
            "manifest records the index actually used for mapping; a partial "
            "set would be recorded as though it were the whole mapping input"
        )

    actual_suffixes = sorted(
        _index_suffix(Path(index).name) for index in bwa_indexes
    )
    if actual_suffixes != sorted(BWA_INDEX_SUFFIXES):
        raise MalformedReferenceBundleError(
            "the BWA-MEM2 index files do not match this pipeline's index "
            f"contract: expected suffixes {sorted(BWA_INDEX_SUFFIXES)}, got "
            f"{actual_suffixes}"
        )

    files_by_role: list[tuple[str, Path]] = [
        ("fasta", Path(fasta)),
        ("fai", Path(fai)),
        ("dict", Path(dict_file)),
        *[("bwa_index", Path(index)) for index in bwa_indexes],
    ]

    seen: set[str] = set()
    rows: list[str] = []
    for role, path in files_by_role:
        name = path.name
        if name in seen:
            raise MalformedReferenceBundleError(
                f"duplicate reference file name: '{name}'. Recording only one of "
                "them would silently drop a checksum"
            )
        seen.add(name)
        rows.append(f"{role}\t{name}\tsha256:{sha256_file(path)}")

    return rows


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the reference-provenance CLI."""
    parser = argparse.ArgumentParser(
        description="Record the reference bundle's file identity as a TSV."
    )
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--dict", required=True, type=Path, dest="dict_file")
    parser.add_argument(
        "--bwa-index",
        action="append",
        default=[],
        type=Path,
        dest="bwa_indexes",
        help="A BWA-MEM2 index file actually used for mapping; repeatable.",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        rows = build_rows(
            fasta=args.fasta,
            fai=args.fai,
            dict_file=args.dict_file,
            bwa_indexes=args.bwa_indexes,
        )
    except OSError as error:
        print(f"hash_reference_bundle.py: error: {error}", file=sys.stderr)
        return 1
    except MalformedReferenceBundleError as error:
        print(f"hash_reference_bundle.py: error: {error}", file=sys.stderr)
        return 1

    args.output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
