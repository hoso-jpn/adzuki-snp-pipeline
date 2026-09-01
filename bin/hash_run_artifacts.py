#!/usr/bin/env python3
"""Checksum one group of produced run artifacts into a small TSV.

Issue #42: the run manifest's `checksums` field covers the run's own
scientific deliverables (per-sample gVCFs, the cohort raw VCF, the
primary SNP/indel PASS VCFs). Hashing them all inside the
manifest-building process would put every gVCF in the cohort -- hundreds
of gigabytes at this pipeline's target scale -- through one task purely
to produce a JSON file.

This process is invoked once per artifact group instead (per sample for
gVCFs, once for the cohort-level VCFs), so hashing fans out with the
rest of the pipeline and each result is cached and resumed
independently. It hashes the artifact a process actually produced, taken
straight from that process's own output channel, rather than looking up
a published copy by path -- the published file is a copy of this exact
file, and reading the channel keeps the checksum tied to the task that
made it.

Output contract (headerless, tab-separated):

    filename  checksum

Filenames only, never paths. `build_run_manifest.py` merges the rows
from every group and fails on a duplicate filename rather than silently
keeping one of two different files that happen to share a basename.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from manifest_utils import sha256_file

ARTIFACT_CHECKSUM_COLUMNS: tuple[str, ...] = ("filename", "checksum")


class MalformedArtifactGroupError(Exception):
    """Raised when an artifact group cannot be recorded unambiguously."""


def build_rows(artifacts: list[Path]) -> list[str]:
    """Build one `filename<TAB>checksum` row per artifact in this group."""
    if not artifacts:
        raise MalformedArtifactGroupError(
            "no artifacts were given; an empty checksum group would silently "
            "record nothing for a stage the manifest claims to cover"
        )

    seen: set[str] = set()
    rows: list[str] = []
    for artifact in artifacts:
        path = Path(artifact)
        if path.name in seen:
            raise MalformedArtifactGroupError(
                f"duplicate artifact file name: '{path.name}'"
            )
        seen.add(path.name)
        rows.append(f"{path.name}\tsha256:{sha256_file(path)}")

    return rows


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the artifact-checksum CLI."""
    parser = argparse.ArgumentParser(
        description="Checksum one group of produced run artifacts into a TSV."
    )
    parser.add_argument(
        "--artifact",
        action="append",
        default=[],
        type=Path,
        dest="artifacts",
        help="A produced run artifact to checksum; repeatable.",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        rows = build_rows(args.artifacts)
    except OSError as error:
        print(f"hash_run_artifacts.py: error: {error}", file=sys.stderr)
        return 1
    except MalformedArtifactGroupError as error:
        print(f"hash_run_artifacts.py: error: {error}", file=sys.stderr)
        return 1

    args.output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
