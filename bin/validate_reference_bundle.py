#!/usr/bin/env python3
"""Bind a reference FASTA, FAI/dictionary and controlled-build BWA index.

The index receipt is written by BWA_MEM2_INDEX in the same task as indexing.
Prebuilt bundles require a previously generated manifest from that trusted
build path. Checksums detect accidental substitution, not a forged receipt.
No claim of reverse-engineering BWA's binary index is made.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

from hash_reference_bundle import BWA_INDEX_SUFFIXES
from manifest_utils import (
    HostMetadataLeakError,
    assert_no_host_metadata,
    sha256_file,
    validate_container_identity,
)

CONTRACT = "reference_bundle_v1"


class BundleError(ValueError):
    """An unverified reference must not reach mapping or variant calling."""


def inspect_fasta(path: Path) -> tuple[list[dict], list[list]]:
    """Read one sequence line at a time; reconstruct complete FAI and M5."""
    contigs, fai = [], []
    seen = set()
    current = None
    with path.open("rb") as handle:
        while True:
            offset = handle.tell()
            raw = handle.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                if current is not None:
                    finish_contig(current, contigs, fai)
                try:
                    name = raw[1:].split()[0].decode("ascii")
                except (IndexError, UnicodeDecodeError) as error:
                    raise BundleError("FASTA has an empty/non-ASCII contig name") from error
                if name in seen:
                    raise BundleError(f"duplicate FASTA contig: {name}")
                seen.add(name)
                current = {
                    "name": name,
                    "length": 0,
                    "offset": handle.tell(),
                    "line_bases": None,
                    "line_width": None,
                    "md5": hashlib.md5(),
                    "last_short": False,
                }
                continue
            if current is None:
                raise BundleError("FASTA sequence appears before its header")
            bases = raw.rstrip(b"\r\n")
            if not bases or not re.fullmatch(b"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", bases):
                raise BundleError(f"unsupported/empty sequence line at byte {offset}")
            if current["last_short"]:
                raise BundleError(
                    "FASTA has a short/nonstandard line before the final sequence line"
                )
            if current["line_bases"] is None:
                # htslib counts an implicit LF at EOF even when the first
                # (and only) sequence line has no terminator. Keep actual
                # CRLF widths and the first-line width of multiline contigs.
                width = len(raw) + (not raw.endswith(b"\n"))
                current["line_bases"], current["line_width"] = len(bases), width
            elif len(bases) > current["line_bases"] or len(raw) > current["line_width"]:
                raise BundleError("FASTA line widths are not indexable")
            current["last_short"] = (
                len(bases) < current["line_bases"] or len(raw) < current["line_width"]
            )
            current["length"] += len(bases)
            current["md5"].update(bases.upper())
    if current is not None:
        finish_contig(current, contigs, fai)
    if not contigs:
        raise BundleError("FASTA contains no contigs")
    return contigs, fai


def finish_contig(current: dict, contigs: list, fai: list) -> None:
    if not current["length"]:
        raise BundleError(f"empty FASTA contig: {current['name']}")
    contigs.append(
        {
            "name": current["name"],
            "length": current["length"],
            "md5": current["md5"].hexdigest(),
        }
    )
    fai.append([current[key] for key in ("name", "length", "offset", "line_bases", "line_width")])


def validate_sidecars(fai_path: Path, dict_path: Path, contigs: list, expected_fai: list) -> None:
    actual_fai = []
    for line in fai_path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 5:
            raise BundleError("FAI must have exactly five fields for uncompressed FASTA")
        try:
            actual_fai.append([fields[0], *map(int, fields[1:])])
        except ValueError as error:
            raise BundleError("FAI has non-integer offsets/lengths") from error
    if actual_fai != expected_fai:
        raise BundleError("FASTA/FAI mismatch (name, length, order, offset or line width)")
    actual_contigs = []
    for line in dict_path.read_text().splitlines():
        if not line.startswith("@SQ\t"):
            continue
        tags = {}
        for field in line.split("\t")[1:]:
            key, sep, value = field.partition(":")
            if not sep or key in tags:
                raise BundleError("dictionary contains malformed/duplicate tags")
            tags[key] = value
        if not {"SN", "LN", "M5"} <= tags.keys():
            raise BundleError("dictionary requires SN, LN and M5; missing M5 is unverified")
        try:
            length = int(tags["LN"])
        except ValueError as error:
            raise BundleError("dictionary LN is not an integer") from error
        actual_contigs.append({"name": tags["SN"], "length": length, "md5": tags["M5"].lower()})
    if actual_contigs != contigs:
        raise BundleError("FASTA/dictionary mismatch (name, length, order or M5)")


def file_record(path: Path) -> dict:
    return {"filename": path.name, "sha256": sha256_file(path)}


def index_records(paths: list[Path], fasta_name: str) -> list[dict]:
    expected = [fasta_name + suffix for suffix in BWA_INDEX_SUFFIXES]
    by_name = {p.name: p for p in paths}
    if len(paths) != len(expected) or set(by_name) != set(expected):
        raise BundleError("BWA index must have the five exact FASTA-prefixed files")
    return [file_record(by_name[name]) for name in expected]


def fingerprint(manifest: dict) -> str:
    content = {key: value for key, value in manifest.items() if key != "fingerprint"}
    return (
        "sha256:"
        + hashlib.sha256(
            json.dumps(content, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
    )


def validate_bundle(
    *,
    fasta: Path,
    fai: Path,
    dictionary: Path,
    indexes: list[Path],
    receipt: Path | None = None,
    tool_version: Path | None = None,
    container: str | None = None,
    prebuilt_manifest: Path | None = None,
) -> dict:
    contigs, expected_fai = inspect_fasta(fasta)
    validate_sidecars(fai, dictionary, contigs, expected_fai)
    fasta_record = file_record(fasta)
    bwa_records = index_records(indexes, fasta.name)
    if (receipt is None) == (prebuilt_manifest is None):
        raise BundleError("require exactly one controlled-build receipt or prebuilt manifest")
    expected_files = {r["filename"]: r["sha256"] for r in [fasta_record, *bwa_records]}
    if receipt is not None:
        seen = {}
        for row in receipt.read_text().splitlines():
            match = re.fullmatch(r"([0-9a-f]{64})  ([^/\\]+)", row)
            if not match or match[2] in seen:
                raise BundleError("invalid/duplicate controlled-build receipt row")
            seen[match[2]] = match[1]
        if seen != expected_files:
            raise BundleError("BWA build receipt does not match FASTA/index bytes")
        version = tool_version.read_text().strip() if tool_version else ""
        if not container or not version:
            raise BundleError("controlled build requires effective container and tool version")
        index_build = {
            "fasta_sha256": fasta_record["sha256"],
            "container": container,
            "tool_version": version,
            "method": "same_task_bwa_mem2_index",
        }
    else:
        try:
            previous = json.loads(prebuilt_manifest.read_text())
        except (ValueError, OSError) as error:
            raise BundleError("prebuilt manifest is missing or invalid JSON") from error
        if not isinstance(previous, dict) or previous.get("contract") != CONTRACT:
            raise BundleError("unsupported prebuilt reference manifest contract")
        index_build = previous.get("index_build")
        if not isinstance(index_build, dict) or set(index_build) != {
            "fasta_sha256",
            "container",
            "tool_version",
            "method",
        }:
            raise BundleError("prebuilt manifest lacks complete controlled-build provenance")
        if (
            index_build["fasta_sha256"] != fasta_record["sha256"]
            or index_build["method"] != "same_task_bwa_mem2_index"
            or not all(isinstance(v, str) and v for v in index_build.values())
        ):
            raise BundleError("prebuilt index is not bound to this FASTA and controlled builder")
    try:
        validate_container_identity("BWA_MEM2_INDEX", index_build["container"])
    except ValueError as error:
        raise BundleError(str(error)) from error
    result = {
        "schema_version": 1,
        "contract": CONTRACT,
        "fasta": fasta_record,
        "fai": file_record(fai),
        "dictionary": file_record(dictionary),
        "contigs": contigs,
        "bwa_index": bwa_records,
        "index_build": index_build,
        "validation": {
            "version": 1,
            "result": "passed",
            "checks": ["fasta_fai", "fasta_dictionary_m5", "index_build_binding"],
            "assurance": "controlled_builder_provenance_not_binary_equivalence",
        },
    }
    try:
        assert_no_host_metadata(result)
    except HostMetadataLeakError as error:
        raise BundleError(str(error)) from error
    result["fingerprint"] = fingerprint(result)
    if prebuilt_manifest is not None and result != previous:
        raise BundleError("prebuilt manifest fingerprint/content differs from supplied bundle")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("fasta", "fai", "dictionary"):
        parser.add_argument("--" + arg, required=True, type=Path)
    parser.add_argument("--bwa-index", action="append", required=True, type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--tool-version", type=Path)
    parser.add_argument("--container")
    parser.add_argument("--prebuilt-manifest", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        result = validate_bundle(
            fasta=args.fasta,
            fai=args.fai,
            dictionary=args.dictionary,
            indexes=args.bwa_index,
            receipt=args.receipt,
            tool_version=args.tool_version,
            container=args.container,
            prebuilt_manifest=args.prebuilt_manifest,
        )
        args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    except (BundleError, OSError) as error:
        print(f"reference bundle rejected: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
