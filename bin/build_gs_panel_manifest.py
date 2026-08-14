#!/usr/bin/env python3
"""Build the reproducibility manifest for a GS panel run.

Mirrors the *shape* of the sibling `genomic-prediction-resnet-hybrid`
repository's own `run_manifest.py` (schema_version, a sortable run_id,
deterministic canonical JSON, filename-only checksums, atomic writes)
without depending on that repository's code -- a bioinformatics
pipeline has no reason to import a model-training repository's Python
package, so the pattern is replicated here in stdlib-only Python
instead.

Software versions are recorded as the pinned container image
references already baked into each process's `container` directive
(`modules/local/gs_normalize_variants.nf`, `gatk_variantfiltration.nf`,
`build_gs_panel.nf`, ...), which are the ground truth for "what
actually ran" in a Nextflow-plus-Docker pipeline; this script does not
shell out to `bcftools --version` or similar at run time; passed in
via `--bcftools-container`/`--gatk-container`/`--python-container`
instead.

The git commit is passed in via `--git-commit`, resolved by the
workflow from Nextflow's own `workflow.commitId` rather than by
running `git rev-parse` inside a task container: the task work
directory holds only staged input files, never the pipeline's `.git`
checkout, so a container-side `git rev-parse` would find no repository
at all, not merely an excluded one. `workflow.commitId` is documented
as populated only when Nextflow itself pulls a git-hosted pipeline
(`nextflow run owner/repo`), so it can legitimately be empty for this
repository's own documented `nextflow run .` local-directory
invocation; an empty value is recorded as `null`, honestly, rather
than worked around.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import uuid
from datetime import UTC, datetime
from pathlib import Path

SCHEMA_VERSION = 1

SNP_FILTER_PARAM_NAMES: tuple[str, ...] = (
    "snp_filter_qd_min",
    "snp_filter_qual_min",
    "snp_filter_sor_max",
    "snp_filter_fs_max",
    "snp_filter_mq_min",
    "snp_filter_mq_rank_sum_min",
    "snp_filter_read_pos_rank_sum_min",
)


class MalformedAccountingError(Exception):
    """Raised when the record-accounting TSV is missing required data."""


def new_run_id(now: datetime | None = None, suffix: str | None = None) -> str:
    """Build a sortable, human-readable run identifier."""
    moment = now if now is not None else datetime.now(UTC)
    timestamp = moment.strftime("%Y%m%dT%H%M%SZ")
    token = suffix if suffix is not None else uuid.uuid4().hex[:8]
    return f"{timestamp}-{token}"


def utc_now_iso(now: datetime | None = None) -> str:
    """Format a UTC timestamp as an ISO-8601 string with a 'Z' suffix."""
    moment = now if now is not None else datetime.now(UTC)
    return moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path, chunk_size: int = 1_048_576) -> str:
    """Compute a file's SHA-256 hex digest without loading it fully into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checksum_files(paths: list[Path]) -> dict[str, str]:
    """Checksum each path, keyed by filename only (never an absolute path).

    Raises if two paths share a filename: silently keeping only the
    last one would drop a checksum without any indication it happened.
    """
    checksums: dict[str, str] = {}
    for path in paths:
        name = Path(path).name
        if name in checksums:
            raise ValueError(f"duplicate checksum file name: '{name}'")
        checksums[name] = f"sha256:{sha256_file(Path(path))}"
    return checksums


def read_panel_status(path: Path) -> str:
    """Read the panel_status metric from reconcile_gs_panel_accounting.py's TSV."""
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            raise MalformedAccountingError(
                f"{path}: row has {len(fields)} tab-separated fields, expected at least 3"
            )

        if fields[1] == "panel_status":
            return fields[2]

    raise MalformedAccountingError(f"{path}: missing required metric: panel_status")


def _json_default(value: object) -> object:
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def canonical_json_hash(payload: dict[str, object]) -> str:
    """Hash a JSON-serializable document deterministically, order-sensitive."""
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
        default=_json_default,
    )
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def build_manifest(
    *,
    cohort_id: str,
    pipeline_version: str,
    git_commit: str,
    bcftools_container: str,
    gatk_container: str,
    python_container: str,
    snp_filter_params: dict[str, float],
    panel_status: str,
    checksums: dict[str, str],
    run_id: str,
    generated_at: str,
) -> dict[str, object]:
    """Build the manifest document, including its own content hash."""
    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "generated_at": generated_at,
        "cohort_id": cohort_id,
        "pipeline_version": pipeline_version,
        "git_commit": git_commit or None,
        "containers": {
            "bcftools": bcftools_container,
            "gatk": gatk_container,
            "python": python_container,
        },
        "parameters": snp_filter_params,
        "panel_status": panel_status,
        "checksums": checksums,
    }
    manifest["manifest_hash"] = canonical_json_hash(manifest)
    return manifest


def write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    """Write JSON via a same-directory temp file, then rename it into place.

    A failure partway through leaves the original file (if any)
    untouched and no partially-written file at the final path.
    """
    text = json.dumps(
        payload,
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
        default=_json_default,
    )
    tmp_path = path.with_name(f".{path.name}.tmp")
    tmp_path.write_text(text + "\n", encoding="utf-8")
    os.replace(tmp_path, path)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the manifest-builder CLI."""
    parser = argparse.ArgumentParser(
        description="Build the reproducibility manifest for a GS panel run."
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--pipeline-version", required=True, help="Pipeline manifest version (workflow.manifest.version)."
    )
    parser.add_argument(
        "--git-commit",
        default="",
        help="Nextflow workflow.commitId; empty string is recorded as null.",
    )
    parser.add_argument("--bcftools-container", required=True, help="Pinned bcftools image reference.")
    parser.add_argument("--gatk-container", required=True, help="Pinned GATK image reference.")
    parser.add_argument("--python-container", required=True, help="Pinned Python image reference.")

    for name in SNP_FILTER_PARAM_NAMES:
        parser.add_argument(f"--{name.replace('_', '-')}", required=True, type=float)

    parser.add_argument(
        "--record-accounting",
        required=True,
        type=Path,
        help="Path to reconcile_gs_panel_accounting.py's output TSV (read for panel_status).",
    )
    parser.add_argument(
        "--checksum-file",
        action="append",
        default=[],
        type=Path,
        help="A GS panel artifact to checksum; repeatable.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Output path for the manifest JSON."
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    try:
        panel_status = read_panel_status(args.record_accounting)
        checksums = checksum_files(args.checksum_file)
    except OSError as error:
        print(f"build_gs_panel_manifest.py: error: {error}", file=sys.stderr)
        return 1
    except (MalformedAccountingError, ValueError) as error:
        print(f"build_gs_panel_manifest.py: error: {error}", file=sys.stderr)
        return 1

    snp_filter_params = {name: getattr(args, name) for name in SNP_FILTER_PARAM_NAMES}

    manifest = build_manifest(
        cohort_id=args.cohort_id,
        pipeline_version=args.pipeline_version,
        git_commit=args.git_commit,
        bcftools_container=args.bcftools_container,
        gatk_container=args.gatk_container,
        python_container=args.python_container,
        snp_filter_params=snp_filter_params,
        panel_status=panel_status,
        checksums=checksums,
        run_id=new_run_id(),
        generated_at=utc_now_iso(),
    )

    write_json_atomic(args.output, manifest)

    return 0


if __name__ == "__main__":
    sys.exit(main())
