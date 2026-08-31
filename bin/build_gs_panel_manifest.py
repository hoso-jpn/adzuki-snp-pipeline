#!/usr/bin/env python3
"""Build the reproducibility manifest for a GS panel run.

Mirrors the *shape* of the sibling `genomic-prediction-resnet-hybrid`
repository's own `run_manifest.py` (schema_version, a sortable run_id,
deterministic canonical JSON, filename-only checksums, atomic writes)
without depending on that repository's code -- a bioinformatics
pipeline has no reason to import a model-training repository's Python
package, so the pattern is replicated here in stdlib-only Python
instead.

Software versions are recorded as each GS-lineage process's own
*effective* container identity (schema v2, Issue #52): Nextflow's
`task.container`, captured from inside that exact task after any
`withName`/alias/fully-qualified-selector/profile override has already
been resolved on top of that process's `container` directive default
(itself sourced from `conf/containers.config`, the single default-
value source of truth). This is the ground truth for "what actually
ran" in a Nextflow-plus-Docker pipeline, and cannot silently drift from
reality the way a hand-copied literal digest could -- this script does
not shell out to `bcftools --version` or similar at run time; the
values are passed in via `--container-<process-name>` (one per
recorded process; see `CONTAINER_PROCESS_NAMES` below), wired from
`workflows/adzuki_snp_pipeline.nf` -- except
`--container-build-gs-panel-manifest`, which `BUILD_GS_PANEL_MANIFEST`
reads from its own `task.container` inside its own script, so the
process that writes this manifest records its own container too rather
than leaving the one entry a reader most needs to reproduce the
manifest itself unrecorded. `validate_container_identity()` fails fast
if any value looks like a host filesystem path (bare or `file://`-
wrapped) or embeds registry credentials, since this manifest is a
shareable artifact.

Schema v1 recorded these under a single `containers.bcftools/gatk/
python` triple, one shared value per *tool category* -- ambiguous the
moment two processes using the same tool are overridden differently
(e.g. only `GS_NORMALIZE_VARIANTS`'s bcftools container, not
`GS_INDEX_CLASSIFIED_VARIANTS`'s). Schema v2's `containers` is instead
keyed by GS process name, one entry per process, so two processes
sharing a tool can never be silently collapsed into one value.

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

To let a reader reconstruct exactly which inputs and reference produced
a given panel -- not just verify the panel's own outputs against
themselves -- the checksums also cover the raw/all cohort VCF and the
reference FASTA/FAI used for normalization (`--checksum-file`, passed
by the workflow alongside every `gs_panel/` artifact; there is nothing
input-specific about `checksum_files()` itself, so no separate code
path was needed for this). `--sample-ploidy` and the genotype-encoding
scheme itself (dosage table, missing token, matrix orientation) are
recorded directly in the manifest so a reader never has to cross-
reference this script's source or `docs/gs_panel_data_contract.md` to
know how to interpret the matrix.

`main()` fails fast (exit 1, no manifest written) if `--sample-ploidy`
is not 2, independently of `build_gs_panel.py`'s own identical check:
in the normal pipeline `build_gs_panel.py` always runs first and would
already have failed on the same input, but this script has no way to
know that when invoked on its own, and a manifest whose
`parameters.sample_ploidy` contradicts its own
`genotype_encoding.ploidy == "diploid_only"` would be a self-
contradictory provenance record.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import uuid
from datetime import UTC, datetime
from pathlib import Path

SCHEMA_VERSION = 2

SNP_FILTER_PARAM_NAMES: tuple[str, ...] = (
    "snp_filter_qd_min",
    "snp_filter_qual_min",
    "snp_filter_sor_max",
    "snp_filter_fs_max",
    "snp_filter_mq_min",
    "snp_filter_mq_rank_sum_min",
    "snp_filter_read_pos_rank_sum_min",
)

# Issue #52: every GS-lineage process whose effective container identity
# (Nextflow's task.container, resolved after any override) this manifest
# records -- one manifest `containers` entry per process, never collapsed
# by shared tool ("bcftools"/"gatk"/"python"), so two processes using the
# same tool that are overridden differently can never be silently merged
# into one value. Order here is call order in
# workflows/adzuki_snp_pipeline.nf's GS lineage, not alphabetical.
#
# The list ends with `build_gs_panel_manifest`: this manifest's own
# generating process runs in a container too, and a provenance record
# that silently omitted the one process that wrote it would be a gap in
# exactly the claim this file makes. Its value is passed the same way as
# every other -- a real `task.container` resolved by Nextflow after any
# override -- read inside BUILD_GS_PANEL_MANIFEST's own script (Nextflow
# resolves `task.container` before rendering a task's script, so a
# process can read its own effective container; verified against
# Nextflow 26.04.6 with both a default and a `withName` override).
CONTAINER_PROCESS_NAMES: tuple[str, ...] = (
    "gs_normalize_variants",
    "classify_normalized_variants",
    "gs_index_classified_variants",
    "gatk_variantfiltration_gs",
    "gatk_selectpassvariants_gs",
    "build_gs_panel",
    "reconcile_gs_panel_accounting",
    "build_gs_panel_manifest",
)

# A container identity is expected to be a plain, shareable image
# reference (e.g. `registry/image:tag` or `...@sha256:...`), never a host
# filesystem path (a Singularity `.sif` path would leak host layout) or a
# URL embedding registry credentials -- this manifest is a shareable
# artifact.
#
# A host path reaches this function in two shapes, not one: bare
# (`/opt/images/gatk.sif`, `./gatk.sif`, `~/gatk.sif`) and wrapped in a
# `file://` URI (`file:///opt/images/gatk.sif`), which is what a
# Singularity/Apptainer `container` directive or an image cache can
# legitimately hold. Both disclose host filesystem layout just as fully,
# so both are rejected; a prefix-only check would have let the `file://`
# form through, since it starts with neither `/` nor `.` nor `~`.
_HOST_PATH_PREFIX_RE = re.compile(r"^(/|\./|\.\./|~)")
_FILE_URI_SCHEME_RE = re.compile(r"^file://", re.IGNORECASE)
_URL_CREDENTIALS_RE = re.compile(r"://[^/@\s]+:[^/@\s]+@")


def validate_container_identity(process_name: str, value: str) -> str:
    """Fail fast if a container identity leaks a host path or credential."""
    if not value or not value.strip():
        raise ValueError(f"container identity for '{process_name}' is empty")
    if any(character.isspace() for character in value):
        raise ValueError(
            f"container identity for '{process_name}' contains whitespace: {value!r}"
        )
    if _FILE_URI_SCHEME_RE.match(value):
        raise ValueError(
            f"container identity for '{process_name}' is a file:// URI, which "
            "records a host filesystem path rather than a shareable image "
            f"reference: {value!r}. This manifest is published as a shareable "
            "artifact and must not disclose the host's image layout; pin the "
            "image by registry reference (e.g. 'registry/image@sha256:...') "
            "instead."
        )
    if _HOST_PATH_PREFIX_RE.match(value):
        raise ValueError(
            f"container identity for '{process_name}' looks like a host filesystem "
            f"path rather than an image reference: {value!r}. This manifest is "
            "published as a shareable artifact and must not disclose the host's "
            "image layout; pin the image by registry reference (e.g. "
            "'registry/image@sha256:...') instead."
        )
    if _URL_CREDENTIALS_RE.search(value):
        raise ValueError(
            f"container identity for '{process_name}' appears to embed registry "
            f"credentials: {value!r}"
        )
    return value

# A fixed description of this schema's genotype encoding (see
# bin/build_gs_panel.py and docs/gs_panel_data_contract.md for the
# full reasoning) -- static, not derived from any file, but recorded
# directly in every manifest so a reader never has to cross-reference
# this script's source to know how to interpret the matrix.
GENOTYPE_ENCODING_SCHEMA: dict[str, object] = {
    "schema": "diploid_additive_dosage_v1",
    "dosage_by_genotype": {"0/0": -1, "0/1_or_1/0": 0, "1/1": 1},
    "phasing": "ignored for dosage; 0|1 encodes identically to 0/1",
    "missing_token": "nan",
    "matrix_orientation": "variant_rows_by_sample_columns",
    "ploidy": "diploid_only",
}


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
    containers: dict[str, str],
    sample_ploidy: int,
    snp_filter_params: dict[str, float],
    panel_status: str,
    checksums: dict[str, str],
    run_id: str,
    generated_at: str,
) -> dict[str, object]:
    """Build the manifest document, including its own content hash."""
    parameters: dict[str, object] = {"sample_ploidy": sample_ploidy, **snp_filter_params}

    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "generated_at": generated_at,
        "cohort_id": cohort_id,
        "pipeline_version": pipeline_version,
        "git_commit": git_commit or None,
        "containers": containers,
        "parameters": parameters,
        "genotype_encoding": GENOTYPE_ENCODING_SCHEMA,
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
    for name in CONTAINER_PROCESS_NAMES:
        parser.add_argument(
            f"--container-{name.replace('_', '-')}",
            required=True,
            help=(
                f"Effective container identity ({name.upper()}'s own task.container, "
                "resolved after any override) -- not that process's `container` "
                "directive default."
            ),
        )
    parser.add_argument(
        "--sample-ploidy",
        required=True,
        type=int,
        help="The pipeline's configured sample ploidy (params.sample_ploidy). "
        "This schema is diploid-only: a value other than 2 fails fast (exit 1, "
        "no manifest written) rather than being merely recorded.",
    )

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

    if args.sample_ploidy != 2:
        print(
            "build_gs_panel_manifest.py: error: this GS panel manifest schema "
            f"(v{SCHEMA_VERSION}) is diploid-only, but --sample-ploidy was "
            f"{args.sample_ploidy}. Writing a manifest with "
            "parameters.sample_ploidy set to this value alongside "
            "genotype_encoding.ploidy=\"diploid_only\" would record self-"
            "contradictory provenance. In the normal pipeline this is unreachable "
            "because build_gs_panel.py already fails fast on the same input, but "
            "this script must not depend on that when invoked on its own.",
            file=sys.stderr,
        )
        return 1

    try:
        panel_status = read_panel_status(args.record_accounting)
        checksums = checksum_files(args.checksum_file)
        containers = {
            name: validate_container_identity(
                name, getattr(args, f"container_{name}")
            )
            for name in CONTAINER_PROCESS_NAMES
        }
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
        containers=containers,
        sample_ploidy=args.sample_ploidy,
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
