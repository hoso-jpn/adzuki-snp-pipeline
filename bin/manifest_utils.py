#!/usr/bin/env python3
"""Utilities shared by this pipeline's two reproducibility manifests.

`bin/build_run_manifest.py` (whole-run provenance) and
`bin/build_gs_panel_manifest.py` (GS panel provenance) are separate
documents with separate schemas, separate parameter sets, and separate
lifecycles -- but they make the *same* promises about how a manifest is
built: a sortable `run_id`, UTC timestamps, SHA-256 checksums keyed by
filename only, deterministic canonical JSON, a self-referential
`manifest_hash`, an atomic write, and a fail-fast refusal to record a
container identity that would leak a host filesystem path or a registry
credential.

Those promises were previously implemented twice, character for
character (Issue #42's own investigation compared the two files' ASTs
and found seven functions identical down to their docstrings). Two
copies of a hashing/serialization contract is exactly the kind of
duplication that drifts silently: a fix to one manifest's canonical
JSON would change its `manifest_hash` while the other kept the old
behavior, and nothing would fail. This module is the single
implementation both import.

Deliberately *not* shared here: either manifest's payload shape, schema
version, CLI, or parameter set. Those differ on purpose -- the run
manifest records indel filter thresholds and per-sample accounting the
GS manifest has no use for, and generalizing them into one "manifest
builder" would force each document to carry the other's fields. Only
the mechanics are common; the contracts are not.

Import contract: this module sits next to its callers in `bin/`, and is
imported as a plain top-level module (`from manifest_utils import ...`)
with no package, no installation, and no PYTHONPATH setup. That works
in both places these scripts run, because in both of them Python puts
the calling script's own directory first on `sys.path`:

  * from a source checkout -- `python3 bin/build_run_manifest.py ...`;
  * inside a Nextflow task container -- Nextflow stages the project's
    whole `bin/` directory and puts it on PATH, so the script that runs
    is `<staged bin>/build_run_manifest.py` and `sys.path[0]` is that
    same staged directory (verified on Nextflow 26.04.6 with Docker
    against this repository's own pinned Python image).

Consequently this module must stay in `bin/` alongside its callers, and
must remain stdlib-only: a task container is this repository's pinned
`python:3.12` image with no `pip install` step.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
from datetime import UTC, datetime
from pathlib import Path

# A container identity is expected to be a plain, shareable image
# reference (e.g. `registry/image:tag` or `...@sha256:...`), never a host
# filesystem path (a Singularity `.sif` path would leak host layout) or a
# URL embedding registry credentials -- a manifest is a shareable
# artifact.
#
# A host path reaches the validator in two shapes, not one: bare
# (`/opt/images/gatk.sif`, `./gatk.sif`, `~/gatk.sif`) and wrapped in a
# `file://` URI (`file:///opt/images/gatk.sif`), which is what a
# Singularity/Apptainer `container` directive or an image cache can
# legitimately hold. Both disclose host filesystem layout just as fully,
# so both are rejected; a prefix-only check would have let the `file://`
# form through, since it starts with neither `/` nor `.` nor `~`.
_HOST_PATH_PREFIX_RE = re.compile(r"^(/|\./|\.\./|~)")
_FILE_URI_SCHEME_RE = re.compile(r"^file://", re.IGNORECASE)
_URL_CREDENTIALS_RE = re.compile(r"://[^/@\s]+:[^/@\s]+@")


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


def write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    """Write JSON via a same-directory temp file, then rename it into place.

    A failure partway through leaves the original file (if any)
    untouched and no partially-written file at the final path.

    Serialization happens before the temp file is created, so a payload
    that cannot be serialized raises without touching the filesystem at
    all. A failure of the write or the rename itself removes the temp
    file rather than leaving a stale `.<name>.tmp` next to the final
    path: a manifest directory is a published artifact directory, and a
    leftover partial file there is indistinguishable, to a later reader,
    from something a run meant to produce.
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
    try:
        tmp_path.write_text(text + "\n", encoding="utf-8")
        os.replace(tmp_path, path)
    except OSError:
        tmp_path.unlink(missing_ok=True)
        raise
