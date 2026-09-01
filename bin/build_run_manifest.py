#!/usr/bin/env python3
"""Build a whole-run reproducibility manifest for a completed pipeline execution.

Issue #26: a commissioned delivery needs more than the GS panel's own
manifest (`build_gs_panel_manifest.py`, scoped to `gs_panel/` artifacts
alone) -- it needs one record covering the whole run: which samples,
which reference, which container/tool versions, which Nextflow
parameters, and what the cohort/variant accounting came out to. This
script builds exactly that, following the same shape as
`build_gs_panel_manifest.py` (schema_version, a sortable run_id,
deterministic canonical JSON, filename-only checksums, atomic writes).
Those shared mechanics live in `bin/manifest_utils.py` (Issue #42) --
one implementation both manifests import, rather than two copies that
could drift into disagreeing about what a `manifest_hash` means. The
two manifests' *payloads* remain separate: neither document's schema,
parameter set, or CLI is generalized into the other's.

This script has two modes, and they build two different schema
versions on purpose.

`--mode legacy-v1` (the default, and the only behavior before Issue
#42) is the original standalone CLI: run once, by hand, against a
completed run's own already-published artifacts (samplesheet, reference
bundle, QC/accounting TSVs, and whichever files the caller wants
checksummed via repeatable `--checksum-file`). It records software
versions as four tool-level container references passed in via
`--bwa-mem2-container`/`--samtools-container`/`--gatk-container`/
`--python-container`. This mode is unchanged: the Issue #26/#33 real
cohort manifests in `docs/` were produced under it, they are not being
regenerated, and a v1 document must go on meaning exactly what it meant
when it was written.

`--mode dag-v2` (Issue #42) is what the Nextflow workflow itself runs,
as an ordinary required process at the end of the DAG, so a successful
run always produces its own provenance rather than depending on someone
remembering to run this afterwards. It differs from v1 in three ways
that could not be expressed as additions to v1:

  * `containers` is keyed by *process*, not by tool. Every containerized
    process emits its own effective `task.container` (the mechanism
    Issue #52 established for the GS panel manifest), so FastQC, fastp,
    MultiQC, bcftools and the two independently-overridable aliases of
    the shared GATK modules are each recorded as themselves. v1's four
    tool-level keys cannot represent that, and reusing the same field
    name for a different meaning would silently change what every
    already-published v1 manifest appears to claim.
  * `git_commit` is required to be a full 40-character SHA rather than
    nullable. A provenance artifact the pipeline produces for itself has
    no excuse for not knowing which code produced it.
  * `sample_accounting` records the per-sample QC the pipeline already
    computes (`summarize_variant_qc.py`'s raw/all `sample_qc.tsv`),
    cross-checked against the cohort-level accounting rather than
    recomputed.

In DAG mode this script hashes nothing large itself. Input FASTQs,
the reference bundle and the run's output artifacts are checksummed by
their own small processes (`hash_input_fastqs.py`,
`hash_reference_bundle.py`, `hash_run_artifacts.py`), which fan out
across the cluster the way the rest of the pipeline does; this script
reads their tiny TSVs. Staging a whole cohort's raw reads and gVCFs
into one terminal task purely to write a JSON file does not scale to
the 327-sample cohort this pipeline is aimed at.

Neither mode records anything host-specific: files are recorded by
basename and checksum, and no working directory, launch directory,
command line or user name is passed in at all (see
`manifest_utils.assert_no_host_metadata`, the backstop for that).

The Nextflow execution engine's own version (`--nextflow-version`,
e.g. Nextflow's `nextflow.version`/`-v`) is recorded separately from
those containers: Nextflow itself runs on the host, outside any pinned
container, and its own runtime semantics (channel/path handling,
scheduling, retry behavior) have previously been the actual root cause
of a real production failure in this pipeline (Issue #11's
single-element List/scalar collapse). Omitting it from a
"reproducibility" manifest would leave out exactly the kind of
execution-engine detail that mattered before.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# Issue #42: the hashing/serialization mechanics this manifest and
# bin/build_gs_panel_manifest.py both depend on live in one module now
# (they were previously duplicated character for character). See
# bin/manifest_utils.py for why a plain sibling import needs no
# packaging or PYTHONPATH in either place these scripts run.
from manifest_utils import (
    HostMetadataLeakError,
    _json_default,
    assert_no_host_metadata,
    canonical_json_hash,
    checksum_files,
    new_run_id,
    sha256_file,
    utc_now_iso,
    validate_container_identity,
    write_json_atomic,
)

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

INDEL_FILTER_PARAM_NAMES: tuple[str, ...] = (
    "indel_filter_qd_min",
    "indel_filter_qual_min",
    "indel_filter_fs_max",
    "indel_filter_read_pos_rank_sum_min",
)

REQUIRED_SAMPLESHEET_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "read_group_id",
    "fastq_1",
    "fastq_2",
    "library_id",
    "platform",
    "platform_unit",
)

# Issue #42: the schema this script writes when the Nextflow DAG builds
# the manifest itself. It is a version bump rather than an extension of
# v1 because `containers` changes meaning: v1 recorded four *tool*
# identities (bwa_mem2/samtools/gatk/python) hand-copied from four
# `container` directives, which cannot express a pipeline where FastQC,
# fastp, MultiQC and bcftools each run in their own image and any
# process can be overridden individually. v2 records one entry per
# process actually executed, read from that process's own effective
# `task.container` (the mechanism Issue #52 established for the GS panel
# manifest). Silently redefining a v1 field's meaning would leave every
# already-published v1 manifest looking like it made a claim it never
# made.
DAG_SCHEMA_VERSION = 2

LEGACY_MODE = "legacy-v1"
DAG_MODE = "dag-v2"

# The 40-character hexadecimal commit the run was executed from. A DAG
# manifest requires it: an abbreviated or absent SHA cannot identify the
# code that produced a delivered result. `workflow.commitId` is null for
# this repository's own documented `nextflow run .` invocation (measured
# on Nextflow 26.04.6), so the workflow resolves it explicitly and this
# check is what stops that resolution from silently degrading.
FULL_GIT_SHA_RE = re.compile(r"^[0-9a-f]{40}$")

# A canonical process key: the lowercase process/alias name the workflow
# assigns each invocation. Constrained so a malformed provenance row
# cannot introduce an arbitrary key into a published document.
PROCESS_KEY_RE = re.compile(r"^[a-z][a-z0-9_]*$")

# Column contracts of the staged TSVs the DAG hands this script. Each is
# written by the correspondingly named bin/ script; the counts are
# asserted on read so a truncated or reshaped row fails loudly instead of
# being silently misread by position.
INPUT_PROVENANCE_FIELD_COUNT = 10
REFERENCE_PROVENANCE_FIELD_COUNT = 3
ARTIFACT_CHECKSUM_FIELD_COUNT = 2
RUNTIME_PROVENANCE_FIELD_COUNT = 2

REFERENCE_SINGLE_FILE_ROLES: tuple[str, ...] = ("fasta", "fai", "dict")
REFERENCE_MULTI_FILE_ROLES: tuple[str, ...] = ("bwa_index",)

# How many files each multi-file reference role must carry. The BWA-MEM2
# index is exactly five files -- the same contract `main.nf` enforces on a
# prebuilt `--bwa_index_prefix` and `hash_reference_bundle.py` enforces on
# the files themselves. Checked again here because this side reads a
# staged TSV: a truncated or partially written provenance file would
# otherwise produce a manifest that silently understates the mapping input
# it claims to describe. (The suffix names are validated by the writer,
# which sees the real files; this side validates the shape of what
# reached it.)
REFERENCE_MULTI_FILE_ROLE_COUNTS: dict[str, int] = {"bwa_index": 5}

# summarize_variant_qc.py's sample_qc.tsv header, in order. The run
# manifest re-uses that existing artifact as its per-sample accounting
# rather than recomputing per-sample statistics: a second implementation
# of the same numbers is a second chance to disagree with the QC the
# pipeline already published.
SAMPLE_QC_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "sample",
    "reference_homozygous",
    "non_reference_homozygous",
    "heterozygous",
    "missing",
    "missingness_rate",
    "average_depth",
    "singletons",
)

SAMPLE_ACCOUNTING_FIELDS: tuple[str, ...] = SAMPLE_QC_HEADER[3:]

VARIANT_QC_HEADER: tuple[str, ...] = (
    "cohort_id",
    "stage",
    "variant_type",
    "metric",
    "value",
)

# The stage/variant-type slice of the QC outputs the run manifest reads.
# Both the cohort and per-sample accounting come from the raw/all run --
# the cohort VCF before any filtering -- so the two describe the same
# variants and can be cross-checked against each other.
COHORT_ACCOUNTING_STAGE = "raw"
COHORT_ACCOUNTING_VARIANT_TYPE = "all"

REQUIRED_COHORT_ACCOUNTING_METRICS: tuple[str, ...] = (
    "number_of_samples",
    "sample_names",
    "cohort_total_genotypes",
    "cohort_missing_genotypes",
    "cohort_missingness_rate",
)

REQUIRED_VARIANT_TYPE_ACCOUNTING_METRICS: tuple[str, ...] = (
    "raw_all_records",
    "raw_snp_records",
    "raw_indel_records",
)

REQUIRED_GS_PANEL_MANIFEST_KEYS: tuple[str, ...] = (
    "schema_version",
    "run_id",
    "cohort_id",
    "panel_status",
    "manifest_hash",
)


class MalformedSamplesheetError(Exception):
    """Raised when the samplesheet is missing a required column or is empty."""


class MalformedAccountingError(Exception):
    """Raised when a metric TSV is missing required data or malformed."""


def parse_samplesheet(path: Path) -> list[dict[str, str]]:
    """Parse the pipeline's own samplesheet contract into a list of row dicts.

    Every row's `fastq_1`/`fastq_2` are also checksummed here (not left
    for the caller to pass separately via `--checksum-file`): the
    samplesheet is the one place that already names every raw input
    file this run consumed, so re-deriving that list independently
    would risk drifting from what was actually run.
    """
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise MalformedSamplesheetError(f"{path}: file is empty")

        missing = [column for column in REQUIRED_SAMPLESHEET_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise MalformedSamplesheetError(f"{path}: missing required column(s): {missing}")

        rows = list(reader)

    if not rows:
        raise MalformedSamplesheetError(f"{path}: no sample rows found")

    return rows


def build_sample_entries(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    """Checksum each row's FASTQ pair and build one manifest entry per read group."""
    fastq_checksums: dict[str, str] = {}
    entries: list[dict[str, object]] = []

    for row in rows:
        for key in ("fastq_1", "fastq_2"):
            fastq_path = Path(row[key])
            if str(fastq_path) not in fastq_checksums:
                fastq_checksums[str(fastq_path)] = f"sha256:{sha256_file(fastq_path)}"

        entries.append(
            {
                "sample_id": row["sample_id"],
                "read_group_id": row["read_group_id"],
                "library_id": row["library_id"],
                "platform": row["platform"],
                "platform_unit": row["platform_unit"],
                "fastq_1": {
                    "filename": Path(row["fastq_1"]).name,
                    "checksum": fastq_checksums[str(Path(row["fastq_1"]))],
                },
                "fastq_2": {
                    "filename": Path(row["fastq_2"]).name,
                    "checksum": fastq_checksums[str(Path(row["fastq_2"]))],
                },
            }
        )

    return entries


def read_metric_tsv(path: Path, required_metrics: tuple[str, ...]) -> dict[str, str]:
    """Read a `cohort_id\\tmetric\\tvalue` TSV (this pipeline's own QC/accounting shape).

    Returns every metric found, not only `required_metrics` -- the
    caller decides which subset is embedded in the manifest; this
    function's own job is only to fail loudly if a metric the caller
    said it needs is absent.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    metrics: dict[str, str] = {}
    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            raise MalformedAccountingError(
                f"{path}: row has {len(fields)} tab-separated fields, expected at least 3"
            )

        metrics[fields[1]] = fields[2]

    missing = [metric for metric in required_metrics if metric not in metrics]
    if missing:
        raise MalformedAccountingError(f"{path}: missing required metric(s): {missing}")

    return metrics


def read_variant_qc_tsv(path: Path, required_metrics: tuple[str, ...]) -> dict[str, str]:
    """Read a `summarize_variant_qc.py`-shaped TSV: `cohort_id\\tstage\\ttype\\tmetric\\tvalue`.

    Unlike `read_metric_tsv`'s simpler 3-column accounting shape, this
    file carries a stage/variant_type dimension per row (e.g. every row
    in `cohort.raw.all.variant_qc.tsv` is already `stage=raw,
    type=all`) -- the caller is expected to pass the single-stage/type
    file it wants (as this pipeline's own `--output` naming convention
    guarantees), not a multi-stage file mixing several.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()

    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    metrics: dict[str, str] = {}
    for line in lines[1:]:
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 5:
            raise MalformedAccountingError(
                f"{path}: row has {len(fields)} tab-separated fields, expected at least 5"
            )

        metrics[fields[3]] = fields[4]

    missing = [metric for metric in required_metrics if metric not in metrics]
    if missing:
        raise MalformedAccountingError(f"{path}: missing required metric(s): {missing}")

    return metrics


def read_gs_panel_manifest(path: Path) -> dict[str, object]:
    """Read an existing GS panel manifest, embedding only its own summary fields.

    The GS panel manifest already checksums its own artifacts; this
    function does not re-hash or duplicate that -- it embeds the GS
    manifest's own `manifest_hash` as a pointer, so a reader can verify
    the two manifests agree about which GS panel run this refers to
    without this script re-deriving anything the GS manifest already
    computed.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = ("schema_version", "run_id", "cohort_id", "panel_status", "manifest_hash")
    missing = [key for key in required if key not in payload]
    if missing:
        raise MalformedAccountingError(f"{path}: GS panel manifest missing key(s): {missing}")

    return {
        "schema_version": payload["schema_version"],
        "run_id": payload["run_id"],
        "cohort_id": payload["cohort_id"],
        "panel_status": payload["panel_status"],
        "manifest_hash": payload["manifest_hash"],
    }


class MalformedProvenanceError(Exception):
    """Raised when a DAG-staged provenance TSV is missing or malformed."""


class ProvenanceInconsistencyError(Exception):
    """Raised when two provenance inputs disagree about the same run."""


def _read_tsv_rows(path: Path, expected_fields: int) -> list[list[str]]:
    """Read a headerless TSV, failing on any row of the wrong width.

    Splitting on tabs and checking the width is what keeps a truncated
    or reshaped row from being read by position into the wrong field --
    which would otherwise put, say, a checksum where a filename belongs
    and record it without complaint.
    """
    rows: list[list[str]] = []
    for line_number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != expected_fields:
            raise MalformedProvenanceError(
                f"{path}:{line_number}: row has {len(fields)} tab-separated "
                f"fields, expected exactly {expected_fields}"
            )
        rows.append(fields)

    if not rows:
        raise MalformedProvenanceError(f"{path}: no rows found")

    return rows


def parse_runtime_provenance(paths: list[Path]) -> dict[str, str]:
    """Collapse per-task container identities into one entry per process.

    Every containerized process emits its own effective `task.container`
    once per task, so a process that ran for eight samples contributes
    eight identical rows; those collapse to a single entry.

    Two *different* identities under one process key is a different
    situation entirely, and this raises rather than picking one: it means
    some tasks of that process ran in one image and some in another (a
    per-task dynamic container, or a selector that matched only part of
    the run), and no single value would be a true statement about the
    run. Schema v2's contract is one effective container per process; a
    run that violates it must fail loudly rather than be recorded
    ambiguously.
    """
    identities_by_process: dict[str, set[str]] = {}

    for path in paths:
        for process_key, container in _read_tsv_rows(
            path, RUNTIME_PROVENANCE_FIELD_COUNT
        ):
            if not PROCESS_KEY_RE.match(process_key):
                raise MalformedProvenanceError(
                    f"{path}: '{process_key}' is not a canonical process key "
                    "(lowercase letters, digits and underscores, starting with a "
                    "letter)"
                )
            validate_container_identity(process_key, container)
            identities_by_process.setdefault(process_key, set()).add(container)

    if not identities_by_process:
        raise MalformedProvenanceError(
            "no runtime container provenance was recorded for this run"
        )

    divergent = {
        process_key: sorted(identities)
        for process_key, identities in identities_by_process.items()
        if len(identities) > 1
    }
    if divergent:
        details = "; ".join(
            f"{process_key}: {identities}" for process_key, identities in sorted(divergent.items())
        )
        raise ProvenanceInconsistencyError(
            "the same process ran in more than one container during this run, so "
            "no single effective container identity describes it "
            f"({details}). Schema v2 records one container per process; recording "
            "either value alone would be a false statement about the other tasks."
        )

    return {
        process_key: identities.pop()
        for process_key, identities in sorted(identities_by_process.items())
    }


def parse_input_provenance(paths: list[Path]) -> list[dict[str, object]]:
    """Restore samplesheet order from the per-read-group provenance rows.

    Nextflow makes no ordering promise across parallel tasks, so each row
    carries the zero-based samplesheet position it came from. That rank
    is used here to sort and then dropped: the manifest presents read
    groups in the order the samplesheet listed them, which is the order a
    reader can compare against their own input, without publishing an
    index that means nothing outside this run.
    """
    ranked: list[tuple[int, dict[str, object]]] = []
    seen_ranks: dict[int, str] = {}

    for path in paths:
        for fields in _read_tsv_rows(path, INPUT_PROVENANCE_FIELD_COUNT):
            (
                raw_rank,
                sample_id,
                read_group_id,
                library_id,
                platform,
                platform_unit,
                fastq_1_filename,
                fastq_1_checksum,
                fastq_2_filename,
                fastq_2_checksum,
            ) = fields

            try:
                rank = int(raw_rank)
            except ValueError as error:
                raise MalformedProvenanceError(
                    f"{path}: samplesheet rank is not an integer: {raw_rank!r}"
                ) from error

            if rank in seen_ranks:
                raise ProvenanceInconsistencyError(
                    f"two read groups claim samplesheet position {rank}: "
                    f"'{seen_ranks[rank]}' and '{read_group_id}'"
                )
            seen_ranks[rank] = read_group_id

            ranked.append(
                (
                    rank,
                    {
                        "sample_id": sample_id,
                        "read_group_id": read_group_id,
                        "library_id": library_id,
                        "platform": platform,
                        "platform_unit": platform_unit,
                        "fastq_1": {
                            "filename": fastq_1_filename,
                            "checksum": fastq_1_checksum,
                        },
                        "fastq_2": {
                            "filename": fastq_2_filename,
                            "checksum": fastq_2_checksum,
                        },
                    },
                )
            )

    return [entry for _rank, entry in sorted(ranked, key=lambda item: item[0])]


def parse_reference_provenance(path: Path) -> dict[str, object]:
    """Read the reference bundle's `role/filename/checksum` rows."""
    single: dict[str, dict[str, str]] = {}
    multi: dict[str, list[dict[str, str]]] = {role: [] for role in REFERENCE_MULTI_FILE_ROLES}

    for role, filename, checksum in _read_tsv_rows(
        path, REFERENCE_PROVENANCE_FIELD_COUNT
    ):
        entry = {"filename": filename, "checksum": checksum}
        if role in REFERENCE_SINGLE_FILE_ROLES:
            if role in single:
                raise ProvenanceInconsistencyError(
                    f"{path}: reference role '{role}' appears more than once"
                )
            single[role] = entry
        elif role in REFERENCE_MULTI_FILE_ROLES:
            multi[role].append(entry)
        else:
            raise MalformedProvenanceError(
                f"{path}: unknown reference role: {role!r}"
            )

    missing = [role for role in REFERENCE_SINGLE_FILE_ROLES if role not in single]
    if missing:
        raise MalformedProvenanceError(
            f"{path}: missing required reference role(s): {missing}"
        )
    for role in REFERENCE_MULTI_FILE_ROLES:
        expected_count = REFERENCE_MULTI_FILE_ROLE_COUNTS[role]
        if len(multi[role]) != expected_count:
            raise MalformedProvenanceError(
                f"{path}: reference role '{role}' must have exactly "
                f"{expected_count} files, found {len(multi[role])}"
            )
        # Sorted by filename so the manifest is byte-identical across
        # runs regardless of the order Nextflow staged the index files.
        multi[role] = sorted(multi[role], key=lambda entry: entry["filename"])

    return {**single, **multi}


def parse_artifact_checksums(paths: list[Path]) -> dict[str, str]:
    """Merge every artifact-checksum group into one filename-keyed map."""
    checksums: dict[str, str] = {}

    for path in paths:
        for filename, checksum in _read_tsv_rows(path, ARTIFACT_CHECKSUM_FIELD_COUNT):
            if filename in checksums and checksums[filename] != checksum:
                raise ProvenanceInconsistencyError(
                    f"two different artifacts were recorded under the same "
                    f"filename '{filename}'"
                )
            if filename in checksums:
                raise ValueError(f"duplicate checksum file name: '{filename}'")
            checksums[filename] = checksum

    if not checksums:
        raise MalformedProvenanceError("no run artifact checksums were recorded")

    return dict(sorted(checksums.items()))


def _read_header_and_rows(
    path: Path, header: tuple[str, ...]
) -> list[list[str]]:
    """Read a headed TSV whose header must match `header` exactly."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise MalformedAccountingError(f"{path}: file is empty")

    actual_header = tuple(lines[0].split("\t"))
    if actual_header != header:
        raise MalformedAccountingError(
            f"{path}: unexpected header {list(actual_header)}, expected {list(header)}"
        )

    rows: list[list[str]] = []
    for line_number, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            raise MalformedAccountingError(
                f"{path}:{line_number}: row has {len(fields)} tab-separated "
                f"fields, expected exactly {len(header)}"
            )
        rows.append(fields)

    if not rows:
        raise MalformedAccountingError(f"{path}: no data rows found")

    return rows


def read_cohort_accounting(path: Path, cohort_id: str) -> dict[str, str]:
    """Read the raw/all cohort accounting, validating every row's identity.

    Unlike the legacy `read_variant_qc_tsv`, which reads metrics by
    position and ignores the cohort/stage/type columns entirely, this
    checks them: a manifest that embedded the filtered-SNP QC under
    `cohort_accounting` while claiming raw/all would be wrong in a way no
    downstream reader could detect.
    """
    metrics: dict[str, str] = {}

    for cohort, stage, variant_type, metric, value in _read_header_and_rows(
        path, VARIANT_QC_HEADER
    ):
        if cohort != cohort_id:
            raise ProvenanceInconsistencyError(
                f"{path}: cohort accounting is for cohort '{cohort}', but this run "
                f"is '{cohort_id}'"
            )
        if (stage, variant_type) != (
            COHORT_ACCOUNTING_STAGE,
            COHORT_ACCOUNTING_VARIANT_TYPE,
        ):
            raise ProvenanceInconsistencyError(
                f"{path}: cohort accounting must come from the "
                f"{COHORT_ACCOUNTING_STAGE}/{COHORT_ACCOUNTING_VARIANT_TYPE} QC "
                f"run, found stage '{stage}' type '{variant_type}'"
            )
        metrics[metric] = value

    missing = [
        metric
        for metric in REQUIRED_COHORT_ACCOUNTING_METRICS
        if metric not in metrics
    ]
    if missing:
        raise MalformedAccountingError(
            f"{path}: missing required metric(s): {missing}"
        )

    return metrics


def read_sample_accounting(path: Path, cohort_id: str) -> list[dict[str, str]]:
    """Read summarize_variant_qc.py's raw/all sample_qc.tsv, in file order.

    Row order is the order that file already publishes (bcftools' own
    per-sample-counts order), preserved rather than re-sorted so the
    manifest and the QC artifact can be compared line by line.
    """
    samples: list[dict[str, str]] = []
    seen: set[str] = set()

    for row in _read_header_and_rows(path, SAMPLE_QC_HEADER):
        cohort, stage, variant_type = row[0], row[1], row[2]
        if cohort != cohort_id:
            raise ProvenanceInconsistencyError(
                f"{path}: sample accounting is for cohort '{cohort}', but this run "
                f"is '{cohort_id}'"
            )
        if (stage, variant_type) != (
            COHORT_ACCOUNTING_STAGE,
            COHORT_ACCOUNTING_VARIANT_TYPE,
        ):
            raise ProvenanceInconsistencyError(
                f"{path}: sample accounting must come from the "
                f"{COHORT_ACCOUNTING_STAGE}/{COHORT_ACCOUNTING_VARIANT_TYPE} QC "
                f"run, found stage '{stage}' type '{variant_type}'"
            )

        sample_name = row[3]
        if sample_name in seen:
            raise ProvenanceInconsistencyError(
                f"{path}: sample '{sample_name}' appears more than once"
            )
        seen.add(sample_name)

        samples.append(dict(zip(SAMPLE_ACCOUNTING_FIELDS, row[3:], strict=True)))

    return samples


def read_variant_type_accounting(path: Path, cohort_id: str) -> dict[str, str]:
    """Read RECONCILE_VARIANT_TYPE_COUNTS' accounting TSV for this cohort."""
    metrics: dict[str, str] = {}

    for cohort, metric, value in _read_header_and_rows(
        path, ("cohort_id", "metric", "value")
    ):
        if cohort != cohort_id:
            raise ProvenanceInconsistencyError(
                f"{path}: variant type accounting is for cohort '{cohort}', but "
                f"this run is '{cohort_id}'"
            )
        metrics[metric] = value

    missing = [
        metric
        for metric in REQUIRED_VARIANT_TYPE_ACCOUNTING_METRICS
        if metric not in metrics
    ]
    if missing:
        raise MalformedAccountingError(
            f"{path}: missing required metric(s): {missing}"
        )

    return metrics


def read_gs_panel_summary(path: Path, cohort_id: str) -> dict[str, object]:
    """Embed the GS panel manifest's own summary, checking it is this run's.

    The GS manifest already checksums its own artifacts; this embeds its
    `manifest_hash` as a pointer rather than re-deriving anything. The
    cohort check is what stops a stale GS manifest from a previous run
    being attached to this one -- the two documents would each be
    internally consistent and jointly wrong.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))

    missing = [key for key in REQUIRED_GS_PANEL_MANIFEST_KEYS if key not in payload]
    if missing:
        raise MalformedAccountingError(
            f"{path}: GS panel manifest missing key(s): {missing}"
        )

    if payload["cohort_id"] != cohort_id:
        raise ProvenanceInconsistencyError(
            f"{path}: GS panel manifest is for cohort '{payload['cohort_id']}', "
            f"but this run is '{cohort_id}'"
        )

    return {key: payload[key] for key in REQUIRED_GS_PANEL_MANIFEST_KEYS}


def cross_validate_accounting(
    *, cohort_accounting: dict[str, str], sample_accounting: list[dict[str, str]]
) -> None:
    """Refuse to record cohort and per-sample accounting that disagree.

    Both come from the same raw/all QC run, so the cohort's own
    `number_of_samples`/`sample_names` must describe exactly the rows in
    the per-sample table. A provenance builder that recorded two
    contradictory statements about the same run without complaint would
    be worse than one that recorded nothing.
    """
    recorded_names = [sample["sample"] for sample in sample_accounting]

    try:
        expected_count = int(cohort_accounting["number_of_samples"])
    except ValueError as error:
        raise MalformedAccountingError(
            "cohort accounting number_of_samples is not an integer: "
            f"{cohort_accounting['number_of_samples']!r}"
        ) from error

    if expected_count != len(recorded_names):
        raise ProvenanceInconsistencyError(
            f"cohort accounting reports {expected_count} samples, but per-sample "
            f"accounting has {len(recorded_names)} rows"
        )

    expected_names = [
        name for name in cohort_accounting["sample_names"].split(",") if name
    ]
    if expected_names != recorded_names:
        raise ProvenanceInconsistencyError(
            f"cohort accounting sample_names {expected_names} do not match the "
            f"per-sample accounting rows {recorded_names}, in order"
        )


def build_dag_manifest(
    *,
    cohort_id: str,
    pipeline_version: str,
    git_commit: str,
    nextflow_version: str,
    containers: dict[str, str],
    reference: dict[str, object],
    parameters: dict[str, object],
    samples: list[dict[str, object]],
    cohort_accounting: dict[str, str],
    sample_accounting: list[dict[str, str]],
    variant_type_accounting: dict[str, str],
    gs_panel: dict[str, object] | None,
    checksums: dict[str, str],
    run_id: str,
    generated_at: str,
) -> dict[str, object]:
    """Build the schema v2 (DAG-generated) manifest, including its own hash."""
    manifest: dict[str, object] = {
        "schema_version": DAG_SCHEMA_VERSION,
        "run_id": run_id,
        "generated_at": generated_at,
        "cohort_id": cohort_id,
        "pipeline_version": pipeline_version,
        "git_commit": git_commit,
        "nextflow_version": nextflow_version,
        "containers": containers,
        "reference": reference,
        "parameters": parameters,
        "samples": samples,
        "cohort_accounting": cohort_accounting,
        "sample_accounting": sample_accounting,
        "variant_type_accounting": variant_type_accounting,
        "gs_panel": gs_panel,
        "checksums": checksums,
    }
    manifest["manifest_hash"] = canonical_json_hash(manifest)
    return manifest


def build_manifest(
    *,
    cohort_id: str,
    pipeline_version: str,
    git_commit: str,
    nextflow_version: str,
    containers: dict[str, str],
    reference: dict[str, object],
    parameters: dict[str, object],
    samples: list[dict[str, object]],
    cohort_accounting: dict[str, str],
    variant_type_accounting: dict[str, str],
    gs_panel: dict[str, object] | None,
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
        "nextflow_version": nextflow_version,
        "containers": containers,
        "reference": reference,
        "parameters": parameters,
        "samples": samples,
        "cohort_accounting": cohort_accounting,
        "variant_type_accounting": variant_type_accounting,
        "gs_panel": gs_panel,
        "checksums": checksums,
    }
    manifest["manifest_hash"] = canonical_json_hash(manifest)
    return manifest


def _peek_mode(argv: list[str] | None) -> str:
    """Resolve --mode before building the real parser.

    The two modes require genuinely different inputs, and expressing
    that as one parser with everything optional would replace argparse's
    own "the following arguments are required" errors with silent
    `None`s discovered much later. Parsing the mode first lets each mode
    declare its own required set -- and, importantly, leaves the legacy
    v1 invocation's argparse behavior exactly as it was before this
    script grew a second mode.
    """
    mode_parser = argparse.ArgumentParser(add_help=False)
    mode_parser.add_argument("--mode", choices=(LEGACY_MODE, DAG_MODE), default=LEGACY_MODE)
    known, _unknown = mode_parser.parse_known_args(argv)
    return known.mode


def _add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the arguments both modes require, with identical meanings."""
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument(
        "--pipeline-version", required=True, help="workflow.manifest.version"
    )
    parser.add_argument(
        "--nextflow-version",
        required=True,
        help="The Nextflow execution engine's own version (e.g. 26.04.6), not a "
        "container-pinned tool version.",
    )
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--reference-species", required=True)
    parser.add_argument("--reference-cultivar", required=True)
    parser.add_argument("--reference-accession", required=True)

    parser.add_argument("--sample-ploidy", required=True, type=int)
    parser.add_argument("--genomicsdb-batch-size", required=True, type=int)
    parser.add_argument("--optical-duplicate-pixel-distance", required=True, type=int)
    parser.add_argument(
        "--enable-gs-panel", required=True, action=argparse.BooleanOptionalAction
    )
    for name in (*SNP_FILTER_PARAM_NAMES, *INDEL_FILTER_PARAM_NAMES):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True, type=float)

    parser.add_argument(
        "--variant-qc-tsv",
        required=True,
        type=Path,
        help="Path to cohort.raw.all.variant_qc.tsv.",
    )
    parser.add_argument(
        "--variant-type-accounting-tsv",
        required=True,
        type=Path,
        help="Path to cohort.variant_type_accounting.tsv.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Output path for the manifest JSON."
    )


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the whole-run manifest-builder CLI."""
    mode = _peek_mode(argv)

    parser = argparse.ArgumentParser(
        description="Build a whole-run reproducibility manifest for an "
        "adzuki-snp-pipeline execution."
    )
    parser.add_argument(
        "--mode",
        choices=(LEGACY_MODE, DAG_MODE),
        default=LEGACY_MODE,
        help=(
            f"'{LEGACY_MODE}' (the default) builds the schema v1 document this "
            "script has always built, from a completed run's published artifacts, "
            f"invoked by hand. '{DAG_MODE}' builds the schema v2 document the "
            "Nextflow workflow builds for itself, from provenance staged by the "
            "DAG (Issue #42). The two record container identity differently and "
            "are therefore different schema versions, not two spellings of one."
        ),
    )
    _add_shared_arguments(parser)

    if mode == LEGACY_MODE:
        parser.add_argument(
            "--git-commit",
            default="",
            help="Nextflow workflow.commitId; empty string is recorded as null.",
        )
        parser.add_argument("--bwa-mem2-container", required=True)
        parser.add_argument("--samtools-container", required=True)
        parser.add_argument("--gatk-container", required=True)
        parser.add_argument("--python-container", required=True)
        parser.add_argument("--reference-fasta", required=True, type=Path)
        parser.add_argument("--reference-fai", required=True, type=Path)
        parser.add_argument("--reference-dict", required=True, type=Path)
        parser.add_argument(
            "--samplesheet",
            required=True,
            type=Path,
            help="Path to the run's samplesheet.csv.",
        )
        parser.add_argument(
            "--gs-panel-manifest",
            type=Path,
            default=None,
            help="Path to an existing gs_panel manifest.json, if the GS panel ran.",
        )
        parser.add_argument(
            "--checksum-file",
            action="append",
            default=[],
            type=Path,
            help="An additional run artifact to checksum (gVCFs, cohort VCF, ...); "
            "repeatable.",
        )
        return parser.parse_args(argv)

    parser.add_argument(
        "--git-commit",
        required=True,
        help="The full 40-character git SHA the run was executed from. Unlike "
        "legacy mode, this is required and must be a full SHA: an abbreviated or "
        "absent commit cannot identify the code that produced a delivered result.",
    )
    parser.add_argument("--reference-name", required=True)
    parser.add_argument(
        "--runtime-provenance",
        required=True,
        action="append",
        type=Path,
        help="A `process_key<TAB>container` TSV of effective task.container "
        "values, one row per executed task; repeatable.",
    )
    parser.add_argument(
        "--input-provenance",
        required=True,
        action="append",
        type=Path,
        help="A hash_input_fastqs.py row for one read group; repeatable.",
    )
    parser.add_argument(
        "--reference-provenance",
        required=True,
        type=Path,
        help="hash_reference_bundle.py's `role/filename/checksum` TSV.",
    )
    parser.add_argument(
        "--artifact-checksums",
        required=True,
        action="append",
        type=Path,
        help="A hash_run_artifacts.py `filename/checksum` TSV; repeatable.",
    )
    parser.add_argument(
        "--sample-qc-tsv",
        required=True,
        type=Path,
        help="Path to cohort.raw.all.sample_qc.tsv (summarize_variant_qc.py's "
        "existing per-sample output; this script never recomputes it).",
    )
    gs_group = parser.add_mutually_exclusive_group(required=True)
    gs_group.add_argument(
        "--gs-panel-manifest",
        type=Path,
        default=None,
        help="The GS panel manifest this run produced.",
    )
    gs_group.add_argument(
        "--no-gs-panel",
        action="store_true",
        help="This run had the GS panel disabled; `gs_panel` is recorded as null. "
        "Required explicitly rather than inferred from a missing argument, so a "
        "wiring mistake cannot quietly produce a manifest that understates what "
        "the run did.",
    )
    return parser.parse_args(argv)


def _build_parameters(args: argparse.Namespace) -> dict[str, object]:
    """Collect the material parameters both schema versions record."""
    return {
        "sample_ploidy": args.sample_ploidy,
        "genomicsdb_batch_size": args.genomicsdb_batch_size,
        "optical_duplicate_pixel_distance": args.optical_duplicate_pixel_distance,
        "enable_gs_panel": args.enable_gs_panel,
        **{name: getattr(args, name) for name in SNP_FILTER_PARAM_NAMES},
        **{name: getattr(args, name) for name in INDEL_FILTER_PARAM_NAMES},
    }


def _run_dag_mode(args: argparse.Namespace) -> int:
    """Build the schema v2 manifest from DAG-staged provenance."""
    if not FULL_GIT_SHA_RE.match(args.git_commit):
        print(
            "build_run_manifest.py: error: --git-commit must be a full 40-character "
            f"git SHA in {DAG_MODE} mode, got {args.git_commit!r}. Nextflow's own "
            "workflow.commitId is null for a local `nextflow run .`, so the "
            "workflow resolves the commit explicitly; an unresolved or abbreviated "
            "value would leave the manifest unable to identify the code that "
            "produced this run.",
            file=sys.stderr,
        )
        return 1

    try:
        containers = parse_runtime_provenance(args.runtime_provenance)
        samples = parse_input_provenance(args.input_provenance)
        reference_files = parse_reference_provenance(args.reference_provenance)
        checksums = parse_artifact_checksums(args.artifact_checksums)
        cohort_accounting = read_cohort_accounting(args.variant_qc_tsv, args.cohort_id)
        sample_accounting = read_sample_accounting(args.sample_qc_tsv, args.cohort_id)
        variant_type_accounting = read_variant_type_accounting(
            args.variant_type_accounting_tsv, args.cohort_id
        )
        cross_validate_accounting(
            cohort_accounting=cohort_accounting, sample_accounting=sample_accounting
        )
        gs_panel = (
            None
            if args.no_gs_panel
            else read_gs_panel_summary(args.gs_panel_manifest, args.cohort_id)
        )
    except OSError as error:
        print(f"build_run_manifest.py: error: {error}", file=sys.stderr)
        return 1
    except (
        MalformedAccountingError,
        MalformedProvenanceError,
        ProvenanceInconsistencyError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        print(f"build_run_manifest.py: error: {error}", file=sys.stderr)
        return 1

    reference: dict[str, object] = {
        "reference_id": args.reference_id,
        "reference_name": args.reference_name,
        "reference_species": args.reference_species,
        "reference_cultivar": args.reference_cultivar,
        "reference_accession": args.reference_accession,
        **reference_files,
    }

    manifest = build_dag_manifest(
        cohort_id=args.cohort_id,
        pipeline_version=args.pipeline_version,
        git_commit=args.git_commit,
        nextflow_version=args.nextflow_version,
        containers=containers,
        reference=reference,
        parameters=_build_parameters(args),
        samples=samples,
        cohort_accounting=cohort_accounting,
        sample_accounting=sample_accounting,
        variant_type_accounting=variant_type_accounting,
        gs_panel=gs_panel,
        checksums=checksums,
        run_id=new_run_id(),
        generated_at=utc_now_iso(),
    )

    # The payload is assembled from basenames, checksums, scientific
    # identifiers and numbers only -- nothing host-specific is passed in.
    # This is the backstop for that discipline, not a substitute for it.
    try:
        assert_no_host_metadata(manifest)
    except HostMetadataLeakError as error:
        print(f"build_run_manifest.py: error: {error}", file=sys.stderr)
        return 1

    write_json_atomic(args.output, manifest)
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    if args.mode == DAG_MODE:
        return _run_dag_mode(args)

    return _run_legacy_mode(args)


def _run_legacy_mode(args: argparse.Namespace) -> int:
    """Build the schema v1 manifest from a completed run's published artifacts.

    Unchanged from before Issue #42 added a second mode: this is the
    contract the historical Issue #26/#33 real-cohort manifests were
    produced under, and those documents are not being regenerated.
    """
    try:
        rows = parse_samplesheet(args.samplesheet)
        samples = build_sample_entries(rows)

        reference_checksums = checksum_files(
            [args.reference_fasta, args.reference_fai, args.reference_dict]
        )
        reference = {
            "reference_id": args.reference_id,
            "reference_species": args.reference_species,
            "reference_cultivar": args.reference_cultivar,
            "reference_accession": args.reference_accession,
            "fasta": {
                "filename": args.reference_fasta.name,
                "checksum": reference_checksums[args.reference_fasta.name],
            },
            "fai": {
                "filename": args.reference_fai.name,
                "checksum": reference_checksums[args.reference_fai.name],
            },
            "dict": {
                "filename": args.reference_dict.name,
                "checksum": reference_checksums[args.reference_dict.name],
            },
        }

        cohort_metrics = read_variant_qc_tsv(
            args.variant_qc_tsv,
            required_metrics=(
                "number_of_samples",
                "sample_names",
                "cohort_total_genotypes",
                "cohort_missing_genotypes",
                "cohort_missingness_rate",
            ),
        )
        variant_type_metrics = read_metric_tsv(
            args.variant_type_accounting_tsv,
            required_metrics=("raw_all_records", "raw_snp_records", "raw_indel_records"),
        )

        gs_panel = (
            read_gs_panel_manifest(args.gs_panel_manifest)
            if args.gs_panel_manifest is not None
            else None
        )

        checksums = checksum_files(args.checksum_file)
    except OSError as error:
        print(f"build_run_manifest.py: error: {error}", file=sys.stderr)
        return 1
    except (
        MalformedSamplesheetError,
        MalformedAccountingError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        print(f"build_run_manifest.py: error: {error}", file=sys.stderr)
        return 1

    manifest = build_manifest(
        cohort_id=args.cohort_id,
        pipeline_version=args.pipeline_version,
        git_commit=args.git_commit,
        nextflow_version=args.nextflow_version,
        containers={
            "bwa_mem2": args.bwa_mem2_container,
            "samtools": args.samtools_container,
            "gatk": args.gatk_container,
            "python": args.python_container,
        },
        reference=reference,
        parameters=_build_parameters(args),
        samples=samples,
        cohort_accounting=cohort_metrics,
        variant_type_accounting=variant_type_metrics,
        gs_panel=gs_panel,
        checksums=checksums,
        run_id=new_run_id(),
        generated_at=utc_now_iso(),
    )

    write_json_atomic(args.output, manifest)

    return 0


if __name__ == "__main__":
    sys.exit(main())
