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

This is a standalone CLI, not wired into the Nextflow workflow: it is
run once, by hand, against a completed run's own already-published
artifacts (samplesheet, reference bundle, QC/accounting TSVs, and
whichever files the caller wants checksummed via repeatable
`--checksum-file`). It does not read raw FASTQ content beyond hashing
the files the samplesheet already points at, does not re-run or
second-guess any GATK/samtools/bcftools computation, and does not
decide which files belong in a delivery package -- that judgment call
(this repository's `analysis/report.py`-equivalent for this pipeline
does not yet exist) is left to whoever invokes this script.

Software versions are recorded as the pinned container image
references already baked into each process's `container` directive
(the ground truth for "what actually ran" in a Nextflow-plus-Docker
pipeline), passed in via `--bwa-mem2-container`/`--samtools-container`/
`--gatk-container`/`--python-container` rather than queried at run
time -- the same reasoning as `build_gs_panel_manifest.py`.

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
import sys
from pathlib import Path

# Issue #42: the hashing/serialization mechanics this manifest and
# bin/build_gs_panel_manifest.py both depend on live in one module now
# (they were previously duplicated character for character). See
# bin/manifest_utils.py for why a plain sibling import needs no
# packaging or PYTHONPATH in either place these scripts run.
from manifest_utils import (
    _json_default,
    canonical_json_hash,
    checksum_files,
    new_run_id,
    sha256_file,
    utc_now_iso,
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


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the whole-run manifest-builder CLI."""
    parser = argparse.ArgumentParser(
        description="Build a whole-run reproducibility manifest for a completed "
        "adzuki-snp-pipeline execution."
    )
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument(
        "--pipeline-version", required=True, help="workflow.manifest.version"
    )
    parser.add_argument(
        "--git-commit",
        default="",
        help="Nextflow workflow.commitId; empty string is recorded as null.",
    )
    parser.add_argument(
        "--nextflow-version",
        required=True,
        help="The Nextflow execution engine's own version (e.g. 26.04.6), not a "
        "container-pinned tool version.",
    )

    parser.add_argument("--bwa-mem2-container", required=True)
    parser.add_argument("--samtools-container", required=True)
    parser.add_argument("--gatk-container", required=True)
    parser.add_argument("--python-container", required=True)

    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--reference-species", required=True)
    parser.add_argument("--reference-cultivar", required=True)
    parser.add_argument("--reference-accession", required=True)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-fai", required=True, type=Path)
    parser.add_argument("--reference-dict", required=True, type=Path)

    parser.add_argument("--sample-ploidy", required=True, type=int)
    parser.add_argument("--genomicsdb-batch-size", required=True, type=int)
    parser.add_argument("--optical-duplicate-pixel-distance", required=True, type=int)
    parser.add_argument(
        "--enable-gs-panel", required=True, action=argparse.BooleanOptionalAction
    )
    for name in (*SNP_FILTER_PARAM_NAMES, *INDEL_FILTER_PARAM_NAMES):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True, type=float)

    parser.add_argument(
        "--samplesheet", required=True, type=Path, help="Path to the run's samplesheet.csv."
    )
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
        help="An additional run artifact to checksum (gVCFs, cohort VCF, ...); repeatable.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Output path for the manifest JSON."
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

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

    parameters: dict[str, object] = {
        "sample_ploidy": args.sample_ploidy,
        "genomicsdb_batch_size": args.genomicsdb_batch_size,
        "optical_duplicate_pixel_distance": args.optical_duplicate_pixel_distance,
        "enable_gs_panel": args.enable_gs_panel,
        **{name: getattr(args, name) for name in SNP_FILTER_PARAM_NAMES},
        **{name: getattr(args, name) for name in INDEL_FILTER_PARAM_NAMES},
    }

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
        parameters=parameters,
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
