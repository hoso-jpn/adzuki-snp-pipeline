#!/usr/bin/env python3
"""Validate 50+ homogeneous gVCFs and prepare an auditable benchmark plan."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from summarize_filter_qc import staged_qc_outputs

MANIFEST_COLUMNS = ("sample_id", "gvcf", "gvcf_index")
PIPELINE_COMMIT = re.compile(r"^[0-9a-f]{40}$")
CONTIG_HEADER = re.compile(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,|>)")

INTERVAL_PLAN_HEADER = (
    "plan_id",
    "interval_group_id",
    "interval_order",
    "intervals_json",
    "total_bp",
    "source_contig_count",
    "strategy",
)


class BenchmarkInputError(Exception):
    """Raised when benchmark inputs cannot support a defensible comparison."""


@dataclass(frozen=True)
class GvcfInput:
    sample_id: str
    gvcf: Path
    gvcf_index: Path
    gvcf_sha256: str
    gvcf_index_sha256: str


@dataclass(frozen=True)
class Contig:
    name: str
    length: int


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_reference_fai(path: Path) -> tuple[Contig, ...]:
    """Read ordered unique contig names and lengths from a reference FAI."""
    contigs: list[Contig] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise BenchmarkInputError(
                    f"{path}: line {line_number}: expected at least 2 tab-separated fields"
                )
            name = fields[0]
            if name in seen:
                raise BenchmarkInputError(f"{path}: duplicate contig: {name}")
            try:
                length = int(fields[1])
            except ValueError as error:
                raise BenchmarkInputError(
                    f"{path}: line {line_number}: contig length is not an integer"
                ) from error
            if length <= 0:
                raise BenchmarkInputError(
                    f"{path}: line {line_number}: contig length must be positive"
                )
            seen.add(name)
            contigs.append(Contig(name, length))
    if not contigs:
        raise BenchmarkInputError(f"{path}: no contigs found")
    return tuple(contigs)


def read_gvcf_identity(path: Path) -> tuple[str, tuple[Contig, ...]]:
    """Read the single sample and ordered contig dictionary from a gVCF header."""
    contigs: list[Contig] = []
    sample_id: str | None = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##contig="):
                match = CONTIG_HEADER.match(line.rstrip("\n"))
                if match is None:
                    raise BenchmarkInputError(
                        f"{path}: unsupported ##contig header: {line.strip()}"
                    )
                contigs.append(Contig(match.group(1), int(match.group(2))))
            elif line.startswith("#CHROM\t"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                if len(samples) != 1:
                    raise BenchmarkInputError(
                        f"{path}: expected exactly one sample column, found {len(samples)}"
                    )
                sample_id = samples[0]
                break
            elif not line.startswith("#"):
                raise BenchmarkInputError(f"{path}: data row appeared before #CHROM header")
    if sample_id is None:
        raise BenchmarkInputError(f"{path}: missing #CHROM header")
    if not contigs:
        raise BenchmarkInputError(f"{path}: no ##contig headers found")
    if len({contig.name for contig in contigs}) != len(contigs):
        raise BenchmarkInputError(f"{path}: duplicate ##contig ID")
    return sample_id, tuple(contigs)


def load_inputs(
    manifest: Path,
    reference_contigs: tuple[Contig, ...],
    minimum_samples: int,
) -> tuple[GvcfInput, ...]:
    """Validate manifest paths, sample names, indexes, and reference identity."""
    inputs: list[GvcfInput] = []
    sample_ids: set[str] = set()
    gvcf_paths: set[Path] = set()
    index_paths: set[Path] = set()
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(MANIFEST_COLUMNS):
            raise BenchmarkInputError(
                f"{manifest}: columns must be exactly {', '.join(MANIFEST_COLUMNS)}"
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(row[name] is None for name in MANIFEST_COLUMNS):
                raise BenchmarkInputError(
                    f"{manifest}: line {line_number}: expected exactly 3 fields"
                )
            sample_id = row["sample_id"]
            if not sample_id or any(character.isspace() for character in sample_id):
                raise BenchmarkInputError(f"{manifest}: line {line_number}: invalid sample_id")
            if sample_id in sample_ids:
                raise BenchmarkInputError(f"{manifest}: duplicate sample_id: {sample_id}")
            paths = []
            for name in ("gvcf", "gvcf_index"):
                value = row[name]
                if not value or any(character in value for character in "\t\r\n"):
                    raise BenchmarkInputError(
                        f"{manifest}: line {line_number}: invalid {name} path"
                    )
                path = Path(value).expanduser()
                paths.append((path if path.is_absolute() else manifest.parent / path).resolve())
            gvcf, gvcf_index = paths
            if gvcf in gvcf_paths or gvcf_index in index_paths:
                raise BenchmarkInputError(
                    f"{manifest}: line {line_number}: gVCF or index path is reused"
                )
            if not gvcf.is_file() or not gvcf_index.is_file():
                raise BenchmarkInputError(
                    f"{manifest}: line {line_number}: gVCF and index must both exist"
                )
            # The index is handed to GenomicsDBImport as the third
            # sample-name-map column, so a row pairing one sample's gVCF
            # with another sample's existing .tbi is accepted by every
            # check above -- both files exist, both paths are unique, and
            # the header sample still matches. That pairing silently
            # attributes a benchmark measurement to the wrong input
            # lineage, which is the one thing this preparation gate exists
            # to prevent, so require the index this pipeline actually
            # emits for that gVCF.
            expected_index = Path(str(gvcf) + ".tbi")
            if gvcf_index != expected_index:
                raise BenchmarkInputError(
                    f"{manifest}: line {line_number}: gvcf_index must be the gVCF's own "
                    f"index {expected_index.name}, not {gvcf_index.name}"
                )
            header_sample, gvcf_contigs = read_gvcf_identity(gvcf)
            if header_sample != sample_id:
                raise BenchmarkInputError(
                    f"{gvcf}: header sample {header_sample!r} does not match manifest {sample_id!r}"
                )
            if gvcf_contigs != reference_contigs:
                raise BenchmarkInputError(
                    f"{gvcf}: ordered contig dictionary does not match reference FAI"
                )

            sample_ids.add(sample_id)
            gvcf_paths.add(gvcf)
            index_paths.add(gvcf_index)
            inputs.append(
                GvcfInput(
                    sample_id=sample_id,
                    gvcf=gvcf,
                    gvcf_index=gvcf_index,
                    gvcf_sha256=_sha256(gvcf),
                    gvcf_index_sha256=_sha256(gvcf_index),
                )
            )

    if len(inputs) < minimum_samples:
        raise BenchmarkInputError(
            f"{manifest}: found {len(inputs)} samples, need at least {minimum_samples} "
            "to exercise 50+ batching"
        )
    return tuple(inputs)


def build_interval_plan_rows(
    contigs: tuple[Contig, ...],
    window_size_bp: int,
    small_scaffold_max_bp: int,
) -> list[list[str]]:
    """Build baseline and split/group candidate plans without selecting a winner."""
    rows: list[list[str]] = []
    for order, contig in enumerate(contigs, start=1):
        rows.append(
            [
                "baseline_per_contig",
                f"baseline_{order:04d}_{contig.name}",
                str(order),
                json.dumps([contig.name], separators=(",", ":")),
                str(contig.length),
                "1",
                "one FAI row per task; current production behavior",
            ]
        )

    candidate_groups: list[tuple[str, list[str], int, int, str]] = []
    small: list[Contig] = []
    small_group_count = 0
    small_bp = 0

    def flush_small() -> None:
        nonlocal small_group_count, small_bp
        if not small:
            return
        small_group_count += 1
        group_id = (
            "small_scaffolds" if small_group_count == 1 else f"small_scaffolds_{small_group_count}"
        )
        candidate_groups.append(
            (
                group_id,
                [contig.name for contig in small],
                sum(contig.length for contig in small),
                len(small),
                "consecutive small scaffolds grouped in reference dictionary order",
            )
        )
        small.clear()
        small_bp = 0

    for contig in contigs:
        if contig.length <= small_scaffold_max_bp:
            if small_bp + contig.length > window_size_bp:
                flush_small()
            small.append(contig)
            small_bp += contig.length
            continue
        flush_small()
        if contig.length <= window_size_bp:
            candidate_groups.append(
                (contig.name, [contig.name], contig.length, 1, "unsplit contig")
            )
            continue
        start = 1
        window = 1
        while start <= contig.length:
            end = min(start + window_size_bp - 1, contig.length)
            interval = f"{contig.name}:{start}-{end}"
            candidate_groups.append(
                (
                    f"{contig.name}_window_{window:03d}",
                    [interval],
                    end - start + 1,
                    1,
                    "fixed-size chromosome window",
                )
            )
            start = end + 1
            window += 1

    flush_small()

    for order, (group_id, intervals, total_bp, source_count, strategy) in enumerate(
        candidate_groups, start=1
    ):
        rows.append(
            [
                "candidate_split_group",
                group_id,
                str(order),
                json.dumps(intervals, separators=(",", ":")),
                str(total_bp),
                str(source_count),
                strategy,
            ]
        )
    split_rows = []
    for row in rows:
        if row[0] != "candidate_split_group":
            continue
        intervals = json.loads(row[3])
        if int(row[5]) == 1:
            split_rows.append(["candidate_split_only", *row[1:]])
        else:
            lengths = {contig.name: contig.length for contig in contigs}
            for name in intervals:
                split_rows.append(
                    [
                        "candidate_split_only",
                        name,
                        "",
                        json.dumps([name]),
                        str(lengths[name]),
                        "1",
                        "ungrouped small scaffold",
                    ]
                )
    for order, row in enumerate(split_rows, start=1):
        row[2] = str(order)
    rows.extend(split_rows)
    return rows


def build_evidence_template(
    args: argparse.Namespace,
    inputs: tuple[GvcfInput, ...],
    reference_contigs: tuple[Contig, ...],
) -> dict[str, object]:
    """Build a machine-readable record with measurements deliberately unset."""
    experiments = [
        ("E0", None, "repeated_variant_args", "baseline_per_contig", False, False),
        ("E1", "E0", "sample_name_map", "baseline_per_contig", False, False),
        ("E2a", "E1", "sample_name_map", "candidate_split_only", False, False),
        ("E2b", "E2a", "sample_name_map", "candidate_split_group", False, False),
        ("E3", "E2b", "sample_name_map", "candidate_split_group", True, False),
        ("E4", "E2b", "sample_name_map", "candidate_split_group", False, True),
    ]
    return {
        "schema_version": 1,
        "issue": 45,
        "status": "PENDING_REAL_BENCHMARK",
        "pipeline_commit": args.pipeline_commit,
        "gatk_container": args.gatk_container,
        "input_validation": {
            "benchmark_ready": False,
            "lineage_verified": False,
            "declared_lineage_source": "CLI arguments; not verified against source run evidence",
            "required_before_execution": [
                "source-run manifests tying every gVCF checksum to pipeline SHA and GATK version",
                "reference FASTA identity; FAI names/lengths/order alone do not identify sequence",
                "index-content validation; matching .tbi filenames alone do not prove pairing",
                "homogeneous ploidy and HaplotypeCaller parameters",
            ],
        },
        "reference": {
            "fai_name": args.reference_fai.name,
            "fai_sha256": _sha256(args.reference_fai),
            "contig_count": len(reference_contigs),
            "total_bp": sum(contig.length for contig in reference_contigs),
        },
        "inputs": {
            "sample_count": len(inputs),
            "samples": [
                {
                    "sample_id": item.sample_id,
                    "gvcf_name": item.gvcf.name,
                    "gvcf_sha256": item.gvcf_sha256,
                    "gvcf_index_name": item.gvcf_index.name,
                    "gvcf_index_sha256": item.gvcf_index_sha256,
                }
                for item in inputs
            ],
        },
        "fixed_parameters": {
            "batch_size": args.batch_size,
            "window_size_bp": args.window_size_bp,
            "small_scaffold_max_bp": args.small_scaffold_max_bp,
        },
        "experiments": [
            {
                "experiment_id": experiment_id,
                "compare_to": compare_to,
                "gvcf_input_mode": input_mode,
                "interval_plan_id": interval_plan,
                "reblock_gvcfs": reblock,
                "consolidate": consolidate,
                "measurements": {
                    "wall_time_seconds": None,
                    "peak_rss_bytes": None,
                    "workspace_bytes": None,
                    "workspace_file_count": None,
                    "genotypegvcfs_wall_time_seconds": None,
                    "gathervcfs_wall_time_seconds": None,
                    "retry_count": None,
                    "swap_delta_bytes": None,
                    "output_sample_count": None,
                    "output_sample_order_sha256": None,
                    "variant_accounting_sha256": None,
                    "failure_count": None,
                    "nextflow_version": None,
                    "gatk_version": None,
                },
            }
            for experiment_id, compare_to, input_mode, interval_plan, reblock, consolidate in experiments
        ],
        "decision": {
            "sample_name_map": "PENDING",
            "interval_strategy": "PENDING",
            "reblock_gvcfs": "PENDING",
            "consolidate": "PENDING",
            "full_327_sample_gate": "PENDING",
            "rationale": None,
        },
    }


def _write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _write_sample_name_map(path: Path, inputs: tuple[GvcfInput, ...]) -> None:
    """Write GATK's headerless sample, gVCF, optional-index three-column form."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        # GATK splits literal tabs; it does not decode CSV quoting. Input
        # validation rejects delimiters/newlines, so quotes remain literal.
        for item in inputs:
            handle.write(f"{item.sample_id}\t{item.gvcf}\t{item.gvcf_index}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gvcf-manifest", required=True, type=Path)
    parser.add_argument("--reference-fai", required=True, type=Path)
    parser.add_argument("--pipeline-commit", required=True)
    parser.add_argument("--gatk-container", required=True)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--window-size-bp", type=int, default=20_000_000)
    parser.add_argument("--small-scaffold-max-bp", type=int, default=1_000_000)
    parser.add_argument("--sample-name-map-output", required=True, type=Path)
    parser.add_argument("--interval-plan-output", required=True, type=Path)
    parser.add_argument("--evidence-template-output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if not PIPELINE_COMMIT.fullmatch(args.pipeline_commit):
            raise BenchmarkInputError("--pipeline-commit must be a full lowercase 40-hex SHA")
        if not re.fullmatch(r"[^\s]+@sha256:[0-9a-f]{64}", args.gatk_container):
            raise BenchmarkInputError(
                "--gatk-container must have a complete lowercase SHA256 digest"
            )
        for name in ("batch_size", "window_size_bp", "small_scaffold_max_bp"):
            if getattr(args, name) <= 0:
                raise BenchmarkInputError(f"--{name.replace('_', '-')} must be positive")
        # The candidate plan classifies a contig as a small scaffold
        # before it considers windowing it. With the two bounds inverted,
        # every contig -- including chromosome-scale ones -- falls into
        # the single grouped "small_scaffolds" task, so E2 would compare
        # a degenerate one-interval plan while the protocol table still
        # calls it "chromosome split + small-scaffold group". Reject the
        # inversion instead of silently benchmarking a different design.
        if args.small_scaffold_max_bp >= args.window_size_bp:
            raise BenchmarkInputError(
                "--small-scaffold-max-bp must be smaller than --window-size-bp; "
                f"got {args.small_scaffold_max_bp} >= {args.window_size_bp}"
            )

        reference_contigs = read_reference_fai(args.reference_fai)
        minimum_samples = max(51, args.batch_size + 1)
        inputs = load_inputs(args.gvcf_manifest, reference_contigs, minimum_samples)
        interval_rows = build_interval_plan_rows(
            reference_contigs,
            args.window_size_bp,
            args.small_scaffold_max_bp,
        )
        evidence = build_evidence_template(args, inputs, reference_contigs)
        outputs = (
            args.sample_name_map_output,
            args.interval_plan_output,
            args.evidence_template_output,
        )
        source_paths = (args.gvcf_manifest, args.reference_fai) + tuple(
            path for item in inputs for path in (item.gvcf, item.gvcf_index)
        )
        with staged_qc_outputs(outputs, inputs=source_paths) as staged:
            _write_sample_name_map(staged[0], inputs)
            _write_tsv(staged[1], INTERVAL_PLAN_HEADER, interval_rows)
            staged[2].write_text(
                json.dumps(evidence, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    except (OSError, EOFError, UnicodeError, ValueError, BenchmarkInputError) as error:
        print(f"prepare_joint_genotyping_benchmark.py: error: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
