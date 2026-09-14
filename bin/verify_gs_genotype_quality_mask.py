#!/usr/bin/env python3
"""Verify that every quality-masked GS artifact describes one masked state (Issue #64).

`build_gs_panel.py` writes the matrix, the metadata, the genotype accounting
and the quality-masked VCF from a single pass. Being written together is not
evidence that they agree, so this reads all of them back -- together with the
original GS-eligible PASS VCF they came from -- and checks, row by row and
cell by cell, that they describe the same masked calls:

* the masked VCF is the original with exactly one provenance header line
  added, the same fixed columns, the same FORMAT, and every sample field
  unchanged except that a masked call's GT alleles are `.`;
* the call is masked exactly when the shared policy evaluator says so for
  the *original* sample field, so neither an over-masked nor an under-masked
  cell can pass;
* the matrix token of every cell is the dosage of the original call, or
  `nan` exactly when the call was non-standard or masked;
* a row with masked calls has AC/AN/AF recomputed from its masked GTs, and a
  row without any keeps its INFO byte for byte;
* the variant metadata, sample metadata and genotype accounting report the
  counts this tool re-derives on its own, including the partition
  `total cells = encoded calls + every missing or masked category`.

This tool does not decide what to mask: it imports the same evaluator the
builder used (`gs_genotype_quality.py`), because a second, independent copy
of the masking rule is exactly the drift the Issue asks to avoid. What it
verifies independently is that every artifact agrees with that rule and with
each other. Every mismatch is a hard error; no output is written.

Memory is one row of each input plus per-sample counters, as in the builder.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path

import gs_genotype_quality as quality
from build_gs_panel import (
    FIXED_COLUMN_COUNT,
    MISSING_CELL_TOKEN,
    SAMPLE_METADATA_HEADER_WITH_QUALITY,
    VARIANT_METADATA_HEADER_WITH_QUALITY,
    classify_genotype,
)

OUTPUT_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")

_DOSAGE_METRIC = {
    "-1": "standard_hom_ref_calls",
    "0": "standard_het_calls",
    "1": "standard_hom_alt_calls",
}
_SHAPE_METRIC = {
    "missing": "missing_calls",
    "non_diploid": "non_diploid_calls_treated_as_missing",
    "non_biallelic_index": "non_biallelic_index_calls_treated_as_missing",
}
_MASKED_BY_DOSAGE_METRIC = {
    "-1": "quality_masked_hom_ref_calls",
    "0": "quality_masked_het_calls",
    "1": "quality_masked_hom_alt_calls",
}


class InconsistentQualityMaskError(Exception):
    """The quality-masked artifacts disagree with each other or with the policy."""


@dataclass
class _Tally:
    """Counts this tool derives on its own, keyed by accounting metric name."""

    metrics: dict[str, int] = field(default_factory=dict)
    sample_missing: list[int] = field(default_factory=list)
    sample_non_standard: list[int] = field(default_factory=list)
    sample_masked: list[int] = field(default_factory=list)
    variant_records: int = 0
    rows_with_masked_calls: int = 0

    def add(self, metric: str, amount: int = 1) -> None:
        self.metrics[metric] = self.metrics.get(metric, 0) + amount


def _data_lines(handle):
    """Yield (line_number, line) for non-empty lines, stripped of the newline."""
    for line_number, line in enumerate(handle, start=1):
        line = line.rstrip("\n")
        if line:
            yield line_number, line


def _fail(message: str) -> None:
    raise InconsistentQualityMaskError(message)


def _read_headers(original, masked, policy: quality.GenotypeQualityPolicy):
    """Consume both VCF headers and return the sample names.

    The masked header must be the original header with the policy's own
    provenance line inserted immediately before #CHROM, and nothing else.
    """
    original_meta: list[str] = []
    original_chrom = None
    for _number, line in original:
        if line.startswith("#CHROM"):
            original_chrom = line
            break
        if line.startswith("#"):
            original_meta.append(line)
            continue
        _fail("original VCF has a data row before #CHROM")
    if original_chrom is None:
        _fail("original VCF has no #CHROM header")

    expected_masked = [*original_meta, quality.mask_header_line(policy), original_chrom]
    for index, expected in enumerate(expected_masked):
        try:
            _number, line = next(masked)
        except StopIteration:
            _fail("quality-masked VCF header ends early")
        if line != expected:
            _fail(
                f"quality-masked VCF header line {index + 1} differs from the original "
                "header plus the policy provenance line"
            )
    return tuple(original_chrom.split("\t")[FIXED_COLUMN_COUNT:])


def verify(
    *,
    cohort_id: str,
    policy: quality.GenotypeQualityPolicy,
    gs_pass_vcf: Path,
    quality_masked_vcf: Path,
    matrix: Path,
    variant_metadata: Path,
    sample_metadata: Path,
    genotype_accounting: Path,
) -> _Tally:
    tally = _Tally()
    with ExitStack() as stack:
        original = _data_lines(stack.enter_context(gzip.open(gs_pass_vcf, "rt", encoding="utf-8")))
        masked = _data_lines(
            stack.enter_context(gzip.open(quality_masked_vcf, "rt", encoding="utf-8"))
        )
        matrix_lines = _data_lines(stack.enter_context(gzip.open(matrix, "rt", encoding="utf-8")))
        variant_lines = _data_lines(stack.enter_context(variant_metadata.open(encoding="utf-8")))

        samples = _read_headers(original, masked, policy)
        sample_count = len(samples)
        tally.sample_missing = [0] * sample_count
        tally.sample_non_standard = [0] * sample_count
        tally.sample_masked = [0] * sample_count

        _n, matrix_header = next(matrix_lines, (0, None))
        if matrix_header != "\t".join(["variant_key", *samples]):
            _fail("matrix header does not list the VCF's samples in order")
        _n, variant_header = next(variant_lines, (0, None))
        if variant_header != "\t".join(VARIANT_METADATA_HEADER_WITH_QUALITY):
            _fail("variant metadata header is not the quality-policy header")
        column = {name: index for index, name in enumerate(VARIANT_METADATA_HEADER_WITH_QUALITY)}

        record_index = -1
        for _number, original_line in original:
            if original_line.startswith("#"):
                # The builder ignores a stray header line after #CHROM too.
                continue
            record_index += 1
            masked_row = next(masked, None)
            matrix_row = next(matrix_lines, None)
            variant_row = next(variant_lines, None)
            if masked_row is None or matrix_row is None or variant_row is None:
                _fail(f"record {record_index}: an artifact ends before the original VCF")
            _verify_row(
                tally,
                policy,
                samples,
                record_index,
                original_line,
                masked_row[1],
                matrix_row[1],
                variant_row[1],
                column,
                cohort_id,
            )

        for label, stream in (
            ("quality-masked VCF", masked),
            ("matrix", matrix_lines),
            ("variant metadata", variant_lines),
        ):
            if next(stream, None) is not None:
                _fail(f"{label} has more rows than the original VCF")

    _verify_sample_metadata(tally, samples, sample_metadata)
    _verify_accounting(tally, policy, genotype_accounting)
    return tally


def _verify_row(
    tally: _Tally,
    policy: quality.GenotypeQualityPolicy,
    samples: tuple[str, ...],
    record_index: int,
    original_line: str,
    masked_line: str,
    matrix_line: str,
    variant_line: str,
    column: dict[str, int],
    cohort_id: str,
) -> None:
    fields = original_line.split("\t")
    masked_fields = masked_line.split("\t")
    where = f"record {record_index} ({fields[0]}:{fields[1] if len(fields) > 1 else '?'})"
    width = FIXED_COLUMN_COUNT + len(samples)
    if len(fields) != width or len(masked_fields) != width:
        _fail(f"{where}: VCF row width differs from its header")
    if masked_fields[:7] != fields[:7] or masked_fields[8] != fields[8]:
        _fail(f"{where}: fixed columns or FORMAT differ between the original and masked VCF")

    variant_key = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
    tokens = matrix_line.split("\t")
    if len(tokens) != 1 + len(samples) or tokens[0] != variant_key:
        _fail(f"{where}: matrix row is not this variant or has the wrong width")

    format_keys = fields[8].split(":")
    if "GT" not in format_keys:
        _fail(f"{where}: FORMAT has no GT")
    gt_index = format_keys.index("GT")
    evaluator = quality.RowQualityEvaluator(policy, fields[8])

    row_missing = row_masked = 0
    masked_genotypes: list[str] = []
    for position, (original_field, masked_field) in enumerate(
        zip(fields[FIXED_COLUMN_COUNT:], masked_fields[FIXED_COLUMN_COUNT:])
    ):
        subfields = original_field.split(":")
        if gt_index >= len(subfields):
            _fail(f"{where}: sample '{samples[position]}' has no GT subfield")
        gt = subfields[gt_index]
        cell = classify_genotype(gt)
        tally.add("total_genotype_cells")
        if cell.is_phased:
            tally.add("phased_genotype_count")

        expected_field = original_field
        expected_token = cell.dosage
        if cell.category != "standard":
            tally.add(_SHAPE_METRIC[cell.category])
            tally.sample_missing[position] += 1
            if cell.category != "missing":
                tally.sample_non_standard[position] += 1
            row_missing += 1
        else:
            try:
                verdict = evaluator.evaluate(subfields)
            except quality.GenotypeQualityPolicyError as error:
                _fail(
                    f"{where}: sample '{samples[position]}': a published panel holds a call "
                    f"its own policy rejects: {error}"
                )
            tally.add("quality_evaluated_calls")
            for (key, _), status in zip(policy.thresholds, verdict.statuses):
                tally.add(quality.status_metric(key, status))
            if verdict.masked:
                tally.add("quality_masked_calls")
                tally.add(_MASKED_BY_DOSAGE_METRIC[cell.dosage])
                tally.add(quality.reason_metric(policy, verdict.reasons))
                tally.sample_missing[position] += 1
                tally.sample_masked[position] += 1
                row_missing += 1
                row_masked += 1
                subfields[gt_index] = quality.masked_genotype(gt)
                expected_field = ":".join(subfields)
                expected_token = MISSING_CELL_TOKEN
            else:
                tally.add(
                    "quality_kept_unevaluated_calls"
                    if verdict.unevaluated
                    else "quality_passed_calls"
                )
                tally.add(_DOSAGE_METRIC[cell.dosage])

        if masked_field != expected_field:
            _fail(
                f"{where}: sample '{samples[position]}': masked VCF field differs from the "
                "policy's verdict on the original call"
            )
        if tokens[1 + position] != expected_token:
            _fail(
                f"{where}: sample '{samples[position]}': matrix token "
                f"{tokens[1 + position]!r} differs from expected {expected_token!r}"
            )
        masked_genotypes.append(masked_field.split(":")[gt_index])

    if row_masked:
        tally.rows_with_masked_calls += 1
        ac, an = quality.allele_counts(masked_genotypes)
        if masked_fields[7] != quality.rewrite_allele_info(fields[7], ac, an):
            _fail(f"{where}: masked VCF INFO is not the original with AC/AN/AF recomputed")
        info = dict(entry.split("=", 1) for entry in masked_fields[7].split(";") if "=" in entry)
        if int(info["AC"]) != ac or int(info["AN"]) != an:
            _fail(f"{where}: AC/AN do not match the masked genotypes")
    elif masked_fields[7] != fields[7]:
        _fail(f"{where}: INFO changed on a row with no masked call")

    values = variant_line.split("\t")
    expected_metadata = {
        "cohort_id": cohort_id,
        "variant_index": str(record_index),
        "variant_key": variant_key,
        "missing_genotype_count": str(row_missing),
        "quality_masked_genotype_count": str(row_masked),
    }
    if len(values) != len(column) or any(
        values[column[name]] != value for name, value in expected_metadata.items()
    ):
        _fail(f"{where}: variant metadata row disagrees with the masked state")
    tally.variant_records += 1


def _verify_sample_metadata(tally: _Tally, samples: tuple[str, ...], path: Path) -> None:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not lines or lines[0] != "\t".join(SAMPLE_METADATA_HEADER_WITH_QUALITY):
        _fail("sample metadata header is not the quality-policy header")
    column = {name: index for index, name in enumerate(SAMPLE_METADATA_HEADER_WITH_QUALITY)}
    if len(lines) - 1 != len(samples):
        _fail("sample metadata does not have one row per sample")
    for position, line in enumerate(lines[1:]):
        values = line.split("\t")
        expected = {
            "sample_index": str(position),
            "sample_id": samples[position],
            "missing_genotype_count": str(tally.sample_missing[position]),
            "non_standard_genotype_count": str(tally.sample_non_standard[position]),
            "quality_masked_genotype_count": str(tally.sample_masked[position]),
        }
        if len(values) != len(column) or any(
            values[column[name]] != value for name, value in expected.items()
        ):
            _fail(
                f"sample metadata row {position} ('{samples[position]}') disagrees with the masked state"
            )


def _verify_accounting(tally: _Tally, policy: quality.GenotypeQualityPolicy, path: Path) -> None:
    reported: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        if line:
            parts = line.split("\t")
            reported[parts[1]] = parts[2]

    derived = dict(tally.metrics)
    counts = derived.get
    derived["total_treated_as_missing"] = (
        counts("missing_calls", 0)
        + counts("non_diploid_calls_treated_as_missing", 0)
        + counts("non_biallelic_index_calls_treated_as_missing", 0)
        + counts("quality_masked_calls", 0)
    )
    expected_metrics = [
        "total_genotype_cells",
        *_DOSAGE_METRIC.values(),
        *_SHAPE_METRIC.values(),
        "total_treated_as_missing",
        "phased_genotype_count",
        "quality_evaluated_calls",
        "quality_passed_calls",
        "quality_kept_unevaluated_calls",
        "quality_masked_calls",
        *_MASKED_BY_DOSAGE_METRIC.values(),
        *quality.reason_metrics(policy),
        *quality.status_metrics(policy),
    ]
    for metric in expected_metrics:
        if reported.get(metric) != str(derived.get(metric, 0)):
            _fail(
                f"genotype accounting {metric}={reported.get(metric)!r}, re-derived {derived.get(metric, 0)}"
            )
    if reported.get("genotype_quality_policy_hash") != policy.policy_hash():
        _fail("genotype accounting records a different policy hash")

    # The partitions, stated explicitly rather than left implied by equal counts.
    encoded = sum(derived.get(metric, 0) for metric in _DOSAGE_METRIC.values())
    if derived.get("total_genotype_cells", 0) != encoded + derived["total_treated_as_missing"]:
        _fail("total cells != encoded calls + every missing or masked category")
    evaluated = derived.get("quality_evaluated_calls", 0)
    if evaluated != encoded + derived.get("quality_masked_calls", 0):
        _fail("evaluated calls != encoded calls + masked calls")
    if evaluated != sum(
        derived.get(m, 0)
        for m in ("quality_passed_calls", "quality_kept_unevaluated_calls", "quality_masked_calls")
    ):
        _fail("evaluated calls != passed + kept unevaluated + masked")
    if derived.get("quality_masked_calls", 0) != sum(
        derived.get(m, 0) for m in quality.reason_metrics(policy)
    ):
        _fail("masked calls != the sum of the reason pairs")
    for key, _ in policy.thresholds:
        if evaluated != sum(
            derived.get(quality.status_metric(key, s), 0) for s in quality.FIELD_STATUSES
        ):
            _fail(f"{key} statuses do not partition the evaluated calls")


def build_output_rows(
    cohort_id: str, policy: quality.GenotypeQualityPolicy, tally: _Tally
) -> list[list[str]]:
    return [
        [cohort_id, "genotype_quality_policy_hash", policy.policy_hash()],
        [cohort_id, "verified_variant_records", str(tally.variant_records)],
        [cohort_id, "verified_genotype_cells", str(tally.metrics.get("total_genotype_cells", 0))],
        [cohort_id, "quality_masked_calls", str(tally.metrics.get("quality_masked_calls", 0))],
        [cohort_id, "rows_with_masked_calls", str(tally.rows_with_masked_calls)],
        [cohort_id, "verification_status", "consistent"],
    ]


def build_summary_text(cohort_id: str, policy: quality.GenotypeQualityPolicy, tally: _Tally) -> str:
    lines = [
        "GS panel genotype quality mask verification",
        f"Cohort ID: {cohort_id}",
        f"Policy hash: {policy.policy_hash()}",
        f"Variant records checked cell by cell: {tally.variant_records}",
        f"Genotype cells checked: {tally.metrics.get('total_genotype_cells', 0)}",
        f"Masked calls: {tally.metrics.get('quality_masked_calls', 0)} "
        f"in {tally.rows_with_masked_calls} row(s)",
        (
            "The original PASS VCF, the quality-masked VCF, the matrix, both metadata files "
            "and the genotype accounting all describe the same masked calls, each call was "
            "masked exactly when the recorded policy masks it, and AC/AN/AF were recomputed "
            "from the masked genotypes on every changed row. This file would not exist otherwise."
        ),
    ]
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--genotype-quality-policy", required=True, type=Path)
    parser.add_argument("--gs-pass-vcf", required=True, type=Path)
    parser.add_argument("--quality-masked-vcf", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--variant-metadata", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--genotype-accounting", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        policy = quality.policy_from_document(
            json.loads(args.genotype_quality_policy.read_text(encoding="utf-8"))
        )
        if policy is None:
            raise quality.InvalidPolicyError("the policy document is the disabled policy")
        tally = verify(
            cohort_id=args.cohort_id,
            policy=policy,
            gs_pass_vcf=args.gs_pass_vcf,
            quality_masked_vcf=args.quality_masked_vcf,
            matrix=args.matrix,
            variant_metadata=args.variant_metadata,
            sample_metadata=args.sample_metadata,
            genotype_accounting=args.genotype_accounting,
        )
    except (
        OSError,
        ValueError,
        KeyError,
        quality.InvalidPolicyError,
        InconsistentQualityMaskError,
    ) as error:
        print(f"verify_gs_genotype_quality_mask.py: error: {error}", file=sys.stderr)
        return 1

    lines = [
        "\t".join(OUTPUT_HEADER),
        *("\t".join(row) for row in build_output_rows(args.cohort_id, policy, tally)),
    ]
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    args.summary_output.write_text(
        build_summary_text(args.cohort_id, policy, tally), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
