#!/usr/bin/env python3
"""Summarize hard-filter annotation distributions and threshold sensitivity.

This is an evidence-generation CLI, not a filter optimizer. It evaluates a
small, versioned set of exploratory threshold scenarios against annotations in
an existing GATK VariantFiltration VCF. No truth set is available, so its output
describes variant-set sensitivity only and must not be interpreted as accuracy.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from summarize_filter_qc import (
    ANNOTATION_NAMES,
    FILTER_TAG_BY_VARIANT_TYPE,
    MalformedVcfError,
    VcfRecord,
    _annotation_value,
    _filter_tags,
    _is_present,
    parse_filtered_vcf,
)


DEFAULT_SCENARIO_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "conf"
    / "hard_filter_sensitivity_scenarios.json"
)
NOT_APPLICABLE = "NA"

FILTERED_ANNOTATIONS = {
    variant_type: tuple(
        annotation for annotation, tag in tags.items() if tag is not None
    )
    for variant_type, tags in FILTER_TAG_BY_VARIANT_TYPE.items()
}

OPERATORS: dict[str, str] = {
    "QD": "<",
    "QUAL": "<",
    "SOR": ">",
    "FS": ">",
    "MQ": "<",
    "MQRankSum": "<",
    "ReadPosRankSum": "<",
}

HISTOGRAM_EDGES: dict[str, tuple[float, ...]] = {
    "QD": (0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0),
    "QUAL": (0.0, 10.0, 20.0, 30.0, 50.0, 100.0, 1000.0),
    "SOR": (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0),
    "FS": (0.0, 10.0, 30.0, 60.0, 100.0, 200.0, 300.0, 500.0),
    "MQ": (0.0, 20.0, 30.0, 40.0, 50.0, 60.0),
    "MQRankSum": (-20.0, -15.0, -12.5, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0),
    "ReadPosRankSum": (-25.0, -20.0, -10.0, -8.0, -5.0, 0.0, 5.0, 10.0, 20.0),
}

DISTRIBUTION_HEADER = (
    "cohort_id",
    "variant_type",
    "annotation",
    "row_type",
    "bin_lower_inclusive",
    "bin_upper_exclusive",
    "total_records",
    "present_records",
    "missing_records",
    "evaluable_rate",
    "value_count",
    "rate_among_present",
    "min_value",
    "max_value",
    "mean_value",
)

SENSITIVITY_HEADER = (
    "cohort_id",
    "variant_type",
    "scenario",
    "annotation",
    "operator",
    "threshold",
    "total_records",
    "present_records",
    "missing_records",
    "hit_records",
    "hit_rate_among_present",
    "hit_rate_among_total",
    "observed_filter_tag_records",
    "predicted_minus_observed",
)


class ScenarioConfigError(Exception):
    """Raised when a scenario file does not match the fixed analysis contract."""


@dataclass
class DistributionAccumulator:
    """Bounded streaming statistics for one annotation."""

    edges: tuple[float, ...]
    counts: list[int]
    present: int = 0
    minimum: float | None = None
    maximum: float | None = None
    total: float = 0.0

    @classmethod
    def create(cls, annotation: str) -> DistributionAccumulator:
        edges = HISTOGRAM_EDGES[annotation]
        return cls(edges=edges, counts=[0] * (len(edges) + 1))

    def add(self, value: float) -> None:
        self.present += 1
        self.total += value
        self.minimum = value if self.minimum is None else min(self.minimum, value)
        self.maximum = value if self.maximum is None else max(self.maximum, value)
        self.counts[bisect.bisect_right(self.edges, value)] += 1


@dataclass(frozen=True)
class AnalysisResult:
    """All bounded counters produced by one VCF scan."""

    total_records: int
    distributions: dict[str, DistributionAccumulator]
    scenario_hits: dict[str, dict[str, int]]
    scenario_any_hits: dict[str, int]
    observed_tag_hits: dict[str, int]
    observed_any_hits: int


def load_scenarios(path: Path) -> dict[str, dict[str, dict[str, float]]]:
    """Load and validate the versioned scenario configuration."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScenarioConfigError(f"cannot read {path}: {error}") from error

    if payload.get("schema_version") != 1:
        raise ScenarioConfigError(f"{path}: schema_version must be 1")

    scenarios = payload.get("scenarios")
    if not isinstance(scenarios, dict) or not scenarios:
        raise ScenarioConfigError(f"{path}: scenarios must be a non-empty object")
    if "current" not in scenarios:
        raise ScenarioConfigError(f"{path}: scenarios.current is required")

    validated: dict[str, dict[str, dict[str, float]]] = {}
    for scenario_name, by_variant_type in scenarios.items():
        if not isinstance(scenario_name, str) or not scenario_name:
            raise ScenarioConfigError(
                f"{path}: scenario names must be non-empty strings"
            )
        if not isinstance(by_variant_type, dict):
            raise ScenarioConfigError(
                f"{path}: scenario {scenario_name} must be an object"
            )

        validated[scenario_name] = {}
        for variant_type, expected_annotations in FILTERED_ANNOTATIONS.items():
            thresholds = by_variant_type.get(variant_type)
            if not isinstance(thresholds, dict):
                raise ScenarioConfigError(
                    f"{path}: scenario {scenario_name}.{variant_type} must be an object"
                )
            if set(thresholds) != set(expected_annotations):
                raise ScenarioConfigError(
                    f"{path}: scenario {scenario_name}.{variant_type} must define exactly "
                    f"{', '.join(expected_annotations)}"
                )

            validated[scenario_name][variant_type] = {}
            for annotation, raw_threshold in thresholds.items():
                if isinstance(raw_threshold, bool) or not isinstance(
                    raw_threshold, (int, float)
                ):
                    raise ScenarioConfigError(
                        f"{path}: {scenario_name}.{variant_type}.{annotation} must be numeric"
                    )
                threshold = float(raw_threshold)
                if not math.isfinite(threshold):
                    raise ScenarioConfigError(
                        f"{path}: {scenario_name}.{variant_type}.{annotation} must be finite"
                    )
                validated[scenario_name][variant_type][annotation] = threshold

    return validated


def _parse_numeric_annotation(record: VcfRecord, annotation: str) -> float | None:
    raw_value = _annotation_value(record, annotation)
    if not _is_present(raw_value):
        return None

    try:
        value = float(raw_value)
    except (TypeError, ValueError) as error:
        raise MalformedVcfError(
            f"{record.chrom}:{record.pos}: annotation {annotation} is not numeric: {raw_value!r}"
        ) from error
    if not math.isfinite(value):
        raise MalformedVcfError(
            f"{record.chrom}:{record.pos}: annotation {annotation} is not finite: {raw_value!r}"
        )
    return value


def _hits(value: float, operator: str, threshold: float) -> bool:
    if operator == "<":
        return value < threshold
    if operator == ">":
        return value > threshold
    raise AssertionError(f"unsupported operator: {operator}")


def analyze_records(
    records: Iterable[VcfRecord],
    variant_type: str,
    scenarios: dict[str, dict[str, dict[str, float]]],
) -> AnalysisResult:
    """Scan records once and retain only fixed-size distribution counters."""
    filtered_annotations = FILTERED_ANNOTATIONS[variant_type]
    distributions = {
        annotation: DistributionAccumulator.create(annotation)
        for annotation in ANNOTATION_NAMES
    }
    scenario_hits = {
        scenario: {annotation: 0 for annotation in filtered_annotations}
        for scenario in scenarios
    }
    scenario_any_hits = {scenario: 0 for scenario in scenarios}
    observed_tag_hits = {annotation: 0 for annotation in filtered_annotations}
    expected_tags = {
        FILTER_TAG_BY_VARIANT_TYPE[variant_type][annotation]
        for annotation in filtered_annotations
    }
    total_records = 0
    observed_any_hits = 0

    for record in records:
        total_records += 1
        values = {
            annotation: _parse_numeric_annotation(record, annotation)
            for annotation in ANNOTATION_NAMES
        }
        tags = set(_filter_tags(record))
        if tags & expected_tags:
            observed_any_hits += 1

        for annotation, value in values.items():
            if value is not None:
                distributions[annotation].add(value)
            if (
                annotation in observed_tag_hits
                and FILTER_TAG_BY_VARIANT_TYPE[variant_type][annotation] in tags
            ):
                observed_tag_hits[annotation] += 1

        for scenario, by_variant_type in scenarios.items():
            any_hit = False
            for annotation, threshold in by_variant_type[variant_type].items():
                value = values[annotation]
                if value is not None and _hits(value, OPERATORS[annotation], threshold):
                    scenario_hits[scenario][annotation] += 1
                    any_hit = True
            if any_hit:
                scenario_any_hits[scenario] += 1

    return AnalysisResult(
        total_records=total_records,
        distributions=distributions,
        scenario_hits=scenario_hits,
        scenario_any_hits=scenario_any_hits,
        observed_tag_hits=observed_tag_hits,
        observed_any_hits=observed_any_hits,
    )


def _format_rate(numerator: int, denominator: int) -> str:
    return NOT_APPLICABLE if denominator == 0 else f"{numerator / denominator:.6f}"


def _format_number(value: float | None) -> str:
    return NOT_APPLICABLE if value is None else format(value, ".12g")


def build_distribution_rows(
    cohort_id: str, variant_type: str, result: AnalysisResult
) -> list[list[str]]:
    rows: list[list[str]] = []
    for annotation in ANNOTATION_NAMES:
        distribution = result.distributions[annotation]
        missing = result.total_records - distribution.present
        common = [
            cohort_id,
            variant_type,
            annotation,
        ]
        rows.append(
            common
            + [
                "summary",
                NOT_APPLICABLE,
                NOT_APPLICABLE,
                str(result.total_records),
                str(distribution.present),
                str(missing),
                _format_rate(distribution.present, result.total_records),
                str(distribution.present),
                "1.000000" if distribution.present else NOT_APPLICABLE,
                _format_number(distribution.minimum),
                _format_number(distribution.maximum),
                _format_number(
                    distribution.total / distribution.present
                    if distribution.present
                    else None
                ),
            ]
        )

        bounds = (None, *distribution.edges, None)
        for index, count in enumerate(distribution.counts):
            rows.append(
                common
                + [
                    "histogram",
                    "-inf" if bounds[index] is None else _format_number(bounds[index]),
                    "+inf"
                    if bounds[index + 1] is None
                    else _format_number(bounds[index + 1]),
                    str(result.total_records),
                    str(distribution.present),
                    str(missing),
                    _format_rate(distribution.present, result.total_records),
                    str(count),
                    _format_rate(count, distribution.present),
                    NOT_APPLICABLE,
                    NOT_APPLICABLE,
                    NOT_APPLICABLE,
                ]
            )
    return rows


def build_sensitivity_rows(
    cohort_id: str,
    variant_type: str,
    scenarios: dict[str, dict[str, dict[str, float]]],
    result: AnalysisResult,
) -> list[list[str]]:
    rows: list[list[str]] = []
    for scenario, by_variant_type in scenarios.items():
        for annotation in FILTERED_ANNOTATIONS[variant_type]:
            distribution = result.distributions[annotation]
            missing = result.total_records - distribution.present
            hit_records = result.scenario_hits[scenario][annotation]
            observed = (
                result.observed_tag_hits[annotation] if scenario == "current" else None
            )
            rows.append(
                [
                    cohort_id,
                    variant_type,
                    scenario,
                    annotation,
                    OPERATORS[annotation],
                    _format_number(by_variant_type[variant_type][annotation]),
                    str(result.total_records),
                    str(distribution.present),
                    str(missing),
                    str(hit_records),
                    _format_rate(hit_records, distribution.present),
                    _format_rate(hit_records, result.total_records),
                    str(observed) if observed is not None else NOT_APPLICABLE,
                    str(hit_records - observed)
                    if observed is not None
                    else NOT_APPLICABLE,
                ]
            )

        any_hits = result.scenario_any_hits[scenario]
        observed_any = result.observed_any_hits if scenario == "current" else None
        rows.append(
            [
                cohort_id,
                variant_type,
                scenario,
                "ANY_FILTER",
                "any",
                NOT_APPLICABLE,
                str(result.total_records),
                str(result.total_records),
                "0",
                str(any_hits),
                _format_rate(any_hits, result.total_records),
                _format_rate(any_hits, result.total_records),
                str(observed_any) if observed_any is not None else NOT_APPLICABLE,
                str(any_hits - observed_any)
                if observed_any is not None
                else NOT_APPLICABLE,
            ]
        )
    return rows


def _write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_summary(
    args: argparse.Namespace,
    result: AnalysisResult,
) -> str:
    return "\n".join(
        [
            "Hard-filter threshold sensitivity summary",
            f"Cohort ID: {args.cohort_id}",
            f"Variant type: {args.variant_type}",
            f"Input file: {args.filtered_vcf.name}",
            f"Input SHA256: {_sha256(args.filtered_vcf)}",
            f"Scenario config: {args.scenario_config.name}",
            f"Scenario config SHA256: {_sha256(args.scenario_config)}",
            f"Total records: {result.total_records}",
            "",
            "Interpretation boundary:",
            "  These counts measure variant-set sensitivity to exploratory thresholds.",
            "  Without a truth set they do not estimate accuracy, precision, recall,",
            "  biological validity, or an optimal threshold.",
            "  QUAL can depend strongly on cohort/sample size and must not be read as",
            "  a cohort-independent quality score.",
            "",
        ]
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--filtered-vcf", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument(
        "--variant-type", required=True, choices=sorted(FILTERED_ANNOTATIONS)
    )
    parser.add_argument("--scenario-config", type=Path, default=DEFAULT_SCENARIO_CONFIG)
    parser.add_argument("--distribution-output", required=True, type=Path)
    parser.add_argument("--sensitivity-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        scenarios = load_scenarios(args.scenario_config)
        result = analyze_records(
            parse_filtered_vcf(args.filtered_vcf),
            args.variant_type,
            scenarios,
        )
        distribution_rows = build_distribution_rows(
            args.cohort_id, args.variant_type, result
        )
        sensitivity_rows = build_sensitivity_rows(
            args.cohort_id,
            args.variant_type,
            scenarios,
            result,
        )
        summary = build_summary(args, result)
    except (OSError, MalformedVcfError, ScenarioConfigError) as error:
        print(f"analyze_hard_filter_sensitivity.py: error: {error}", file=sys.stderr)
        return 1

    _write_tsv(args.distribution_output, DISTRIBUTION_HEADER, distribution_rows)
    _write_tsv(args.sensitivity_output, SENSITIVITY_HEADER, sensitivity_rows)
    args.summary_output.write_text(summary, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
