"""Evidence classes, evaluation records and the delivery-support rules (Issue #65).

The point of this module is to make one mistake impossible to write down: a
concordance, a stability or a descriptive number presented as accuracy.

Each evidence class names the metrics it may carry. Precision, recall, F1 and
the TP/FP/FN vocabulary exist only for `independent_truth` and for the
`synthetic_truth_fixture` that tests the harness; every other class is refused
if it tries to carry them. Each class also bounds the delivery status it can
support: only independent truth can make a scope `supported`, while
replicate, cross-platform, caller, downsampling and descriptive evidence can at
most make it `supported_with_caveat`. `not_evaluated` is a first-class record
with a required reason and no metrics, not an absent row.
"""

from __future__ import annotations

import hashlib
import json

EVALUATION_SCHEMA_VERSION = 1

TRUTH_METRICS = frozenset(
    {
        "true_positive",
        "false_positive",
        "false_negative",
        "query_nocall_at_truth_variant",
        "truth_nocall_at_query_variant",
        "precision",
        "recall",
        "f1",
        "genotype_concordant",
        "genotype_discordant",
        "genotype_concordance",
        "matched_after_normalization",
    }
)
CONCORDANCE_METRICS = frozenset(
    {
        "shared_variants",
        "only_in_query",
        "only_in_comparator",
        "query_nocall_at_comparator_variant",
        "comparator_nocall_at_query_variant",
        "site_agreement_rate",
        "genotype_concordant",
        "genotype_discordant",
        "genotype_concordance",
        "matched_after_normalization",
    }
)
STABILITY_METRICS = frozenset(
    {
        "shared_variants",
        "only_in_full_depth",
        "only_in_downsampled",
        "downsampled_nocall_at_full_depth_variant",
        "full_depth_nocall_at_downsampled_variant",
        "call_retention",
        "new_call_fraction",
        "genotype_concordant",
        "genotype_discordant",
        "genotype_concordance",
        "matched_after_normalization",
    }
)
DESCRIPTIVE_METRICS = frozenset(
    {
        "bases",
        "fraction_of_universe",
        "variant_records",
        "variant_records_per_mb",
        "snp_records",
        "indel_records",
        "genotype_cells",
        "missing_genotype_cells",
        "missing_genotype_fraction",
        "masked_genotype_cells",
        "masked_genotype_fraction",
        "pass_records",
        "failed_filter_records",
        "unfiltered_records",
        "mean_depth",
        "median_depth",
        "non_reference_calls",
        "heterozygous_calls",
        "homozygous_alt_calls",
        "heterozygous_fraction_of_non_reference_calls",
        "non_reference_calls_per_mb",
        "callable_bases",
        "samples",
    }
)

EVIDENCE_CLASSES: dict[str, dict[str, object]] = {
    "independent_truth": {
        "metrics": TRUTH_METRICS,
        "max_status": "supported",
        "claim": "accuracy (precision/recall) against the named independent truth, within its region",
    },
    "technical_replicate_concordance": {
        "metrics": CONCORDANCE_METRICS,
        "max_status": "supported_with_caveat",
        "claim": "reproducibility between replicates of one sample; not accuracy",
    },
    "cross_platform_concordance": {
        "metrics": CONCORDANCE_METRICS | DESCRIPTIVE_METRICS,
        "max_status": "supported_with_caveat",
        "claim": "agreement with another platform or method on the same sample; not accuracy",
    },
    "caller_concordance": {
        "metrics": CONCORDANCE_METRICS,
        "max_status": "supported_with_caveat",
        "claim": "agreement between callers on the same reads; not accuracy",
    },
    "downsampling_stability": {
        "metrics": STABILITY_METRICS,
        "max_status": "supported_with_caveat",
        "claim": "how calls change as depth falls, relative to the full-depth calls; not accuracy",
    },
    "reference_sample_self_consistency": {
        "metrics": DESCRIPTIVE_METRICS,
        "max_status": "supported_with_caveat",
        "claim": (
            "where short reads of the reference's own BioSample disagree with the reference "
            "consensus; not truth, and not callset-to-callset cross-platform concordance"
        ),
    },
    "descriptive_stratification": {
        "metrics": DESCRIPTIVE_METRICS,
        "max_status": "supported_with_caveat",
        "claim": "counts and rates by region, depth and variant type; no quality claim",
    },
    "synthetic_truth_fixture": {
        "metrics": TRUTH_METRICS,
        "max_status": "validated_test_harness",
        "claim": "the evaluation harness computes truth metrics correctly on constructed data",
    },
}

# Classes whose whole meaning is a comparison of two callsets. An evaluated
# record of one of these must name both datasets: a "concordance" with nothing
# to be concordant with is a description, and belongs to a descriptive class.
COMPARISON_CLASSES = frozenset(
    {
        "independent_truth",
        "technical_replicate_concordance",
        "cross_platform_concordance",
        "caller_concordance",
        "downsampling_stability",
        "synthetic_truth_fixture",
    }
)

STATUSES = (
    "supported",
    "supported_with_caveat",
    "not_evaluated",
    "unsupported",
    "validated_test_harness",
)
_STATUS_RANK = {"unsupported": 0, "supported_with_caveat": 1, "supported": 2}

REQUIRED_FIELDS = (
    "evaluation_schema_version",
    "evaluation_id",
    "reference",
    "region_set",
    "region_definition_hash",
    "evidence_class",
    "query_dataset",
    "comparator_dataset",
    "comparison_unit",
    "eligible_denominator",
    "excluded_denominator",
    "exclusion_reasons",
    "not_evaluated_reason",
    "metrics",
    "limitations",
)


class InvalidEvaluationError(ValueError):
    """An evaluation record breaks the evidence contract."""


def validate_evaluation(record: dict[str, object]) -> None:
    missing = [name for name in REQUIRED_FIELDS if name not in record]
    if missing:
        raise InvalidEvaluationError(f"evaluation lacks {missing}")
    if record["evaluation_schema_version"] != EVALUATION_SCHEMA_VERSION:
        raise InvalidEvaluationError("unknown evaluation_schema_version")
    cls = record["evidence_class"]
    if cls not in EVIDENCE_CLASSES:
        raise InvalidEvaluationError(f"unknown evidence_class {cls!r}")
    if not isinstance(record["limitations"], list) or not record["limitations"]:
        raise InvalidEvaluationError("limitations must be a non-empty list")

    if record["not_evaluated_reason"] is not None:
        if (
            not isinstance(record["not_evaluated_reason"], str)
            or not record["not_evaluated_reason"].strip()
        ):
            raise InvalidEvaluationError("not_evaluated_reason must be a non-empty string or null")
        if record["metrics"]:
            raise InvalidEvaluationError("a not-evaluated record cannot carry metrics")
        return

    if cls in COMPARISON_CLASSES:
        for side in ("query_dataset", "comparator_dataset"):
            if not isinstance(record[side], dict) or not record[side]:
                raise InvalidEvaluationError(
                    f"{cls} compares two callsets, so an evaluated record must name its {side}"
                )
    for name in ("eligible_denominator", "excluded_denominator"):
        value = record[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise InvalidEvaluationError(f"{name} must be a non-negative integer")
    reasons = record["exclusion_reasons"]
    if not isinstance(reasons, dict) or sum(reasons.values()) != record["excluded_denominator"]:
        raise InvalidEvaluationError("exclusion_reasons must sum to excluded_denominator")
    allowed = EVIDENCE_CLASSES[cls]["metrics"]
    for stratum, metrics in _metric_groups(record["metrics"]):
        illegal = sorted(set(metrics) - allowed)
        if illegal:
            raise InvalidEvaluationError(
                f"{cls} may not report {illegal} (stratum {stratum!r}); "
                "accuracy vocabulary is reserved for independent truth"
            )
    if not isinstance(record["region_definition_hash"], str) or not record[
        "region_definition_hash"
    ].startswith("sha256:"):
        raise InvalidEvaluationError("region_definition_hash must be a sha256: digest")


def _metric_groups(metrics: object):
    if not isinstance(metrics, dict):
        raise InvalidEvaluationError("metrics must be an object")
    if "by_stratum" in metrics:
        yield "all", {k: v for k, v in metrics.items() if k != "by_stratum"}
        for stratum, values in metrics["by_stratum"].items():
            yield stratum, values
    else:
        yield "all", metrics


def canonical_hash(document: object) -> str:
    text = json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def delivery_row(
    *,
    scope: str,
    evidence_class: str | None,
    status: str,
    evidence_refs: list[str],
    rationale: str,
) -> dict[str, object]:
    """One delivery-support row, refusing a status stronger than its evidence allows."""
    if status not in STATUSES:
        raise InvalidEvaluationError(f"unknown status {status!r}")
    if status == "not_evaluated":
        if evidence_refs:
            raise InvalidEvaluationError(f"{scope}: a not-evaluated scope cites no evidence")
        claim = "no claim; not evaluated"
    else:
        if evidence_class not in EVIDENCE_CLASSES:
            raise InvalidEvaluationError(f"{scope}: status {status} needs a known evidence class")
        max_status = EVIDENCE_CLASSES[evidence_class]["max_status"]
        if status == "validated_test_harness" or max_status == "validated_test_harness":
            if status != max_status:
                raise InvalidEvaluationError(
                    f"{scope}: synthetic fixtures only validate the harness, and only they do"
                )
        elif _STATUS_RANK[status] > _STATUS_RANK[max_status]:
            raise InvalidEvaluationError(
                f"{scope}: {evidence_class} cannot support status {status} (at most {max_status})"
            )
        if not evidence_refs:
            raise InvalidEvaluationError(f"{scope}: status {status} must cite evidence")
        claim = (
            EVIDENCE_CLASSES[evidence_class]["claim"]
            if status != "unsupported"
            else ("no delivery claim; evidence shows this scope is not supported")
        )
    return {
        "scope": scope,
        "evidence_class": evidence_class,
        "status": status,
        "claim_allowed": claim,
        "evidence": evidence_refs,
        "rationale": rationale,
    }
