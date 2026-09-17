"""Compare two callsets for one sample each, with explicit accounting (Issue #65).

The same engine serves every comparison-shaped evidence class; only the names
of the resulting metrics change, and `evidence_model.validate_evaluation`
refuses a class that tries to borrow another's vocabulary. `query` is the
callset being characterized; `comparator` is the independent truth, the other
replicate, the other platform or caller, or the full-depth calls.

Comparison unit
    One normalized ALT allele, keyed by `(contig, pos, ref, alt)` after
    multi-allelic splitting and left-alignment (see `variants.py`). A sample is
    *positive* for an allele when every called allele is known and at least
    one is that ALT (its dosage is how many); *no-call* when any allele of its
    GT is missing; *reference* otherwise.

Accounting, for units with at least one positive side inside the region
    shared          both positive (genotype concordant when dosages match)
    only_in_query   query positive, comparator reference or absent
    only_in_comp.   comparator positive, query reference or absent
    query_nocall    comparator positive, query no-call
    comp_nocall     query positive, comparator no-call

    Units where neither side is positive are not counted. True negatives are
    not defined: a VCF does not enumerate the confident homozygous-reference
    positions a TN would need, and counting every reference base would
    swamp every rate with a denominator the evaluation never examined.

Exclusions, each counted by reason and never folded into a metric
    coordinate_mismatch, ref_mismatch, symbolic_allele, outside_region and
    region_boundary (span only partly inside).

    An exclusion is counted per *comparison unit*, not per side. The same
    unusable allele in both callsets -- an `ALT=*` both call, a REF both get
    wrong -- is one excluded unit, counted once. A unit takes the first reason
    in `EXCLUSION_PRECEDENCE` that applies to it on either side, so the reasons
    always sum to `excluded_denominator`. Which side each excluded unit came
    from is reported under `definitions`, where it cannot be mistaken for an
    additive count.

    Only a unit that some side is positive or no-call for can be excluded: an
    allele every side calls homozygous reference is not a comparison unit.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from evidence_model import EVALUATION_SCHEMA_VERSION, validate_evaluation
from regions import RegionSet, Stratum
from variants import Reference, classify_record, is_symbolic_allele, normalize_allele, read_vcf

POSITIVE, NOCALL, REFERENCE = "positive", "nocall", "reference"
_PRECEDENCE = {POSITIVE: 2, NOCALL: 1, REFERENCE: 0}

# A unit that more than one reason could exclude is counted under the first of
# these: coordinates that are wrong leave nothing to check REF against, a wrong
# REF leaves nothing to normalize, and an allele that cannot be normalized
# cannot be placed in the region.
EXCLUSION_PRECEDENCE = (
    "coordinate_mismatch",
    "ref_mismatch",
    "symbolic_allele",
    "outside_region",
    "region_boundary",
)
_EXCLUSION_RANK = {reason: rank for rank, reason in enumerate(EXCLUSION_PRECEDENCE)}

METRIC_NAMES = {
    "truth": {
        "shared": "true_positive",
        "only_query": "false_positive",
        "only_comparator": "false_negative",
        "query_nocall": "query_nocall_at_truth_variant",
        "comparator_nocall": "truth_nocall_at_query_variant",
    },
    "concordance": {
        "shared": "shared_variants",
        "only_query": "only_in_query",
        "only_comparator": "only_in_comparator",
        "query_nocall": "query_nocall_at_comparator_variant",
        "comparator_nocall": "comparator_nocall_at_query_variant",
    },
    "stability": {
        "shared": "shared_variants",
        "only_query": "only_in_downsampled",
        "only_comparator": "only_in_full_depth",
        "query_nocall": "downsampled_nocall_at_full_depth_variant",
        "comparator_nocall": "full_depth_nocall_at_downsampled_variant",
    },
}
CLASS_FAMILY = {
    "independent_truth": "truth",
    "synthetic_truth_fixture": "truth",
    "technical_replicate_concordance": "concordance",
    "cross_platform_concordance": "concordance",
    "caller_concordance": "concordance",
    "downsampling_stability": "stability",
}


@dataclass
class _Side:
    states: dict[tuple, tuple[str, int, bool]] = field(default_factory=dict)
    spans: dict[tuple, int] = field(default_factory=dict)
    types: dict[tuple, str] = field(default_factory=dict)
    # (contig, pos, ref, alt) exactly as the record wrote it -> exclusion reason
    invalid: dict[tuple, str] = field(default_factory=dict)
    duplicate_keys: int = 0
    records: int = 0


def _load_side(path, sample: str, reference: Reference) -> _Side:
    side = _Side()
    for _header, record in read_vcf(path, [sample]):
        side.records += 1
        gt = record.genotypes[sample]
        if not record.alts:
            continue
        site_reason = classify_record(reference, record.contig, record.pos, record.ref)
        for alt_index, alt in enumerate(record.alts, start=1):
            if None in gt:
                state, dosage = NOCALL, 0
            else:
                dosage = sum(1 for allele in gt if allele == alt_index)
                state = POSITIVE if dosage else REFERENCE
            reason = site_reason or ("symbolic_allele" if is_symbolic_allele(alt) else None)
            if reason is not None:
                if state != REFERENCE:
                    unit = (record.contig, record.pos, record.ref, alt)
                    seen = side.invalid.get(unit)
                    if seen is None or _EXCLUSION_RANK[reason] < _EXCLUSION_RANK[seen]:
                        side.invalid[unit] = reason
                continue
            allele = normalize_allele(
                reference, record.contig, record.pos, record.ref, alt, alt_index
            )
            key = allele.key
            if key in side.states:
                side.duplicate_keys += 1
                if _PRECEDENCE[state] <= _PRECEDENCE[side.states[key][0]]:
                    continue
            side.states[key] = (state, dosage, allele.changed_by_normalization)
            side.spans[key] = allele.span
            side.types[key] = allele.variant_type
    return side


def _rate(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else round(numerator / denominator, 6)


def _metrics(counts: dict[str, int], family: str) -> dict[str, object]:
    names = METRIC_NAMES[family]
    out: dict[str, object] = {names[key]: counts[key] for key in names}
    shared, only_q, only_c = counts["shared"], counts["only_query"], counts["only_comparator"]
    q_nocall, c_nocall = counts["query_nocall"], counts["comparator_nocall"]
    if family == "truth":
        precision = _rate(shared, shared + only_q)
        recall = _rate(shared, shared + only_c + q_nocall)
        out["precision"] = precision
        out["recall"] = recall
        out["f1"] = (
            None
            if not precision or not recall
            else round(2 * precision * recall / (precision + recall), 6)
        )
    elif family == "concordance":
        out["site_agreement_rate"] = _rate(shared, shared + only_q + only_c + q_nocall + c_nocall)
    else:
        out["call_retention"] = _rate(shared, shared + only_c + q_nocall)
        out["new_call_fraction"] = _rate(only_q, shared + only_q)
    out["genotype_concordant"] = counts["gt_concordant"]
    out["genotype_discordant"] = counts["gt_discordant"]
    out["genotype_concordance"] = _rate(counts["gt_concordant"], shared)
    out["matched_after_normalization"] = counts["normalized_match"]
    return out


def _empty_counts() -> dict[str, int]:
    return dict.fromkeys(
        (
            "shared",
            "only_query",
            "only_comparator",
            "query_nocall",
            "comparator_nocall",
            "gt_concordant",
            "gt_discordant",
            "normalized_match",
        ),
        0,
    )


def evaluate(
    *,
    evaluation_id: str,
    evidence_class: str,
    reference: Reference,
    reference_description: dict[str, object],
    region: RegionSet,
    strata: list[Stratum],
    query: dict[str, object],
    comparator: dict[str, object],
    limitations: list[str],
) -> dict[str, object]:
    """Compare `query` against `comparator` and return a validated evaluation record.

    `query` and `comparator` are `{"path", "sample", "dataset": {...}}`.
    """
    if evidence_class not in CLASS_FAMILY:
        raise ValueError(f"{evidence_class} is not a comparison-shaped evidence class")
    family = CLASS_FAMILY[evidence_class]
    reference_sha = str(reference_description["sha256"])
    q = _load_side(query["path"], str(query["sample"]), reference)
    c = _load_side(comparator["path"], str(comparator["sample"]), reference)

    exclusions: dict[str, int] = {}
    invalid_units: dict[tuple, str] = {}
    invalid_sides: dict[tuple, set[str]] = {}
    for side_name, side in (("query", q), ("comparator", c)):
        for unit, reason in side.invalid.items():
            seen = invalid_units.get(unit)
            if seen is None or _EXCLUSION_RANK[reason] < _EXCLUSION_RANK[seen]:
                invalid_units[unit] = reason
            invalid_sides.setdefault(unit, set()).add(side_name)
    by_side: dict[str, dict[str, int]] = {}
    for unit, reason in invalid_units.items():
        exclusions[reason] = exclusions.get(reason, 0) + 1
        sides = invalid_sides[unit]
        where = "both" if len(sides) == 2 else f"{next(iter(sides))}_only"
        by_side.setdefault(reason, {"query_only": 0, "comparator_only": 0, "both": 0})[where] += 1

    totals = _empty_counts()
    by_stratum = {s.name: _empty_counts() for s in strata}
    by_stratum.update(
        {
            name: _empty_counts()
            for name in (
                "snp",
                "indel",
                "het",
                "hom_alt",
                "snp:het",
                "snp:hom_alt",
                "indel:het",
                "indel:hom_alt",
            )
        }
    )
    boundary_by_stratum = {s.name: 0 for s in strata}
    eligible = 0

    for key in sorted(set(q.states) | set(c.states)):
        q_state = q.states.get(key, (REFERENCE, 0, False))
        c_state = c.states.get(key, (REFERENCE, 0, False))
        if POSITIVE not in (q_state[0], c_state[0]):
            continue
        span = q.spans.get(key) or c.spans[key]
        placement = region.contains_span(key[0], key[1], span)
        if placement != "inside":
            reason = "outside_region" if placement == "outside" else "region_boundary"
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        eligible += 1
        if q_state[0] == POSITIVE and c_state[0] == POSITIVE:
            bucket = "shared"
        elif c_state[0] == POSITIVE:
            bucket = "query_nocall" if q_state[0] == NOCALL else "only_comparator"
        else:
            bucket = "comparator_nocall" if c_state[0] == NOCALL else "only_query"
        variant_type = q.types.get(key) or c.types[key]
        dosage = c_state[1] if c_state[0] == POSITIVE else q_state[1]
        genotype_class = "hom_alt" if dosage >= 2 else "het"
        targets = [
            totals,
            by_stratum[variant_type],
            by_stratum[genotype_class],
            by_stratum[f"{variant_type}:{genotype_class}"],
        ]
        for stratum in strata:
            placement = stratum.region.contains_span(key[0], key[1], span)
            if placement == "inside":
                targets.append(by_stratum[stratum.name])
                targets.append(
                    by_stratum.setdefault(f"{stratum.name}:{variant_type}", _empty_counts())
                )
                targets.append(
                    by_stratum.setdefault(f"{stratum.name}:{genotype_class}", _empty_counts())
                )
            elif placement == "boundary":
                boundary_by_stratum[stratum.name] += 1
        for counts in targets:
            counts[bucket] += 1
            if bucket == "shared":
                counts["gt_concordant" if q_state[1] == c_state[1] else "gt_discordant"] += 1
                if q_state[2] or c_state[2]:
                    counts["normalized_match"] += 1

    metrics = _metrics(totals, family)
    for stratum in strata:
        for split in ("snp", "indel", "het", "hom_alt"):
            by_stratum.setdefault(f"{stratum.name}:{split}", _empty_counts())
    metrics["by_stratum"] = {
        name: _metrics(by_stratum[name], family) for name in sorted(by_stratum)
    }
    record = {
        "evaluation_schema_version": EVALUATION_SCHEMA_VERSION,
        "evaluation_id": evaluation_id,
        "reference": reference_description,
        "region_set": region.name,
        "region_definition_hash": region.definition_hash(reference_sha),
        "evidence_class": evidence_class,
        "query_dataset": query["dataset"],
        "comparator_dataset": comparator["dataset"],
        "comparison_unit": "normalized biallelic ALT allele (contig, pos, ref, alt) per sample",
        "eligible_denominator": eligible,
        "excluded_denominator": sum(exclusions.values()),
        "exclusion_reasons": dict(sorted(exclusions.items())),
        "not_evaluated_reason": None,
        "metrics": metrics,
        "strata": [
            {
                "name": s.name,
                "kind": s.kind,
                "group": s.group,
                "bases": s.region.bases(),
                "region_definition_hash": s.region.definition_hash(reference_sha),
                "spanning_boundary_units_not_counted": boundary_by_stratum[s.name],
            }
            for s in strata
        ]
        + [
            {"name": "snp", "kind": "partition", "group": "variant_type"},
            {"name": "indel", "kind": "partition", "group": "variant_type"},
            {"name": "het", "kind": "partition", "group": "alt_dosage"},
            {"name": "hom_alt", "kind": "partition", "group": "alt_dosage"},
        ],
        "strata_note": (
            "tag strata overlap and their counts are not additive; partition strata within one "
            "group are disjoint and add up to the total for units wholly inside a member; "
            "'<stratum>:snp'/'<stratum>:indel' split each stratum by variant type and '<stratum>:het'/'<stratum>:hom_alt' "
            "by ALT dosage; 'het' and 'hom_alt' "
            "split units by the comparator's ALT dosage (1 or 2), or the query's when the comparator is not positive"
        ),
        "definitions": {
            "true_negatives": "not counted: a VCF does not enumerate confident reference positions",
            "precision_or_rate_denominator_zero": "reported as null, never as 0 or 1",
            "excluded_units": (
                "one unit per (contig, pos, ref, alt) as written; a unit both callsets carry is "
                "counted once"
            ),
            "exclusion_precedence": list(EXCLUSION_PRECEDENCE),
            "excluded_units_by_side": {reason: by_side[reason] for reason in sorted(by_side)},
            "duplicate_normalized_keys": {
                "query": q.duplicate_keys,
                "comparator": c.duplicate_keys,
            },
            "input_records": {"query": q.records, "comparator": c.records},
        },
        "limitations": limitations,
    }
    validate_evaluation(record)
    return record


def not_evaluated(
    *,
    evaluation_id: str,
    evidence_class: str,
    reference_description: dict[str, object],
    region_set: str,
    reason: str,
    limitations: list[str],
) -> dict[str, object]:
    record = {
        "evaluation_schema_version": EVALUATION_SCHEMA_VERSION,
        "evaluation_id": evaluation_id,
        "reference": reference_description,
        "region_set": region_set,
        "region_definition_hash": None,
        "evidence_class": evidence_class,
        "query_dataset": None,
        "comparator_dataset": None,
        "comparison_unit": None,
        "eligible_denominator": None,
        "excluded_denominator": None,
        "exclusion_reasons": None,
        "not_evaluated_reason": reason,
        "metrics": {},
        "limitations": limitations,
    }
    validate_evaluation(record)
    return record
