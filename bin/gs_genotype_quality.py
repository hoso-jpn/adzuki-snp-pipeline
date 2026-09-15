"""Genotype-level DP/GQ quality policy for the GS panel (Issue #64).

This is the single source of truth for *which* genotype calls a quality
policy masks and *why*. `build_gs_panel.py` uses it to write the matrix, the
accounting and the quality-masked VCF in one pass, and
`verify_gs_genotype_quality_mask.py` uses the same evaluator to check that
every published artifact agrees with it. Neither script carries its own
masking rule, so the matrix and the VCF cannot drift apart by construction.

The policy is optional and off by default. It is a separate stage from the
site-level hard filter: a PASS site can still carry individual calls with
little read support, and this stage lets a caller turn such calls into
missing genotypes without touching the site, the sample list, or the site
FILTER. No threshold is recommended here; `min_dp` and `min_gq` are only
ever what the caller explicitly asked for.

## What is evaluated

Only calls the GS encoding would otherwise turn into a dosage -- diploid,
biallelic-index, no missing allele (`0/0`, `0/1`, `1/0`, `1/1`, either
phasing) -- are evaluated. A call that is already missing or non-standard
is `nan` in the matrix regardless, and is counted under its existing reason
rather than being evaluated a second time. Phasing never changes how a call
is evaluated: `0|1` and `0/1` meet the same thresholds.

## Per-field status

For each configured field (`DP` when `min_dp` is set, `GQ` when `min_gq` is
set) a call has exactly one status:

| status | meaning |
| --- | --- |
| `pass` | an integer value `>= threshold` (the threshold itself passes) |
| `below_threshold` | an integer value `< threshold` |
| `format_field_absent` | the row's FORMAT has no such key |
| `value_missing` | the key is present and the sample's value is `.` |
| `value_truncated` | the key is present but the sample carries too few subfields |
| `value_malformed` | anything that is not a non-negative decimal integer |

A missing or absent value is never read as `0` and never read as a pass.
What happens to it is the caller's explicit choice:

* `below_threshold` always masks the call;
* `format_field_absent` follows `missing_format_field`;
* `value_missing` and `value_truncated` follow `missing_value`;
* `value_malformed` follows `malformed_value`;

where each choice is `reject` (stop the build: the input cannot be judged
under this policy), `unevaluated` (keep the call and count that the field
could not be evaluated), or `mask` (treat the call as failing the policy).

## Why a masked call has a reason pair, not a single reason

A call can fail on DP and GQ at once, so the masked calls are partitioned
by the pair of per-field reasons, each one of `none`, `low`, `unavailable`
(absent, missing or truncated, masked by policy) or `malformed`. Every
masked call falls in exactly one pair, so the pair counts sum to the masked
total with no double counting. The per-field status counts are a second,
independent partition of the evaluated calls, one per configured field.
"""

from __future__ import annotations

from dataclasses import dataclass

# The same canonical hash the GS and run manifests use, so a policy hash can
# be compared across all three documents without a second convention.
from manifest_utils import canonical_json_hash

POLICY_SCHEMA = "gs_genotype_quality_policy_v1"

REJECT = "reject"
UNEVALUATED = "unevaluated"
MASK = "mask"
POLICY_ACTIONS: tuple[str, ...] = (REJECT, UNEVALUATED, MASK)

PASS = "pass"
BELOW_THRESHOLD = "below_threshold"
FORMAT_FIELD_ABSENT = "format_field_absent"
VALUE_MISSING = "value_missing"
VALUE_TRUNCATED = "value_truncated"
VALUE_MALFORMED = "value_malformed"
FIELD_STATUSES: tuple[str, ...] = (
    PASS,
    BELOW_THRESHOLD,
    FORMAT_FIELD_ABSENT,
    VALUE_MISSING,
    VALUE_TRUNCATED,
    VALUE_MALFORMED,
)

REASON_NONE = "none"
REASON_LOW = "low"
REASON_UNAVAILABLE = "unavailable"
REASON_MALFORMED = "malformed"
MASK_REASONS: tuple[str, ...] = (REASON_NONE, REASON_LOW, REASON_UNAVAILABLE, REASON_MALFORMED)

#: FORMAT keys a policy can evaluate, in the order they appear everywhere.
QUALITY_FIELDS: tuple[str, ...] = ("DP", "GQ")

#: The INFO keys recomputed from the masked genotypes on rows that changed.
RECOMPUTED_INFO_KEYS: tuple[str, ...] = ("AC", "AN", "AF")

DISABLED_POLICY_DOCUMENT: dict[str, object] = {
    "schema": POLICY_SCHEMA,
    "enabled": False,
}


class GenotypeQualityPolicyError(Exception):
    """A call cannot be judged under the chosen policy, and the policy says reject."""


class InvalidPolicyError(ValueError):
    """The policy itself is not a valid, explicit policy."""


def _validate_threshold(name: str, value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise InvalidPolicyError(f"{name} must be a non-negative integer or null, got {value!r}")
    return value


def _validate_action(name: str, value: object) -> str:
    if value not in POLICY_ACTIONS:
        raise InvalidPolicyError(f"{name} must be one of {list(POLICY_ACTIONS)}, got {value!r}")
    return value  # type: ignore[return-value]


@dataclass(frozen=True)
class GenotypeQualityPolicy:
    """An enabled, explicit genotype quality policy."""

    min_dp: int | None
    min_gq: int | None
    missing_format_field: str = REJECT
    missing_value: str = REJECT
    malformed_value: str = REJECT

    def __post_init__(self) -> None:
        _validate_threshold("min_dp", self.min_dp)
        _validate_threshold("min_gq", self.min_gq)
        if self.min_dp is None and self.min_gq is None:
            raise InvalidPolicyError(
                "an enabled genotype quality policy needs min_dp, min_gq, or both; "
                "with neither there is nothing to evaluate"
            )
        _validate_action("missing_format_field", self.missing_format_field)
        _validate_action("missing_value", self.missing_value)
        _validate_action("malformed_value", self.malformed_value)

    @property
    def thresholds(self) -> tuple[tuple[str, int], ...]:
        """The configured (FORMAT key, threshold) pairs, DP before GQ."""
        pairs = (("DP", self.min_dp), ("GQ", self.min_gq))
        return tuple((key, value) for key, value in pairs if value is not None)

    def document(self) -> dict[str, object]:
        """The machine-readable policy, as published and hashed."""
        return {
            "schema": POLICY_SCHEMA,
            "enabled": True,
            "min_dp": self.min_dp,
            "min_gq": self.min_gq,
            "comparison": "a value >= its threshold passes; a value < its threshold masks the call",
            "evaluated_calls": (
                "diploid biallelic-index calls with no missing allele, phased or unphased; "
                "calls already missing or non-standard are not evaluated"
            ),
            "missing_format_field": self.missing_format_field,
            "missing_value": self.missing_value,
            "missing_value_covers": [VALUE_MISSING, VALUE_TRUNCATED],
            "malformed_value": self.malformed_value,
            "valid_value": "a non-negative decimal integer",
            "masked_genotype": (
                "GT alleles replaced by '.', keeping the call's own separator; "
                "every other FORMAT subfield is retained unchanged"
            ),
            "recomputed_info": list(RECOMPUTED_INFO_KEYS),
            "recomputed_info_scope": (
                "only on rows where at least one call was masked; every other INFO "
                "annotation describes the original call set and is not recomputed"
            ),
        }

    def policy_hash(self) -> str:
        return canonical_json_hash(self.document())


def disabled_policy_hash() -> str:
    return canonical_json_hash(DISABLED_POLICY_DOCUMENT)


def policy_from_document(document: dict[str, object]) -> GenotypeQualityPolicy | None:
    """Rebuild a policy from its published document, refusing anything else.

    Returns None for the disabled document. The document must be exactly
    what `GenotypeQualityPolicy.document()` produces for the same values, so
    a hand-edited or partially written policy file cannot be recorded as if
    it were the policy that ran.
    """
    if not isinstance(document, dict) or document.get("schema") != POLICY_SCHEMA:
        raise InvalidPolicyError(f"not a {POLICY_SCHEMA} document")
    if document == DISABLED_POLICY_DOCUMENT:
        return None
    if document.get("enabled") is not True:
        raise InvalidPolicyError("a policy document must be the disabled document or enabled=true")
    try:
        policy = GenotypeQualityPolicy(
            min_dp=document.get("min_dp"),  # type: ignore[arg-type]
            min_gq=document.get("min_gq"),  # type: ignore[arg-type]
            missing_format_field=document.get("missing_format_field"),  # type: ignore[arg-type]
            missing_value=document.get("missing_value"),  # type: ignore[arg-type]
            malformed_value=document.get("malformed_value"),  # type: ignore[arg-type]
        )
    except TypeError as error:
        raise InvalidPolicyError(str(error)) from error
    if policy.document() != document:
        raise InvalidPolicyError("policy document does not match its own canonical form")
    return policy


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class CallQuality:
    """The policy's verdict on one evaluated call."""

    masked: bool
    unevaluated: bool
    #: One status per configured field, in `policy.thresholds` order.
    statuses: tuple[str, ...]
    #: One mask reason per configured field, in the same order.
    reasons: tuple[str, ...]


def _field_status(value: str | None, key_present: bool, threshold: int) -> str:
    if not key_present:
        return FORMAT_FIELD_ABSENT
    if value is None:
        return VALUE_TRUNCATED
    if value == ".":
        return VALUE_MISSING
    if not value.isascii() or not value.isdigit():
        return VALUE_MALFORMED
    return PASS if int(value) >= threshold else BELOW_THRESHOLD


class RowQualityEvaluator:
    """Evaluate every call of one VCF row against a policy.

    Built once per row, because FORMAT -- and therefore where DP and GQ sit
    -- is a per-row property of VCF.
    """

    def __init__(self, policy: GenotypeQualityPolicy, format_field: str) -> None:
        keys = format_field.split(":")
        self._policy = policy
        self._fields = tuple(
            (key, threshold, keys.index(key) if key in keys else None)
            for key, threshold in policy.thresholds
        )

    def evaluate(self, subfields: list[str]) -> CallQuality:
        """Judge one call from its already-split sample subfields.

        A rejection names only the field and why; the caller adds the line
        and sample when it catches it, so no per-call context string is
        built on the path that almost never raises.
        """
        policy = self._policy
        statuses: list[str] = []
        reasons: list[str] = []
        masked = False
        unevaluated = False
        for key, threshold, index in self._fields:
            if index is None:
                status = _field_status(None, False, threshold)
            else:
                value = subfields[index] if index < len(subfields) else None
                status = _field_status(value, True, threshold)
            statuses.append(status)

            if status == PASS:
                reasons.append(REASON_NONE)
                continue
            if status == BELOW_THRESHOLD:
                reasons.append(REASON_LOW)
                masked = True
                continue

            if status == FORMAT_FIELD_ABSENT:
                action, reason = policy.missing_format_field, REASON_UNAVAILABLE
            elif status == VALUE_MALFORMED:
                action, reason = policy.malformed_value, REASON_MALFORMED
            else:
                action, reason = policy.missing_value, REASON_UNAVAILABLE

            if action == REJECT:
                raise GenotypeQualityPolicyError(
                    f"{key} is {status} and the policy rejects it "
                    f"({_policy_option_for(status)}=reject)"
                )
            if action == MASK:
                reasons.append(reason)
                masked = True
            else:
                reasons.append(REASON_NONE)
                unevaluated = True

        return CallQuality(
            masked=masked,
            unevaluated=unevaluated and not masked,
            statuses=tuple(statuses),
            reasons=tuple(reasons),
        )


def _policy_option_for(status: str) -> str:
    if status == FORMAT_FIELD_ABSENT:
        return "missing_format_field"
    if status == VALUE_MALFORMED:
        return "malformed_value"
    return "missing_value"


# --------------------------------------------------------------------------
# Accounting names
# --------------------------------------------------------------------------


def reason_metric(policy: GenotypeQualityPolicy, reasons: tuple[str, ...]) -> str:
    """`quality_masked_reason.dp_low.gq_none` for a masked call's reason pair."""
    parts = [f"{key.lower()}_{reason}" for (key, _), reason in zip(policy.thresholds, reasons)]
    return "quality_masked_reason." + ".".join(parts)


def reason_metrics(policy: GenotypeQualityPolicy) -> tuple[str, ...]:
    """Every reason-pair metric a masked call can fall in, in a fixed order."""
    combos: list[tuple[str, ...]] = [()]
    for _ in policy.thresholds:
        combos = [(*combo, reason) for combo in combos for reason in MASK_REASONS]
    return tuple(
        reason_metric(policy, combo) for combo in combos if any(r != REASON_NONE for r in combo)
    )


def status_metric(key: str, status: str) -> str:
    return f"{key.lower()}_status.{status}"


def status_metrics(policy: GenotypeQualityPolicy) -> tuple[str, ...]:
    return tuple(
        status_metric(key, status) for key, _ in policy.thresholds for status in FIELD_STATUSES
    )


# --------------------------------------------------------------------------
# VCF rewriting
# --------------------------------------------------------------------------


def masked_genotype(gt: str) -> str:
    """`./.` for an unphased call, `.|.` for a phased one."""
    return ".|." if "|" in gt else "./."


def allele_counts(genotypes: list[str]) -> tuple[int, int]:
    """(AC, AN) over every called allele of every genotype, biallelic ALT only."""
    ac = an = 0
    for gt in genotypes:
        for allele in gt.replace("|", "/").split("/"):
            if allele in (".", ""):
                continue
            an += 1
            if allele == "1":
                ac += 1
    return ac, an


def format_allele_frequency(ac: int, an: int) -> str:
    return "." if an == 0 else f"{ac / an:.6f}"


def rewrite_allele_info(info: str, ac: int, an: int) -> str:
    """Replace AC, AN and AF in an INFO column, keeping every other entry in place."""
    replacements = {"AC": str(ac), "AN": str(an), "AF": format_allele_frequency(ac, an)}
    seen: set[str] = set()
    entries = []
    for entry in info.split(";"):
        key = entry.split("=", 1)[0]
        if key in replacements:
            if key in seen:
                raise ValueError(f"INFO carries {key} more than once")
            seen.add(key)
            entries.append(f"{key}={replacements[key]}")
        else:
            entries.append(entry)
    missing = [key for key in RECOMPUTED_INFO_KEYS if key not in seen]
    if missing:
        raise ValueError(f"INFO lacks {', '.join(missing)}, which masking must recompute")
    return ";".join(entries)


def mask_header_line(policy: GenotypeQualityPolicy) -> str:
    """The provenance line the quality-masked VCF carries before #CHROM."""

    def show(value: int | None) -> str:
        return "." if value is None else str(value)

    return (
        f"##gs_genotype_quality_mask=<Schema={POLICY_SCHEMA},PolicyHash={policy.policy_hash()},"
        f"MinDP={show(policy.min_dp)},MinGQ={show(policy.min_gq)},"
        f"MissingFormatField={policy.missing_format_field},MissingValue={policy.missing_value},"
        f"MalformedValue={policy.malformed_value},"
        'Description="Calls failing this policy have GT set to missing; FORMAT is retained; '
        'AC, AN and AF are recomputed on rows with masked calls">'
    )
