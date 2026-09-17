#!/usr/bin/env python3
"""Descriptive stratification of a cohort callset by region (Issue #65).

This is the `descriptive_stratification` evidence class: counts and rates per
region stratum, with no quality claim. It reads the text that

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT]\\n' <vcf>

prints, so the VCF parsing is bcftools' own. With `--masked-calls` (the same
query over the Issue #64 genotype-masked VCF derived from the same callset) it
also reports, per stratum, how many genotype cells the mask removed; the two
streams must list the identical records in the identical order, or it fails.

That comparison is "site filter only" versus "site filter + genotype mask" on
the *same* records and samples. It does not mix in sample QC, and a difference
between strata describes where the mask acts, not which stratum is accurate.

FILTER is reported in three buckets, because a VCF distinguishes them: `PASS`
(filters applied and passed), a named code (failed), and `.` (no filter
applied, as in a raw single-sample callset).

Every stratum is also split by variant type, as `<stratum>:snp` and
`<stratum>:indel`, with the same base count: a callset's SNP density and its
indel density are different statements about one region, and a matrix row for
SNPs must not quote a figure that also counts indels. Records are not split by
genotype dosage, because a cohort record has one dosage per sample.

A record is placed in a stratum by its REF span. A record that straddles a
stratum's edge is counted for neither side of it and reported as a boundary
count for that stratum. A record outside the universe is excluded with a reason.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

from evidence_model import EVALUATION_SCHEMA_VERSION, validate_evaluation
from regions import MalformedRegionError, RegionSet, Stratum, read_fai, validate_partition


class MalformedCallsError(ValueError):
    """The bcftools query stream is not what this tool reads."""


def _open_text(path: str):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii")
    return open(path, encoding="ascii")


def _parse(line: str, number: int, samples: int | None):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        raise MalformedCallsError(f"line {number}: expected CHROM POS REF ALT FILTER and genotypes")
    contig, pos, ref, alt, flt = fields[:5]
    if not pos.isdigit() or not ref or not alt:
        raise MalformedCallsError(f"line {number}: malformed site fields")
    gts = fields[5:]
    if samples is not None and len(gts) != samples:
        raise MalformedCallsError(f"line {number}: {len(gts)} genotypes, expected {samples}")
    missing = het = hom_alt = 0
    for gt in gts:
        alleles = gt.replace("|", "/").split("/")
        if "." in alleles:
            missing += 1
        elif len(set(alleles)) > 1:
            het += 1
        elif alleles[0] != "0":
            hom_alt += 1
    return contig, int(pos), ref, alt, flt, len(gts), missing, (het, hom_alt)


def _filter_bucket(value: str) -> str:
    """VCF FILTER: `PASS` passed the filters and `.` means none were applied."""
    if value == "PASS":
        return "pass_records"
    if value in (".", ""):
        return "unfiltered_records"
    return "failed_filter_records"


def _empty() -> dict[str, int]:
    return dict.fromkeys(
        (
            "variant_records",
            "snp_records",
            "indel_records",
            "pass_records",
            "failed_filter_records",
            "unfiltered_records",
            "genotype_cells",
            "missing_genotype_cells",
            "non_reference_calls",
            "heterozygous_calls",
            "homozygous_alt_calls",
            "masked_genotype_cells",
        ),
        0,
    )


def _finish(
    counts: dict[str, int], bases: int, universe_bases: int, masked: bool
) -> dict[str, object]:
    out: dict[str, object] = {
        k: v for k, v in counts.items() if masked or k != "masked_genotype_cells"
    }
    out["bases"] = bases
    out["fraction_of_universe"] = round(bases / universe_bases, 6) if universe_bases else None
    cells = counts["genotype_cells"]
    out["missing_genotype_fraction"] = (
        round(counts["missing_genotype_cells"] / cells, 6) if cells else None
    )
    if masked:
        out["masked_genotype_fraction"] = (
            round(counts["masked_genotype_cells"] / cells, 6) if cells else None
        )
    non_ref = counts["non_reference_calls"]
    out["heterozygous_fraction_of_non_reference_calls"] = (
        round(counts["heterozygous_calls"] / non_ref, 6) if non_ref else None
    )
    out["variant_records_per_mb"] = (
        round(counts["variant_records"] / bases * 1e6, 3) if bases else None
    )
    return out


def stratify(
    *,
    calls,
    masked_calls,
    universe: RegionSet,
    strata: list[Stratum],
    reference_sha256: str,
    evaluation_id: str,
    query_dataset: dict[str, object],
    limitations: list[str],
    evidence_class: str = "descriptive_stratification",
) -> dict[str, object]:
    if evidence_class not in ("descriptive_stratification", "reference_sample_self_consistency"):
        raise ValueError(f"{evidence_class} does not carry descriptive per-stratum counts")
    groups: dict[str, list[Stratum]] = {}
    for stratum in strata:
        if stratum.kind == "partition":
            groups.setdefault(stratum.group, []).append(stratum)
    for members in groups.values():
        validate_partition([m.region for m in members], universe)

    totals = _empty()
    by_stratum = {s.name: _empty() for s in strata}
    by_stratum.update({"snp": _empty(), "indel": _empty()})
    by_stratum.update({f"{s.name}:{t}": _empty() for s in strata for t in ("snp", "indel")})
    boundary = {s.name: 0 for s in strata}
    exclusions: dict[str, int] = {}
    samples = None
    eligible = 0
    masked_iter = iter(masked_calls) if masked_calls is not None else None
    for number, line in enumerate(calls, start=1):
        contig, pos, ref, alt, flt, n, missing, non_ref = _parse(line, number, samples)
        samples = n
        masked_cells = 0
        if masked_iter is not None:
            other = next(masked_iter, None)
            if other is None:
                raise MalformedCallsError(f"masked stream ended before line {number}")
            m_contig, m_pos, m_ref, m_alt, m_flt, _n, m_missing, _nr = _parse(
                other, number, samples
            )
            if (m_contig, m_pos, m_ref, m_alt, m_flt) != (contig, pos, ref, alt, flt):
                raise MalformedCallsError(f"line {number}: masked stream lists a different record")
            if m_missing < missing:
                raise MalformedCallsError(f"line {number}: masked stream has fewer missing calls")
            masked_cells = m_missing - missing
        placement = universe.contains_span(contig, pos, len(ref))
        if placement != "inside":
            reason = "outside_universe" if placement == "outside" else "universe_boundary"
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        eligible += 1
        alts = alt.split(",")
        is_snp = len(ref) == 1 and all(len(a) == 1 for a in alts)
        variant_type = "snp" if is_snp else "indel"
        targets = [totals, by_stratum[variant_type]]
        for stratum in strata:
            where = stratum.region.contains_span(contig, pos, len(ref))
            if where == "inside":
                targets.append(by_stratum[stratum.name])
                targets.append(by_stratum[f"{stratum.name}:{variant_type}"])
            elif where == "boundary":
                boundary[stratum.name] += 1
        for counts in targets:
            counts["variant_records"] += 1
            counts["snp_records" if is_snp else "indel_records"] += 1
            counts[_filter_bucket(flt)] += 1
            counts["genotype_cells"] += n
            counts["missing_genotype_cells"] += missing
            counts["non_reference_calls"] += non_ref[0] + non_ref[1]
            counts["heterozygous_calls"] += non_ref[0]
            counts["homozygous_alt_calls"] += non_ref[1]
            counts["masked_genotype_cells"] += masked_cells
    if masked_iter is not None and next(masked_iter, None) is not None:
        raise MalformedCallsError("masked stream has more records than the callset")

    masked = masked_calls is not None
    metrics = _finish(totals, universe.bases(), universe.bases(), masked)
    metrics["samples"] = samples or 0
    bases_of = {s.name: s.region.bases() for s in strata}
    bases_of.update({"snp": universe.bases(), "indel": universe.bases()})
    bases_of.update({f"{s.name}:{t}": bases_of[s.name] for s in strata for t in ("snp", "indel")})
    metrics["by_stratum"] = {
        name: _finish(by_stratum[name], bases_of[name], universe.bases(), masked)
        for name in sorted(by_stratum)
    }
    record = {
        "evaluation_schema_version": EVALUATION_SCHEMA_VERSION,
        "evaluation_id": evaluation_id,
        "reference": {"sha256": reference_sha256},
        "region_set": universe.name,
        "region_definition_hash": universe.definition_hash(reference_sha256),
        "evidence_class": evidence_class,
        "query_dataset": query_dataset,
        "comparator_dataset": None,
        "comparison_unit": "VCF record placed by its REF span; genotype cell (record x sample)",
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
                "region_definition_hash": s.region.definition_hash(reference_sha256),
                "records_straddling_edge_not_counted": boundary[s.name],
            }
            for s in strata
        ]
        + [
            {"name": "snp", "kind": "partition", "group": "variant_type"},
            {"name": "indel", "kind": "partition", "group": "variant_type"},
        ],
        "strata_note": (
            "tag strata overlap and are not additive; members of one partition group are disjoint, "
            "cover the universe, and their counts add up to the total less the edge-straddling records"
        ),
        "definitions": {
            "masked_genotype_cells": "missing in the masked stream but not in the unmasked stream, same record and sample",
            "missing_genotype_cells": "a GT with any missing allele in the unmasked stream",
            "filter_buckets": (
                "pass_records: FILTER=PASS; failed_filter_records: a named filter code; "
                "unfiltered_records: FILTER='.', meaning no filter was applied"
            ),
            "variant_type_split": (
                "'<stratum>:snp' and '<stratum>:indel' count the same region's records by type; "
                "their bases are the stratum's bases, so each density is per type"
            ),
            "non_reference_calls": "a fully called GT with at least one non-reference allele",
            "heterozygous_calls": "a fully called GT with two different alleles",
            "homozygous_alt_calls": "a fully called GT whose alleles are one and the same non-reference allele",
        },
        "limitations": limitations,
    }
    validate_evaluation(record)
    return record


def load_strata(
    config: dict[str, object], base: Path, contigs, universe: RegionSet
) -> list[Stratum]:
    strata = []
    for entry in config["strata"]:
        region = RegionSet.from_bed(base / entry["bed"], entry["name"], contigs)
        if entry.get("restrict_to_universe", True):
            region = region.intersect(universe, entry["name"])
            region.provenance = {"bed": entry["bed"], "intersected_with": universe.name}
        strata.append(Stratum(entry["name"], entry["kind"], entry["group"], region))
    return strata


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="JSON: evaluation_id, universe (bed or 'whole_reference'), strata, limitations, query_dataset",
    )
    parser.add_argument("--reference-fai", required=True, type=Path)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--masked-calls")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    config = json.loads(args.config.read_text())
    base = args.config.parent
    contigs = read_fai(args.reference_fai)
    if config["universe"] == "whole_reference":
        universe = RegionSet.whole_reference("whole_reference", contigs)
    else:
        universe = RegionSet.from_bed(
            base / config["universe"]["bed"], config["universe"]["name"], contigs
        )
    try:
        strata = load_strata(config, base, contigs, universe)
    except MalformedRegionError as error:
        raise SystemExit(f"region error: {error}") from error
    calls = _open_text(args.calls)
    masked = _open_text(args.masked_calls) if args.masked_calls else None
    record = stratify(
        calls=calls,
        masked_calls=masked,
        universe=universe,
        strata=strata,
        reference_sha256=args.reference_sha256,
        evaluation_id=config["evaluation_id"],
        query_dataset=config["query_dataset"],
        limitations=config["limitations"],
        evidence_class=config.get("evidence_class", "descriptive_stratification"),
    )
    record["config_sha256"] = hashlib.sha256(args.config.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(record, indent=2, sort_keys=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
