#!/usr/bin/env python3
"""Assemble the Issue #65 evaluation records and the delivery-support matrix.

    assemble_evidence.py derive-strata --workdir W --reference-fasta ref.fna
    assemble_evidence.py assemble      --workdir W --reference-fasta ref.fna --git-sha SHA

`derive-strata` writes two composite regions over the evaluation window and
checks that they partition it:

* `core`      = cohort_callable AND high_mappability, minus WindowMasker repeats
                and homopolymers >= 10;
* `difficult` = the rest of the window.

`assemble` reads what the run scripts left under the work directory (see
docs/issue65_callable_region_quality_benchmark.md for the layout), computes the
comparison-shaped evaluations with `compare.evaluate`, loads the descriptive
ones written by `stratify_callset.py`, adds the not-evaluated records, and
derives the matrix. Every record passes `evidence_model.validate_evaluation`
and every matrix cell passes `evidence_model.delivery_row`, so a status
stronger than its evidence cannot be written. Output contains no timestamps:
running it twice on the same inputs gives byte-identical files.

Matrix rule (no numeric threshold is introduced):

* a cell for `independent_truth` or `technical_replicate_concordance` is
  `not_evaluated` -- the audit found no usable asset for either;
* a comparison or descriptive cell is `supported_with_caveat` when its
  evaluation has at least one unit in the scope, and `not_evaluated` otherwise;
* a scope's overall status is `supported_with_caveat` when any cell is, and
  `not_evaluated` otherwise. No scope can be `supported` without an independent
  truth cell.

No benchmark stratification parameter becomes a delivery gate here. The
callable rule of `callable_regions.py` was chosen to stratify this benchmark and
was never calibrated against an outcome, so a scope that fails it is *flagged* --
the row carries `fails_benchmark_callable_rule` and a caveat naming the rule --
and not marked `unsupported`, which would adopt that rule as a threshold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import compare
import evidence_model as model
from regions import RegionSet, Stratum, read_fai, validate_partition
from variants import Reference

WINDOW = ("NC_068975.1", 0, 20_000_000)
WINDOW_NAME = "NC_068975.1:1-20000000"
DOWNSAMPLING_SAMPLE = "SRR29908806"
DOWNSAMPLES = {
    "f050_s65": {"fraction": 0.5, "seed": 65},
    "f050_s66": {"fraction": 0.5, "seed": 66},
    "f025_s65": {"fraction": 0.25, "seed": 65},
    "f025_s66": {"fraction": 0.25, "seed": 66},
}
GC = ["gc_00_25", "gc_25_30", "gc_30_35", "gc_35_40", "gc_40_45", "gc_45_100", "gc_undefined"]
MAPPABILITY = ["high_mappability", "low_mappability", "mappability_undefined"]
DEPTH = [
    "median_depth_00_03",
    "median_depth_03_06",
    "median_depth_06_10",
    "median_depth_10_15",
    "median_depth_15_25",
    "median_depth_25_plus",
]
CONTAINERS = {
    "bcftools": "quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3",
    "gatk": "broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb",
    "samtools": "quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561",
    "bwa_mem2_samtools": "community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113@sha256:5ebd1290d9680195817ce75915b79ae2e608834c017824b7e2bc7b141509b242",
    "fastp": "quay.io/biocontainers/fastp:1.3.6--h43da1c4_0@sha256:cbbe2402b6b6704df470d7d77dcb498eefd5bcd01f4c38be0ec69899e79ac134",
    "python": "python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _window(contigs) -> RegionSet:
    return RegionSet(WINDOW_NAME, contigs, {WINDOW[0]: [(WINDOW[1], WINDOW[2])]})


def _load(workdir: Path, relative: str, name: str, contigs, window: RegionSet) -> RegionSet:
    region = RegionSet.from_bed(workdir / relative, name, contigs).intersect(window, name)
    region.provenance = {"bed": relative, "intersected_with": WINDOW_NAME}
    return region


def window_strata(workdir: Path, contigs) -> tuple[RegionSet, list[Stratum]]:
    window = _window(contigs)
    strata: list[Stratum] = []
    for name in ("cohort_callable", "cohort_non_callable"):
        strata.append(
            Stratum(
                name,
                "partition",
                "cohort_callable",
                _load(workdir, f"callable/{name}.bed", name, contigs, window),
            )
        )
    for name in DEPTH:
        strata.append(
            Stratum(
                name,
                "partition",
                "cohort_median_depth",
                _load(workdir, f"callable/{name}.bed", name, contigs, window),
            )
        )
    for name in MAPPABILITY:
        strata.append(
            Stratum(
                name,
                "partition",
                "mappability",
                _load(workdir, f"mappability/{name}.bed", name, contigs, window),
            )
        )
    for name in GC:
        strata.append(
            Stratum(
                name,
                "partition",
                "gc_1kb",
                _load(workdir, f"regions/{name}.bed", name, contigs, window),
            )
        )
    for name, group in (("repeat_windowmasker", "repeat"), ("homopolymer_ge10", "low_complexity")):
        strata.append(
            Stratum(
                name, "tag", group, _load(workdir, f"regions/{name}.bed", name, contigs, window)
            )
        )
    for name in ("core", "difficult"):
        strata.append(
            Stratum(
                name,
                "partition",
                "core_vs_difficult",
                _load(workdir, f"derived/{name}.bed", name, contigs, window),
            )
        )
    groups: dict[str, list[RegionSet]] = {}
    for stratum in strata:
        if stratum.kind == "partition":
            groups.setdefault(stratum.group, []).append(stratum.region)
    for members in groups.values():
        validate_partition(members, window)
    return window, strata


def derive_strata(workdir: Path, reference_fasta: Path) -> dict[str, object]:
    contigs = read_fai(Path(str(reference_fasta) + ".fai"))
    window = _window(contigs)
    callable_ = _load(workdir, "callable/cohort_callable.bed", "cohort_callable", contigs, window)
    high = _load(workdir, "mappability/high_mappability.bed", "high_mappability", contigs, window)
    repeat = _load(
        workdir, "regions/repeat_windowmasker.bed", "repeat_windowmasker", contigs, window
    )
    homopolymer = _load(
        workdir, "regions/homopolymer_ge10.bed", "homopolymer_ge10", contigs, window
    )
    core = callable_.intersect(high, "core").subtract(repeat, "core").subtract(homopolymer, "core")
    difficult = window.subtract(core, "difficult")
    validate_partition([core, difficult], window)
    out = workdir / "derived"
    out.mkdir(exist_ok=False)
    (out / "core.bed").write_text(core.to_bed())
    (out / "difficult.bed").write_text(difficult.to_bed())
    summary = {
        "window": WINDOW_NAME,
        "definition": {
            "core": "cohort_callable AND high_mappability AND NOT repeat_windowmasker AND NOT homopolymer_ge10",
            "difficult": "window AND NOT core",
        },
        "bases": {"core": core.bases(), "difficult": difficult.bases(), "window": window.bases()},
        "sha256": {
            "core.bed": sha256(out / "core.bed"),
            "difficult.bed": sha256(out / "difficult.bed"),
        },
    }
    (out / "derived_strata.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def _sha_table(*paths: Path) -> dict[str, str]:
    table = {}
    for path in paths:
        if path.exists():
            for line in path.read_text().splitlines():
                digest, name = line.split(maxsplit=1)
                table[Path(name).name] = digest
    return table


def _ref_description(reference_fasta: Path, reference_sha: str) -> dict[str, object]:
    return {
        "accession": "GCF_016808095.1",
        "assembly": "ASM1680809v1",
        "fasta": Path(reference_fasta).name,
        "sha256": reference_sha,
    }


def comparison_evaluations(
    workdir: Path, reference_fasta: Path, reference_sha: str
) -> list[dict[str, object]]:
    contigs = read_fai(Path(str(reference_fasta) + ".fai"))
    window, strata = window_strata(workdir, contigs)
    d = workdir / "downsampling"
    shas = _sha_table(d / "outputs.sha256", d / "inputs.sha256")
    ref_desc = _ref_description(reference_fasta, reference_sha)

    def dataset(vcf: str, **extra) -> dict[str, object]:
        return {
            "path": f"downsampling/{vcf}",
            "sha256": shas.get(vcf),
            "sample": DOWNSAMPLING_SAMPLE,
            **extra,
        }

    def reads(id_: str) -> dict[str, object]:
        return {
            "reads_in_window": int((d / f"{id_}.reads.txt").read_text().split()[0]),
            "mean_depth_mapq20_bq10": float((d / f"{id_}.mean_depth.txt").read_text().split()[0]),
        }

    common_limits = [
        f"one sample ({DOWNSAMPLING_SAMPLE}, the most deeply sequenced of the 51) on one 20 Mb window; other samples and the rest of the genome are not represented",
        "single-sample GenotypeGVCFs without the cohort hard filter: these are raw calls, not the delivered GS panel",
        "no independent truth: nothing here measures accuracy",
    ]
    full = dataset(
        "full.vcf.gz",
        description="HaplotypeCaller (production arguments, --intervals window) on the full-depth markdup BAM",
        **reads("full"),
    )
    evaluations = []
    reference = Reference(reference_fasta)
    try:
        for id_, spec in DOWNSAMPLES.items():
            evaluations.append(
                compare.evaluate(
                    evaluation_id=f"downsampling.{DOWNSAMPLING_SAMPLE}.{id_}",
                    evidence_class="downsampling_stability",
                    reference=reference,
                    reference_description=ref_desc,
                    region=window,
                    strata=strata,
                    query=dataset_with_path(
                        workdir,
                        dataset(
                            f"{id_}.vcf.gz",
                            subsample=spec,
                            method="samtools view --subsample (read-name hash, mates kept together)",
                            **reads(id_),
                        ),
                    ),
                    comparator=dataset_with_path(workdir, full),
                    limitations=common_limits
                    + [
                        "the full-depth calls are the comparator, not truth: retention measures agreement with them",
                        "the two seeds at one fraction draw overlapping read sets; their difference shows sampling variation, not independent replication",
                    ],
                )
            )
        evaluations.append(
            compare.evaluate(
                evaluation_id=f"caller.{DOWNSAMPLING_SAMPLE}.bcftools_vs_haplotypecaller",
                evidence_class="caller_concordance",
                reference=reference,
                reference_description=ref_desc,
                region=window,
                strata=strata,
                query=dataset_with_path(
                    workdir,
                    dataset(
                        "full.bcftools.vcf.gz",
                        description="bcftools mpileup -q 20 -Q 10 | bcftools call -m -v (bcftools 1.24), same BAM and window",
                    ),
                ),
                comparator=dataset_with_path(workdir, full),
                limitations=common_limits
                + [
                    "both callers read the same alignments, so shared mapping errors agree with each other; agreement is not accuracy",
                    "bcftools call output is unfiltered and capped at 250 reads per position (mpileup default)",
                ],
            )
        )
    finally:
        reference.close()
    return evaluations


def dataset_with_path(workdir: Path, dataset: dict[str, object]) -> dict[str, object]:
    return {
        "path": workdir / str(dataset["path"]),
        "sample": dataset["sample"],
        "dataset": {k: v for k, v in dataset.items() if k != "sample"},
    }


def not_evaluated_records(reference_desc: dict[str, object]) -> list[dict[str, object]]:
    rows = [
        (
            "independent_truth.all",
            "independent_truth",
            "all",
            "no Vigna angularis small-variant truth set was found (ENA analyses for taxon 3914 are all REFERENCE_ALIGNMENT; EVA lists no Vigna angularis assembly); the reference BioSample's PacBio reads are not independent of the reference",
            ["see public_asset_inventory.json: vigna_angularis_truth_set, SRR11787766"],
        ),
        (
            "technical_replicate.all",
            "technical_replicate_concordance",
            "all",
            "no identity-confirmed technical replicate: every PRJNA1138464 BioSample has one run; SAMN03488483 aggregates a resequencing series and its same-library run pairs are mate-pair libraries or 4x-vs-20x 500 bp runs",
            ["see public_asset_inventory.json: PRJNA1138464, SAMN03488483"],
        ),
        (
            "cross_method.cohort_rad_vs_wgs",
            "cross_platform_concordance",
            "all",
            "PRJNA1138464 has no BioSample with both RAD-Seq and WGS; samples were not matched by name",
            ["see public_asset_inventory.json: PRJNA1138464"],
        ),
        (
            "long_read_calls.SAMN14776547",
            "cross_platform_concordance",
            "all",
            "no long-read caller is pinned by the pipeline, and the reference BioSample's long reads are what the assembly was built from",
            ["see public_asset_inventory.json: SRR11787766"],
        ),
        (
            "descriptive.gs_panel.indel",
            "descriptive_stratification",
            "all",
            "the GS panel is SNP-only by design, so there is no delivered indel callset to stratify",
            [
                "indel strata are evaluated only for single-sample raw calls (downsampling and caller evaluations)"
            ],
        ),
        (
            "callable.outside_window",
            "descriptive_stratification",
            "genome outside NC_068975.1:1-20000000",
            "callable and depth assets were built only for the Issue #45 scaling window from the existing 51 BAMs; extending them genome-wide was not run",
            [
                "whole-genome sequence-derived assets (N, repeat, homopolymer, GC, mappability) do exist"
            ],
        ),
        (
            "downsampling.other_samples_and_genome",
            "downsampling_stability",
            "genome outside the window, samples other than SRR29908806",
            "downsampling was run on one existing BAM and one window, not on further samples or the whole genome",
            ["one-sample, one-window result only"],
        ),
        (
            "repeat_family.all",
            "descriptive_stratification",
            "repeat families",
            "no RepeatMasker family annotation is published for GCF_016808095.1; WindowMasker lower-case marks repeats without family",
            ["see public_asset_inventory.json: refseq_repeatmasker_output"],
        ),
    ]
    return [
        compare.not_evaluated(
            evaluation_id=eid,
            evidence_class=cls,
            reference_description=reference_desc,
            region_set=region,
            reason=reason,
            limitations=limits,
        )
        for eid, cls, region, reason, limits in rows
    ]


# (scope id, compare stratum or None for the whole window, descriptive stratum); each scope is
# split into snp, indel (variant type) and het, hom_alt (ALT dosage of the comparator call)
SCOPES = [
    ("window", None, None),
    ("core", "core", "core"),
    ("difficult", "difficult", "difficult"),
    ("cohort_callable", "cohort_callable", "cohort_callable"),
    ("cohort_non_callable", "cohort_non_callable", "cohort_non_callable"),
    ("high_mappability", "high_mappability", "high_mappability"),
    ("low_mappability", "low_mappability", "low_mappability"),
    ("repeat_windowmasker", "repeat_windowmasker", "repeat_windowmasker"),
    ("homopolymer_ge10", "homopolymer_ge10", "homopolymer_ge10"),
    *[(name, name, name) for name in DEPTH],
    *[(name, name, name) for name in GC if name != "gc_undefined"],
]
# Flagged on the row, never turned into a status: these scopes fail the
# benchmark's own callable rule (see callable_regions.py), which stratifies this
# benchmark and was never calibrated as a delivery gate.
CALLABLE_RULE_CAVEAT = {
    "cohort_non_callable": (
        "fails the benchmark callable rule (depth 5 to 2.5x the sample's median, in at least 80% "
        "of the 51 samples); that rule stratifies this benchmark and is not calibrated as a "
        "delivery threshold"
    )
}


def _by(
    record: dict[str, object], stratum: str | None, variant_type: str
) -> dict[str, object] | None:
    by = record["metrics"]["by_stratum"]
    key = variant_type if stratum is None else f"{stratum}:{variant_type}"
    return by.get(key)


def _units(metrics: dict[str, object] | None, family: str) -> int:
    if metrics is None:
        return 0
    names = compare.METRIC_NAMES[family]
    return sum(int(metrics[names[k]]) for k in names)


def matrix(evaluations: list[dict[str, object]]) -> dict[str, object]:
    by_id = {e["evaluation_id"]: e for e in evaluations}
    downs = [
        e
        for e in evaluations
        if e["evidence_class"] == "downsampling_stability" and e["not_evaluated_reason"] is None
    ]
    caller = [
        e
        for e in evaluations
        if e["evidence_class"] == "caller_concordance" and e["not_evaluated_reason"] is None
    ]
    cross = [
        e
        for e in evaluations
        if e["evidence_class"] == "reference_sample_self_consistency"
        and e["not_evaluated_reason"] is None
    ]
    descriptive = by_id.get("gs_snp_pass_51.NC_068975.1_1-20000000.descriptive")
    rows, cells = [], []
    for scope, stratum, desc_stratum in SCOPES:
        for variant_type in ("snp", "indel", "het", "hom_alt"):
            scope_id = f"{scope}:{variant_type}"
            row_cells = []
            row_cells.append(
                model.delivery_row(
                    scope=scope_id,
                    evidence_class="independent_truth",
                    status="not_evaluated",
                    evidence_refs=[],
                    rationale="independent_truth.all",
                )
            )
            row_cells.append(
                model.delivery_row(
                    scope=scope_id,
                    evidence_class="technical_replicate_concordance",
                    status="not_evaluated",
                    evidence_refs=[],
                    rationale="technical_replicate.all",
                )
            )
            key_metrics: dict[str, object] = {}
            for cls, family, records in (
                ("downsampling_stability", "stability", downs),
                ("caller_concordance", "concordance", caller),
            ):
                present = [(e, _by(e, stratum, variant_type)) for e in records]
                present = [(e, m) for e, m in present if _units(m, family)]
                if present:
                    row_cells.append(
                        model.delivery_row(
                            scope=scope_id,
                            evidence_class=cls,
                            status="supported_with_caveat",
                            evidence_refs=[e["evaluation_id"] for e, _ in present],
                            rationale="comparison units present in scope",
                        )
                    )
                    for e, m in present:
                        short = e["evaluation_id"].split(".")[-1]
                        if family == "stability":
                            key_metrics[f"call_retention.{short}"] = m["call_retention"]
                            key_metrics[f"genotype_concordance.{short}"] = m["genotype_concordance"]
                        else:
                            key_metrics["caller_site_agreement_rate"] = m["site_agreement_rate"]
                            key_metrics["caller_genotype_concordance"] = m["genotype_concordance"]
                            key_metrics["caller_shared_units"] = m["shared_variants"]
                else:
                    row_cells.append(
                        model.delivery_row(
                            scope=scope_id,
                            evidence_class=cls,
                            status="not_evaluated",
                            evidence_refs=[],
                            rationale="no comparison units in scope"
                            if records
                            else "evaluation not run",
                        )
                    )
            cross_metrics = None
            if cross and variant_type == "snp":
                by = cross[0]["metrics"]["by_stratum"]
                cross_metrics = (
                    cross[0]["metrics"] if desc_stratum is None else by.get(desc_stratum)
                )
            if cross_metrics and cross_metrics.get("variant_records"):
                row_cells.append(
                    model.delivery_row(
                        scope=scope_id,
                        evidence_class="reference_sample_self_consistency",
                        status="supported_with_caveat",
                        evidence_refs=[cross[0]["evaluation_id"]],
                        rationale="Illumina calls of the reference BioSample against its own PacBio-derived assembly",
                    )
                )
                key_metrics["reference_sample_illumina_non_reference_records_per_mb"] = (
                    cross_metrics["variant_records_per_mb"]
                )
                key_metrics["reference_sample_illumina_heterozygous_fraction"] = cross_metrics[
                    "heterozygous_fraction_of_non_reference_calls"
                ]
            else:
                reason = (
                    "cross-platform calls were stratified by region, not by variant type or dosage; they are reported on the snp row"
                    if cross and variant_type != "snp"
                    else ("no calls in scope" if cross else "evaluation not run")
                )
                row_cells.append(
                    model.delivery_row(
                        scope=scope_id,
                        evidence_class="reference_sample_self_consistency",
                        status="not_evaluated",
                        evidence_refs=[],
                        rationale=reason,
                    )
                )
            desc_metrics = None
            if descriptive and variant_type == "snp":
                desc_metrics = (
                    descriptive["metrics"]
                    if desc_stratum is None
                    else descriptive["metrics"]["by_stratum"].get(desc_stratum)
                )
            if desc_metrics and desc_metrics.get("variant_records"):
                row_cells.append(
                    model.delivery_row(
                        scope=scope_id,
                        evidence_class="descriptive_stratification",
                        status="supported_with_caveat",
                        evidence_refs=[descriptive["evaluation_id"]],
                        rationale="GS SNP PASS panel counts in scope",
                    )
                )
                for name in (
                    "bases",
                    "variant_records",
                    "variant_records_per_mb",
                    "missing_genotype_fraction",
                    "masked_genotype_fraction",
                    "heterozygous_fraction_of_non_reference_calls",
                ):
                    key_metrics[f"gs_panel_{name}"] = desc_metrics.get(name)
            else:
                reason = {
                    "indel": "descriptive.gs_panel.indel",
                    "het": "GS panel counts are per record, reported on the snp row (heterozygous_fraction_of_non_reference_calls)",
                    "hom_alt": "GS panel counts are per record, reported on the snp row",
                }.get(variant_type, "no GS panel records in scope")
                row_cells.append(
                    model.delivery_row(
                        scope=scope_id,
                        evidence_class="descriptive_stratification",
                        status="not_evaluated",
                        evidence_refs=[],
                        rationale=reason,
                    )
                )

            evaluated = [c for c in row_cells if c["status"] == "supported_with_caveat"]
            if evaluated:
                classes = sorted({c["evidence_class"] for c in evaluated})
                overall = model.delivery_row(
                    scope=scope_id,
                    evidence_class=classes[0],
                    status="supported_with_caveat",
                    evidence_refs=sorted({r for c in evaluated for r in c["evidence"]}),
                    rationale="evaluated by "
                    + ", ".join(classes)
                    + "; no independent truth, so no accuracy claim",
                )
            else:
                overall = model.delivery_row(
                    scope=scope_id,
                    evidence_class=None,
                    status="not_evaluated",
                    evidence_refs=[],
                    rationale="no evaluation has units in this scope",
                )
            rows.append(
                {
                    "scope": scope_id,
                    "region": WINDOW_NAME if stratum is None else f"{stratum} within {WINDOW_NAME}",
                    "split": variant_type,
                    "overall_status": overall["status"],
                    "claim_allowed": overall["claim_allowed"],
                    "overall_rationale": overall["rationale"],
                    "evidence_classes_evaluated": sorted({c["evidence_class"] for c in evaluated}),
                    "fails_benchmark_callable_rule": scope in CALLABLE_RULE_CAVEAT,
                    "caveats": [CALLABLE_RULE_CAVEAT[scope]]
                    if scope in CALLABLE_RULE_CAVEAT
                    else [],
                    "key_metrics": key_metrics,
                }
            )
            cells.extend(row_cells)
    rows.append(
        {
            "scope": "genome_outside_window:snp+indel",
            "region": "GCF_016808095.1 outside " + WINDOW_NAME,
            "split": "snp+indel",
            "overall_status": "not_evaluated",
            "claim_allowed": "no claim; not evaluated",
            "overall_rationale": "callable.outside_window; downsampling.other_samples_and_genome (genome-wide GS panel descriptive counts exist in gs_snp_pass_51.genome_wide.descriptive)",
            "evidence_classes_evaluated": [],
            "fails_benchmark_callable_rule": False,
            "caveats": [],
            "key_metrics": {},
        }
    )
    return {"rule": __doc__.split("Matrix rule")[1].strip(), "rows": rows, "cells": cells}


def matrix_tsv(document: dict[str, object]) -> str:
    columns = [
        "scope",
        "split",
        "overall_status",
        "evidence_classes_evaluated",
        "fails_benchmark_callable_rule",
        "call_retention.f050_s65",
        "call_retention.f050_s66",
        "call_retention.f025_s65",
        "call_retention.f025_s66",
        "caller_site_agreement_rate",
        "caller_genotype_concordance",
        "reference_sample_illumina_non_reference_records_per_mb",
        "reference_sample_illumina_heterozygous_fraction",
        "gs_panel_bases",
        "gs_panel_variant_records_per_mb",
        "gs_panel_missing_genotype_fraction",
        "gs_panel_masked_genotype_fraction",
        "gs_panel_heterozygous_fraction_of_non_reference_calls",
    ]
    lines = ["\t".join(columns)]
    for row in document["rows"]:
        values = []
        for column in columns:
            if column in row:
                value = row[column]
            else:
                value = row["key_metrics"].get(column)
            if isinstance(value, list):
                value = ",".join(value)
            values.append("NA" if value is None else str(value))
        lines.append("\t".join(values))
    return "\n".join(lines) + "\n"


def assemble(
    workdir: Path, reference_fasta: Path, git_sha: str, reference_sha: str
) -> dict[str, object]:
    ref_desc = _ref_description(reference_fasta, reference_sha)
    evaluations = comparison_evaluations(workdir, reference_fasta, reference_sha)
    for path in (
        workdir / "stratify" / "genome_wide.json",
        workdir / "stratify" / "window.json",
        workdir / "crossplatform" / "stratify.json",
    ):
        if path.exists():
            record = json.loads(path.read_text())
            model.validate_evaluation(record)
            evaluations.append(record)
    evaluations.extend(not_evaluated_records(ref_desc))
    for record in evaluations:
        for side in ("query_dataset", "comparator_dataset"):
            if isinstance(record.get(side), dict) and "path" in record[side]:
                record[side]["path"] = str(record[side]["path"]).replace(str(workdir) + "/", "")
        model.validate_evaluation(record)
    out = workdir / "evidence"
    out.mkdir(exist_ok=True)
    (out / "evaluations.json").write_text(json.dumps(evaluations, indent=2) + "\n")
    document = matrix(evaluations)
    (out / "delivery_support_matrix.json").write_text(json.dumps(document, indent=2) + "\n")
    (out / "delivery_support_matrix.tsv").write_text(matrix_tsv(document))
    inputs = {}
    for name in (
        "regions/region_assets.json",
        "mappability/mappability.json",
        "callable/callable_summary.json",
        "derived/derived_strata.json",
        "stratify/inputs.sha256",
        "stratify/outputs.sha256",
        "downsampling/inputs.sha256",
        "downsampling/outputs.sha256",
        "crossplatform/partial.sha256",
        "crossplatform/fastq_used.sha256",
        "crossplatform/outputs.sha256",
    ):
        if (workdir / name).exists():
            inputs[name] = sha256(workdir / name)
    manifest = {
        "evaluation_schema_version": model.EVALUATION_SCHEMA_VERSION,
        "git_sha": git_sha,
        "reference": ref_desc,
        "window": WINDOW_NAME,
        "containers": CONTAINERS,
        "downsampling": {"sample": DOWNSAMPLING_SAMPLE, "subsamples": DOWNSAMPLES},
        "input_manifests_sha256": inputs,
        "outputs_sha256": {
            name: sha256(out / name)
            for name in (
                "evaluations.json",
                "delivery_support_matrix.json",
                "delivery_support_matrix.tsv",
            )
        },
    }
    (out / "evidence_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("derive-strata", "assemble"))
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-sha256")
    parser.add_argument("--git-sha")
    args = parser.parse_args(argv)
    if args.command == "derive-strata":
        derive_strata(args.workdir, args.reference_fasta)
    else:
        if not args.git_sha or not args.reference_sha256:
            parser.error("assemble needs --git-sha and --reference-sha256")
        assemble(args.workdir, args.reference_fasta, args.git_sha, args.reference_sha256)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
