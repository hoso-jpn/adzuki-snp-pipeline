#!/usr/bin/env python3
"""Build a genomic selection (GS) SNP panel from the GS-eligible PASS VCF.

Reads `cohort_gs.snp.pass.vcf.gz` directly (no bcftools re-invocation, to
keep an independent `bcftools query` cross-check meaningful) and emits,
from a single pass over the file:

- a genotype dosage matrix (variant rows x sample columns, matching the
  on-disk convention of the sibling `genomic-prediction-resnet-hybrid`
  repository's own SoyNAM genotype loader, which reads a marker-rows x
  sample-columns TSV and transposes it after loading -- see
  `docs/gs_panel_data_contract.md` for the full reasoning);
- sample metadata (concern 4, separate from the matrix itself);
- variant metadata (concern 4, separate from the matrix itself);
- a genotype-encoding accounting report, broken down by *why* a cell
  became missing (concern 5's genotype-level counterpart -- never
  folding every non-standard call into one undifferentiated bucket).

Dosage encoding is -1 / 0 / +1 for homozygous-reference / heterozygous /
homozygous-alternate, matching the sibling repository's own
`GENOTYPE_ENCODING` convention (`soynam_data.py`): a real SNP/indel
dosage is `allele_count - 1`, so `0/0` -> -1, `0/1` (or `1/0`) -> 0, and
`1/1` -> +1. Missing is IEEE754 NaN, never a sentinel integer, because
a sentinel could collide with a real dosage value or require every
downstream consumer to special-case it.

Phasing (`|` vs `/`) does not affect dosage: per the VCF specification,
the separator only records whether the call is phased, not which or how
many alleles are present, so `0|1` carries exactly the same allele
count as `0/1` and must resolve to the same dosage (0). Earlier
revisions of this script treated any phased call as missing, which was
wrong -- phase is orthogonal to additive dosage. Every genotype is
still checked for how many phased calls it contained
(`phased_genotype_count` in the accounting output), but that count is
informational and never removes a cell from the dosage matrix by
itself.

A genotype is encoded as a dosage as long as it is diploid with a
biallelic-index call (`0/0`, `0/1`, `1/0`, `1/1`, in either phasing).
Every other shape -- missing (`.`, `./.`, or any allele position that
is `.`), non-diploid (haploid, triploid, ...), or an allele index
outside `{0, 1}` (defensive: the input is already biallelic-only by
construction) -- is treated as missing in the matrix, but counted
under its own specific reason so "never silently coerce" is checkable
with real numbers, not just asserted in prose.

This encoding is diploid-only by design (schema v1): `--sample-ploidy`
must equal 2, checked before any other work, because a non-diploid
ploidy would make every genotype call "non-diploid-shaped" by
definition, silently producing an all-missing (but successfully
completing) panel rather than a meaningful error.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

DOSAGE_BY_ALT_COUNT: dict[int, str] = {0: "-1", 1: "0", 2: "1"}
MISSING_CELL_TOKEN = "nan"
NOT_APPLICABLE = "NA"

GENOTYPE_ACCOUNTING_HEADER: tuple[str, ...] = ("cohort_id", "metric", "value")
SAMPLE_METADATA_HEADER: tuple[str, ...] = (
    "cohort_id",
    "sample_index",
    "sample_id",
    "missing_genotype_count",
    "missing_genotype_rate",
    "non_standard_genotype_count",
)
VARIANT_METADATA_HEADER: tuple[str, ...] = (
    "cohort_id",
    "variant_index",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "missing_genotype_count",
    "missing_genotype_rate",
)


class MalformedVcfError(Exception):
    """Raised when a GS-eligible PASS VCF cannot be turned into a panel safely."""


@dataclass(frozen=True)
class GsPassRecord:
    """The subset of one GS-eligible PASS VCF row this tool needs."""

    chrom: str
    pos: str
    ref: str
    alt: str
    qual: str
    sample_genotypes: tuple[str, ...]

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


@dataclass(frozen=True)
class GsPassVcf:
    """Every sample name and record read from a GS-eligible PASS VCF."""

    sample_names: tuple[str, ...]
    records: tuple[GsPassRecord, ...]


@dataclass(frozen=True)
class GenotypeCell:
    """One sample's classified genotype at one variant.

    ``is_phased`` is orthogonal to ``category``/``dosage``: a phased
    call that is otherwise a clean diploid biallelic-index genotype is
    ``category="standard"`` with a real dosage, exactly like its
    unphased counterpart -- phasing is tracked for informational
    accounting only, never as a reason to treat a cell as missing.
    """

    category: str
    dosage: str
    is_phased: bool


def _locate_gt_index(format_field: str, path: Path) -> int:
    keys = format_field.split(":")
    if "GT" not in keys:
        raise MalformedVcfError(f"{path}: FORMAT field has no GT subfield: {format_field}")
    return keys.index("GT")


def parse_gs_pass_vcf(path: Path) -> GsPassVcf:
    """Parse the CHROM/POS/REF/ALT/QUAL/GT columns of a bgzipped VCF."""
    sample_names: tuple[str, ...] | None = None
    records: list[GsPassRecord] = []

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) <= 9:
                    raise MalformedVcfError(
                        f"{path}: #CHROM header has {len(fields)} fields, expected "
                        "at least 10 (9 fixed columns plus one or more samples)"
                    )
                sample_names = tuple(fields[9:])
                continue

            if line.startswith("#"):
                continue

            if sample_names is None:
                raise MalformedVcfError(f"{path}: data row seen before #CHROM header")

            fields = line.split("\t")
            if len(fields) < 10:
                raise MalformedVcfError(
                    f"{path}: line {line_number}: data row has {len(fields)} "
                    "tab-separated fields, expected at least 10"
                )

            chrom, pos, _id, ref, alt, qual = fields[:6]
            format_field = fields[8]
            sample_fields = fields[9:]

            gt_index = _locate_gt_index(format_field, path)
            genotypes = tuple(
                sample_field.split(":")[gt_index] for sample_field in sample_fields
            )

            records.append(
                GsPassRecord(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    sample_genotypes=genotypes,
                )
            )

    if sample_names is None:
        raise MalformedVcfError(f"{path}: no #CHROM header line found")

    return GsPassVcf(sample_names=sample_names, records=tuple(records))


def classify_genotype(gt: str) -> GenotypeCell:
    """Classify one raw GT string into a category and its matrix dosage token.

    Phasing (``|`` vs ``/``) never changes the category or dosage: it
    only changes ``is_phased``. The VCF specification defines ``|``/``/``
    as recording phase, not allele identity or count, so a phased call
    is resolved exactly like its unphased counterpart.
    """
    is_phased = "|" in gt
    alleles = gt.split("|") if is_phased else gt.split("/")

    if any(allele in (".", "") for allele in alleles):
        return GenotypeCell(category="missing", dosage=MISSING_CELL_TOKEN, is_phased=is_phased)

    if len(alleles) != 2:
        return GenotypeCell(category="non_diploid", dosage=MISSING_CELL_TOKEN, is_phased=is_phased)

    if any(allele not in ("0", "1") for allele in alleles):
        return GenotypeCell(
            category="non_biallelic_index", dosage=MISSING_CELL_TOKEN, is_phased=is_phased
        )

    alt_count = alleles.count("1")
    return GenotypeCell(
        category="standard", dosage=DOSAGE_BY_ALT_COUNT[alt_count], is_phased=is_phased
    )


def _format_rate(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return NOT_APPLICABLE

    return f"{numerator / denominator:.6f}"


def build_matrix_rows(vcf: GsPassVcf) -> list[list[str]]:
    """Build the genotype matrix's data rows: one row per variant."""
    rows: list[list[str]] = []

    for record in vcf.records:
        cells = [classify_genotype(gt) for gt in record.sample_genotypes]
        rows.append([record.variant_key, *(cell.dosage for cell in cells)])

    return rows


def build_sample_metadata_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the sample metadata rows: one row per sample, in header order."""
    total_variants = len(vcf.records)
    missing_counts = [0] * len(vcf.sample_names)
    non_standard_counts = [0] * len(vcf.sample_names)

    for record in vcf.records:
        for sample_position, gt in enumerate(record.sample_genotypes):
            cell = classify_genotype(gt)
            if cell.category == "missing":
                missing_counts[sample_position] += 1
            elif cell.category != "standard":
                missing_counts[sample_position] += 1
                non_standard_counts[sample_position] += 1

    rows: list[list[str]] = []
    for sample_index, sample_id in enumerate(vcf.sample_names):
        rows.append(
            [
                cohort_id,
                str(sample_index),
                sample_id,
                str(missing_counts[sample_index]),
                _format_rate(missing_counts[sample_index], total_variants),
                str(non_standard_counts[sample_index]),
            ]
        )

    return rows


def build_variant_metadata_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the variant metadata rows: one row per variant, in file order."""
    total_samples = len(vcf.sample_names)
    rows: list[list[str]] = []

    for variant_index, record in enumerate(vcf.records):
        missing_count = sum(
            1
            for gt in record.sample_genotypes
            if classify_genotype(gt).category != "standard"
        )
        rows.append(
            [
                cohort_id,
                str(variant_index),
                record.variant_key,
                record.chrom,
                record.pos,
                record.ref,
                record.alt,
                record.qual,
                str(missing_count),
                _format_rate(missing_count, total_samples),
            ]
        )

    return rows


def build_genotype_accounting_rows(cohort_id: str, vcf: GsPassVcf) -> list[list[str]]:
    """Build the cohort-wide genotype-encoding accounting rows.

    ``phased_genotype_count`` is reported separately from every other
    metric here: it counts calls that were phased, regardless of
    whether they were standard, missing, or otherwise non-standard, and
    is never added into ``total_treated_as_missing`` -- a phased call
    that resolves to a real dosage is not missing.
    """
    counts = {
        "standard_hom_ref_calls": 0,
        "standard_het_calls": 0,
        "standard_hom_alt_calls": 0,
        "missing_calls": 0,
        "non_diploid_calls_treated_as_missing": 0,
        "non_biallelic_index_calls_treated_as_missing": 0,
    }
    dosage_metric_by_token = {
        "-1": "standard_hom_ref_calls",
        "0": "standard_het_calls",
        "1": "standard_hom_alt_calls",
    }
    non_standard_metric_by_category = {
        "missing": "missing_calls",
        "non_diploid": "non_diploid_calls_treated_as_missing",
        "non_biallelic_index": "non_biallelic_index_calls_treated_as_missing",
    }

    total_genotype_cells = 0
    phased_genotype_count = 0
    for record in vcf.records:
        for gt in record.sample_genotypes:
            total_genotype_cells += 1
            cell = classify_genotype(gt)
            if cell.is_phased:
                phased_genotype_count += 1
            if cell.category == "standard":
                counts[dosage_metric_by_token[cell.dosage]] += 1
            else:
                counts[non_standard_metric_by_category[cell.category]] += 1

    total_treated_as_missing = (
        counts["missing_calls"]
        + counts["non_diploid_calls_treated_as_missing"]
        + counts["non_biallelic_index_calls_treated_as_missing"]
    )

    rows = [[cohort_id, "total_genotype_cells", str(total_genotype_cells)]]
    for metric, value in counts.items():
        rows.append([cohort_id, metric, str(value)])
    rows.append([cohort_id, "total_treated_as_missing", str(total_treated_as_missing)])
    rows.append([cohort_id, "phased_genotype_count", str(phased_genotype_count)])

    return rows


def build_genotype_accounting_summary_text(cohort_id: str, vcf: GsPassVcf) -> str:
    """Build the human-readable genotype-encoding accounting summary."""
    accounting = {row[1]: row[2] for row in build_genotype_accounting_rows(cohort_id, vcf)}

    lines = [
        "GS panel genotype encoding summary",
        f"Cohort ID: {cohort_id}",
        f"Variants: {len(vcf.records)}",
        f"Samples: {len(vcf.sample_names)}",
        f"Total genotype cells: {accounting['total_genotype_cells']}",
        f"  standard hom-ref (-1): {accounting['standard_hom_ref_calls']}",
        f"  standard het (0): {accounting['standard_het_calls']}",
        f"  standard hom-alt (+1): {accounting['standard_hom_alt_calls']}",
        f"  missing (nan): {accounting['missing_calls']}",
        (
            "  non-diploid, treated as missing (nan): "
            f"{accounting['non_diploid_calls_treated_as_missing']}"
        ),
        (
            "  non-biallelic-index, treated as missing (nan): "
            f"{accounting['non_biallelic_index_calls_treated_as_missing']}"
        ),
        f"Total cells treated as missing: {accounting['total_treated_as_missing']}",
        (
            f"Phased genotype calls: {accounting['phased_genotype_count']} "
            "(informational only -- phasing does not affect dosage or "
            "missingness; a phased call that is otherwise a clean diploid "
            "biallelic-index genotype is encoded exactly like its unphased "
            "counterpart, per the VCF specification's definition of "
            "'|' as recording phase, not allele identity)."
        ),
        (
            "Every non-standard genotype shape is counted under its own "
            "reason rather than a single undifferentiated 'missing' "
            "bucket; see the sample and variant metadata files for the "
            "same breakdown at finer granularity."
        ),
    ]

    return "\n".join(lines) + "\n"


def write_tsv(path: Path, header: tuple[str, ...], rows: list[list[str]]) -> None:
    """Write a tab-separated file with the given header followed by rows."""
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_matrix(path: Path, vcf: GsPassVcf) -> None:
    """Write the gzipped genotype matrix: variant rows x sample columns.

    The header is always written, even with zero variants, so an empty
    panel never loses the sample list -- only the data rows are absent.

    Compressed with ``mtime=0`` so that identical logical content always
    produces byte-identical compressed output: ``gzip.open`` embeds the
    current wall-clock time in the gzip header by default, which would
    otherwise make two runs over the same input produce different
    checksums in the manifest even though nothing about the data
    changed -- undermining the "same input reproduces the same panel"
    guarantee this contract requires.
    """
    header = ["variant_key", *vcf.sample_names]
    rows = build_matrix_rows(vcf)

    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    text = "\n".join(lines) + "\n"

    path.write_bytes(gzip.compress(text.encode("utf-8"), mtime=0))


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments for the GS panel builder CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a GS genotype matrix, sample/variant metadata, and a "
            "genotype-encoding accounting report from the GS-eligible PASS VCF."
        )
    )
    parser.add_argument(
        "--gs-pass-vcf",
        required=True,
        type=Path,
        help="Path to cohort_gs.snp.pass.vcf.gz.",
    )
    parser.add_argument("--cohort-id", required=True, help="Cohort identifier.")
    parser.add_argument(
        "--sample-ploidy",
        required=True,
        type=int,
        help="The pipeline's configured sample ploidy (params.sample_ploidy); "
        "this schema is diploid-only, so any value other than 2 is a hard error.",
    )
    parser.add_argument(
        "--matrix-output",
        required=True,
        type=Path,
        help="Output path for the gzipped genotype matrix TSV.",
    )
    parser.add_argument(
        "--sample-metadata-output",
        required=True,
        type=Path,
        help="Output path for the sample metadata TSV.",
    )
    parser.add_argument(
        "--variant-metadata-output",
        required=True,
        type=Path,
        help="Output path for the variant metadata TSV.",
    )
    parser.add_argument(
        "--genotype-accounting-output",
        required=True,
        type=Path,
        help="Output path for the genotype-encoding accounting TSV.",
    )
    parser.add_argument(
        "--genotype-accounting-summary-output",
        required=True,
        type=Path,
        help="Output path for the human-readable genotype-encoding summary.",
    )

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI end to end and return a process exit code."""
    args = parse_args(argv)

    if args.sample_ploidy != 2:
        print(
            "build_gs_panel.py: error: this GS panel schema (v1) is diploid-only, "
            f"but --sample-ploidy was {args.sample_ploidy}. Every genotype call "
            "would be classified as non-diploid-shaped and encoded as missing, "
            "which would silently produce an all-missing panel rather than a "
            "meaningful error. A generalized encoding is tracked as future work; "
            "until then, this schema cannot be used for non-diploid cohorts.",
            file=sys.stderr,
        )
        return 1

    try:
        vcf = parse_gs_pass_vcf(args.gs_pass_vcf)
    except OSError as error:
        print(
            f"build_gs_panel.py: error: cannot read {args.gs_pass_vcf}: {error}",
            file=sys.stderr,
        )
        return 1
    except MalformedVcfError as error:
        print(f"build_gs_panel.py: error: {error}", file=sys.stderr)
        return 1

    write_matrix(args.matrix_output, vcf)
    write_tsv(
        args.sample_metadata_output,
        SAMPLE_METADATA_HEADER,
        build_sample_metadata_rows(args.cohort_id, vcf),
    )
    write_tsv(
        args.variant_metadata_output,
        VARIANT_METADATA_HEADER,
        build_variant_metadata_rows(args.cohort_id, vcf),
    )
    write_tsv(
        args.genotype_accounting_output,
        GENOTYPE_ACCOUNTING_HEADER,
        build_genotype_accounting_rows(args.cohort_id, vcf),
    )
    args.genotype_accounting_summary_output.write_text(
        build_genotype_accounting_summary_text(args.cohort_id, vcf),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
