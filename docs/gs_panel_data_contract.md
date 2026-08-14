# GS panel data contract

This document defines the data contract for the genomic selection (GS) SNP
panel produced under `gs_panel/`. It exists so that a downstream consumer —
today, potentially [`genomic-prediction-resnet-hybrid`](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid)
— can read the panel without guessing at file formats, encodings, or what a
given value does and does not mean.

## Scope and what this is not

GS is a conditional, higher-level service on top of the variant-calling
pipeline: it only produces a usable panel once per-individual genotype data
has been called, filtered, and normalized. This document, and the panel
itself, make no claim that:

- the configured hard-filter thresholds are appropriate for a real adzuki
  bean cohort (see the main [README](../README.md#variant-calling-and-filtering-parameters)
  and Issue #13's scope notes — real-data threshold validation is tracked
  separately, under Issue #1's `#I`);
- the panel is ready for GS model training without further QC (MAF
  filtering, call-rate filtering, LD pruning, and imputation are all
  explicitly out of scope for this pipeline);
- `genomic-prediction-resnet-hybrid` can currently read this panel. As of
  this writing that repository has no adzuki- or VCF-specific ingestion
  code at all; its only verified data path reads SoyNAM soybean genotype
  data in an unrelated gzip-TSV format (`soynam_data.py`). This document's
  conventions were chosen to match that repository's own tested,
  already-working conventions (dosage encoding, dtype, on-disk shape,
  manifest structure) specifically so that a *future* adzuki loader in
  that repository has a well-specified contract to build against — not
  because a matching loader already exists. See "Downstream compatibility"
  below.

## Input contract

The GS panel lineage starts from `variants/raw/cohort.raw.vcf.gz`
(raw/all), not from any of the type-selected, filtered, or PASS VCFs the
primary lineage produces. This is a deliberate, verified choice, not an
oversight: `GATK_SELECTVARIANTS` classifies each `VariantContext` as a
whole into exactly one overall type (`SNP`, `INDEL`, `MIXED`, `MNP`, ...)
and selects it only on an exact match (confirmed against the pinned GATK
4.6.2.0 container in `tests/modules/gatk_selectvariants.nf.test`). A
`MIXED`-type record — one site with both a SNP-type and an indel-type ALT
allele — is therefore excluded from *both* `cohort.snp.vcf.gz` and
`cohort.indel.vcf.gz`, and by extension from every VCF derived from them.
raw/all is the only stage that can still contain such a record, which is
why it is the only valid starting point for a normalization step that is
actually supposed to split and reclassify one.

`FILTER=PASS` in the GS lineage's own PASS VCF
(`variants/gs_pass/cohort_gs.snp.pass.vcf.gz`) means exactly what it means
in the primary lineage: the record passed the configured SNP hard filters
(`snp_filter_*` parameters). It is not a claim that the record is
biologically validated, and it is not a claim that the hard-filter
thresholds themselves are well-calibrated for adzuki bean.

Sample order is whatever order the samples appear in the VCF's `#CHROM`
header row, which is stable across every stage of this pipeline (no
process reorders or drops sample columns; only records are added, split,
or removed). Variant order is genomic-coordinate order, inherited from
`GATK_GatherVcfs`'s reference-contig ordering.

Reference FASTA consistency is enforced by `bcftools norm`'s own
`--check-ref e` behavior: a record whose declared REF does not match the
actual reference base at that position is a fatal error (confirmed
directly: `tests/modules/gs_normalize_variants.nf.test` includes a
dedicated case with a deliberately wrong REF, and the pinned bcftools 1.24
container correctly refuses to normalize it).

## Normalization (concern 1)

`GS_NORMALIZE_VARIANTS` runs `bcftools norm --fasta-ref <reference>
--multiallelics -both --check-ref e` (pinned
`quay.io/biocontainers/bcftools:1.24--h118bc1c_2`) directly on raw/all.
The following behaviors were confirmed empirically against the pinned
container, using a hand-built fixture containing a pure SNP, a pure
multiallelic SNP, a pure indel, and a genuine MIXED record, before this
module was written — not assumed:

- A MIXED record splits into one row per original ALT allele, each now
  trivially biallelic. For example, `REF=C ALT=T,CAAA` becomes two rows:
  `REF=C ALT=T` and `REF=C ALT=CAAA`, at the same position.
- Each split row's genotypes are **recoded relative to that row's own
  allele** (`0`/`1` semantics per split), not the original multiallelic
  allele index. Allele order is not canonicalized (`1/0` can appear, not
  only `0/1`); any genotype parser reading normalized output must treat
  the two alleles as an unordered pair.
- INFO fields declared `Number=A` (one value per ALT allele, e.g. `AC`)
  are correctly de-multiplexed to a single scalar value per split row.
  INFO fields declared `Number=1` (site-level — this includes every GATK
  hard-filter annotation: `QD`, `SOR`, `FS`, `MQ`, `MQRankSum`,
  `ReadPosRankSum`) are **duplicated verbatim, unrecomputed**, to every
  split row. A split row's quality annotations describe the *original,
  possibly multi-type* site, not a value freshly computed for that row's
  allele alone. This is a known, accepted limitation (see below), not
  something this pipeline corrects.
- FILTER is likewise duplicated verbatim to every split row regardless of
  that row's new shape, which can leave a semantically mismatched tag
  (e.g. a SNP-named filter surviving on a row now shaped like an indel).
- An empty (header-only) input VCF completes normally (exit 0, valid
  empty output) — no special-casing needed for this pipeline's default
  synthetic-fixture case where the input can legitimately be small.

`CLASSIFY_NORMALIZED_VARIANTS` (`bin/classify_normalized_variants.py`)
then reclassifies every post-split row from scratch, by its own REF/ALT
shape — **not** by trusting any pre-split GATK type label, which no
longer applies once a site has been split:

| Post-split shape | Class | Eligible for the panel? |
| --- | --- | --- |
| `len(REF) == 1` and `len(ALT) == 1` | `snp` | yes, subject to the duplicate-key check below |
| `len(REF) == len(ALT) > 1` | `mnp` | no — excluded, counted |
| `len(REF) != len(ALT)` (and not symbolic) | `indel` | no — excluded, counted |
| ALT contains `<`, `[`, `]`, or is `*` | `symbolic_or_star` | no — excluded, counted |

A row whose ALT still contains a comma after `bcftools norm -m-` is a
hard error (`MalformedVcfError`): this should never happen post-split, so
silently reclassifying it as some shape would hide a real anomaly.

**Duplicate `(CHROM, POS, REF, ALT)` keys** among `snp`-classified rows
are excluded **in full** — every occurrence of a colliding key, not just
the extras — because there is no automatic way to decide which of two
identical records is "correct"; keeping one at random would be a silent,
unreviewable choice. The count and the actual colliding keys are reported
in `cohort_gs.classification_accounting.summary.txt`.

Every kept record's FILTER is reset to `.` before this stage's output is
handed to hard filtering (see below): this is currently a no-op in
practice, since raw/all's FILTER column is always `.` already (GATK never
runs `VariantFiltration` before this point) — the reset only matters if a
future change feeds this step a VCF that has already been filtered
elsewhere. It is kept anyway as a defensive, explicit statement of intent
rather than relying on the current pipeline's specific shape.

The Python-classified-but-still-plain-text VCF is compressed and indexed
by a separate `GS_INDEX_CLASSIFIED_VARIANTS` step (`bcftools view
--output-type z` + `bcftools index --tbi`) rather than by
`classify_normalized_variants.py` itself: the `python:3.12` image has
neither `bgzip` nor `tabix` (confirmed directly), and the standard
library's `gzip` module writes plain DEFLATE, which `tabix` cannot index.
This mirrors the existing `BCFTOOLS_STATS` → `SUMMARIZE_VARIANT_QC`
two-container split established for Issue #16.

## GS panel eligibility (concern 2)

The classified, biallelic-SNP-only VCF is re-filtered using the **same**
SNP hard-filter expressions and thresholds as the primary lineage
(`snpHardFilters()` in `workflows/adzuki_snp_pipeline.nf`, shared between
both call sites so the thresholds are never duplicated), via a second,
distinctly-named invocation of the existing `GATK_VARIANTFILTRATION` and
`GATK_SELECTPASSVARIANTS` modules (`GATK_VARIANTFILTRATION_GS` /
`GATK_SELECTPASSVARIANTS_GS`). No new filtering criteria (MAF, call rate,
LD) are introduced at this stage — that is explicitly out of scope.

Because GATK's hard-filter expressions read only `Number=1` INFO
annotations, and those are the *unmodified, duplicated* values described
above, a record that was already a pure biallelic SNP before
normalization receives **exactly the same PASS/FAIL verdict** through the
GS lineage as through the primary lineage — verified directly against the
shared synthetic fixture (both lineages tag the same two records
`SNP_SOR_HIGH` with identical `SOR` values, and both lineages' PASS VCFs
end up empty). For a genuinely MIXED site's SNP-shaped child, this
re-filtering step is the *only* time that data is ever evaluated against
SNP thresholds at all — GATK itself never evaluates a MIXED site against
pure-SNP filter expressions.

**Known limitation:** a MIXED-derived row's hard-filter verdict is based
on annotations computed at the original, combined (MIXED) site, not
recomputed from the split allele alone. Recomputing annotations
per-split-allele would require re-running variant calling and genotyping
against the decomposed representation, which is out of scope here.

## Genotype matrix contract (concern 3)

`BUILD_GS_PANEL` (`bin/build_gs_panel.py`) parses the GS-eligible PASS VCF
directly with the standard library's `gzip` module — not by re-invoking
`bcftools` — specifically so that an independent `bcftools query`-based
cross-check of the same file remains a meaningful verification, matching
this repository's established convention for its other QC scripts.

### Dosage encoding

| Genotype | Dosage |
| --- | ---: |
| `0/0` (homozygous reference) | `-1` |
| `0/1` or `1/0` (heterozygous) | `0` |
| `1/1` (homozygous alternate) | `+1` |
| missing or non-standard (see below) | `nan` |

This is the **same** `-1`/`0`/`+1` scheme `genomic-prediction-resnet-hybrid`'s
own `soynam_data.py` already uses for its (unrelated) SoyNAM data
(`GENOTYPE_ENCODING`), confirmed algebraically equivalent to a standard
`0`/`1`/`2` additive allele-count dosage shifted by `-1`
(`gblup_baseline.py` recovers allele frequency as `(mean + 1) / 2`). This
was a deliberate choice to match that repository's tested, working
convention, not an arbitrary default — the two repositories' encodings
are directly comparable without a translation step, if a future adzuki
loader is ever written there.

`int8` (or any sentinel-integer missing value) was explicitly **not**
used, for two reasons: an integer type cannot represent a missing value
as `NaN` without a sentinel that could collide with a real dosage value
or force every consumer to special-case it, and the sibling repository's
own downstream computation (mean imputation, VanRaden relationship
matrices, PCA, standardization) operates in floating point throughout. A
missing cell is always the literal token `nan`; `-1`, `0`, `1`, and `nan`
are the *only* four tokens ever written, and all four round-trip through
`float()` without special-casing.

### Genotype classification (missing and non-standard handling)

A genotype is encoded as a dosage only if it is a clean, unphased,
diploid, biallelic-index call (`0/0`, `0/1`, `1/0`, or `1/1`). Every other
shape is treated as missing in the matrix (`nan`), but counted under its
own specific reason — never folded into one undifferentiated "missing"
bucket:

| Reason | Example | Counted as |
| --- | --- | --- |
| Missing | `.`, `./.`, `0/.` (any allele position is `.`) | `missing_calls` |
| Phased | `0\|1` | `phased_calls_treated_as_missing` |
| Non-diploid | `0` (haploid), `0/0/1` (triploid) | `non_diploid_calls_treated_as_missing` |
| Non-biallelic allele index | `0/2` (defensive — should not occur given the input is already biallelic-only) | `non_biallelic_index_calls_treated_as_missing` |

`cohort.gs_panel.genotype_encoding_accounting.tsv` reports the cohort-wide
total for each reason (plus `total_treated_as_missing`, the sum of all
four); the sample and variant metadata files (below) report the same
breakdown at finer granularity.

### On-disk shape vs. in-memory shape

`cohort.gs_panel.genotype_matrix.tsv.gz` is a gzip-compressed TSV with
**variant rows and sample columns**: the header row is
`variant_key`, followed by one column per sample (in `#CHROM` header
order); each data row is one variant's `variant_key` (`CHROM:POS:REF:ALT`,
computed from the post-split, post-classification record) followed by
one dosage cell per sample.

This on-disk shape matches `genomic-prediction-resnet-hybrid`'s own
`_load_genotype_frame` (`soynam_data.py`): that function reads a
marker-rows-by-sample-columns TSV and only transposes it to
sample-rows-by-marker-columns *after* loading, to build the in-memory
`SoynamDataset.genotypes` array. A future adzuki loader is expected to do
the same transpose after reading this file, rather than expecting
sample-rows-by-variant-columns on disk.

A gzipped plain-text TSV was chosen over a binary format (e.g. NumPy's
`.npz`) for two reasons: this repository's `bin/*.py` scripts are
standard-library-only by established convention, with no dependency-
pinning infrastructure (no `requirements.txt`, `pyproject.toml`, or
`Dockerfile` for Python dependencies exists anywhere in the repository);
and `.npz` would not even match the sibling repository's own real
precedent — its actual, verified genotype ingestion path reads gzip TSV,
not a binary array format. Panel scales here (thousands of markers,
tens-to-low-hundreds of samples) are far below where TSV parsing becomes
a practical bottleneck.

The matrix header row is always written, even when there are zero
variant rows (see "Empty panel contract" below) — an empty panel must
never lose the sample list.

## Metadata (concern 4)

`cohort.gs_panel.sample_metadata.tsv` — one row per sample, in `#CHROM`
header order:

| Column | Meaning |
| --- | --- |
| `cohort_id` | Cohort identifier |
| `sample_index` | 0-indexed column position in the GS-eligible PASS VCF header (and in the matrix) |
| `sample_id` | Sample ID |
| `missing_genotype_count` | Count of `nan` cells for this sample, any reason |
| `missing_genotype_rate` | `missing_genotype_count / total_variants`, or `NA` if there are zero variants |
| `non_standard_genotype_count` | Subset of `missing_genotype_count` caused by phased/non-diploid/non-biallelic-index calls specifically (excludes plain missing) |

`cohort.gs_panel.variant_metadata.tsv` — one row per variant, in
genomic-coordinate (file) order:

| Column | Meaning |
| --- | --- |
| `cohort_id` | Cohort identifier |
| `variant_index` | 0-indexed row position in the matrix |
| `variant_key` | `CHROM:POS:REF:ALT` |
| `chrom`, `pos`, `ref`, `alt`, `qual` | The post-split, post-classification record's own fixed VCF columns |
| `missing_genotype_count` | Count of `nan` cells for this variant, any reason |
| `missing_genotype_rate` | `missing_genotype_count / total_samples` |

## Record accounting (concern 5)

`cohort.gs_panel.record_accounting.tsv` (`bin/reconcile_gs_panel_accounting.py`)
reconciles the full lineage:

```
raw_all_records
  -> normalized_records            (post bcftools norm -m- splitting)
  -> classified_biallelic_snp_records   (post shape reclassification and duplicate-key exclusion)
  -> gs_pass_records                    (post GS-specific hard filtering)
  -> final_matrix_variant_records       (should equal gs_pass_records)
```

`raw_all_records` and `normalized_records` are independently re-counted
directly from their VCFs (not trusted from any other script's output),
matching `bin/reconcile_variant_type_counts.py`'s established
"independently re-verify" convention. `classified_biallelic_snp_records`
is read from `classify_normalized_variants.py`'s own accounting rather
than re-derived, since re-implementing its shape-classification and
duplicate-key logic a second time would duplicate that logic rather than
verify it.

`gs_hard_filter_excluded_records = classified_biallelic_snp_records -
gs_pass_records` is expected to be `>= 0` (hard filtering only removes
records); a negative value, or `final_matrix_variant_records` not
matching `gs_pass_records`, is reported as a `WARNING` in the summary
text rather than hidden — the same "never silently hide a surprising
value" convention as `bin/reconcile_variant_type_counts.py`.

## Reproducibility manifest (concern 6)

`cohort.gs_panel.manifest.json` (`bin/build_gs_panel_manifest.py`) is
schema-versioned (`schema_version: 1`) and mirrors the *shape* of the
sibling repository's own `run_manifest.py` — a sortable `run_id`
(`<UTC-timestamp>-<uuid4-hex8>`), deterministic canonical JSON, a
self-referential `manifest_hash`, and filename-only checksums (never an
absolute path) — without depending on that repository's code: a
Nextflow bioinformatics pipeline has no reason to import a model-training
repository's Python package, so the pattern is reimplemented in
stdlib-only Python instead.

```json
{
  "schema_version": 1,
  "run_id": "20260814T074225Z-c22d6e72",
  "generated_at": "2026-08-14T07:42:25Z",
  "cohort_id": "cohort",
  "pipeline_version": "0.2.0-dev",
  "git_commit": null,
  "containers": { "bcftools": "...", "gatk": "...", "python": "..." },
  "parameters": { "snp_filter_qd_min": 2.0, "...": "..." },
  "panel_status": "empty",
  "checksums": { "cohort.gs_panel.genotype_matrix.tsv.gz": "sha256:...", "...": "..." },
  "manifest_hash": "sha256:..."
}
```

Software versions are recorded as the pinned container image references
already baked into each process's `container` directive — the ground
truth for "what actually ran" in a Nextflow-plus-Docker pipeline — rather
than by shelling out to e.g. `bcftools --version` at run time.

`git_commit` is Nextflow's own `workflow.commitId`, resolved by the
workflow and passed in as a plain argument rather than computed inside a
task container (a task's work directory holds only staged input files,
never the pipeline's `.git` checkout, so a container-side `git
rev-parse` would find no repository at all). `workflow.commitId` is
documented as populated only when Nextflow itself pulls a git-hosted
pipeline (`nextflow run owner/repo`); it is legitimately `null` for this
repository's own documented `nextflow run .` local-directory invocation,
and is reported as such honestly rather than worked around.

## Empty panel contract

Zero GS-eligible PASS records is a normal, expected outcome — not an
error — and is in fact what the pipeline's own default synthetic fixture
produces today (both synthetic SNPs fail `SNP_SOR_HIGH`, matching the
primary lineage's already-documented behavior). In this case:

- The matrix header (sample list) is still written; there are simply
  zero data rows.
- Sample metadata still lists every sample, with `missing_genotype_rate`
  reported as `NA` (not `0`, and not a division-by-zero crash) since the
  rate is undefined when there are zero variants to measure it against.
- `cohort.gs_panel.record_accounting.tsv` reports `panel_status: empty`,
  a machine-readable status any downstream reader can branch on, plus a
  human-readable explanation in the accompanying `.summary.txt` that
  this is not a failure.

A **zero-sample** header, by contrast, is treated as a hard error: this
pipeline always has at least one sample by construction, so a VCF with no
sample columns at all indicates a genuine anomaly upstream, not a valid
(if unusual) GS outcome.

## Downstream compatibility

| This pipeline | `genomic-prediction-resnet-hybrid` (SoyNAM path) |
| --- | --- |
| Dosage: `-1`/`0`/`+1`, `nan` for missing | Same (`GENOTYPE_ENCODING`) |
| dtype: text tokens, float64 semantics | float64 (`FloatArray`) throughout |
| On-disk: variant rows x sample columns (gzip TSV) | Marker rows x sample columns (gzip TSV), transposed after load |
| Sample/variant metadata: separate files from the matrix | `SoynamDataset` separates `genotypes` / `family_ids`+`sample_names` / `marker_names` |
| Manifest: schema-versioned JSON, run_id, checksums, atomic write | `run_manifest.py`: same shape, for GBLUP/ResNet run artifacts |

As stated above, no code in `genomic-prediction-resnet-hybrid` currently
reads this panel — the table above documents convention alignment, not an
existing integration. Adding an adzuki/VCF-panel loader to that
repository is out of scope for this pipeline and is expected to be
tracked as a separate issue on that repository.
