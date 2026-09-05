# GS panel data contract

This document defines the data contract for the genomic selection (GS) SNP
panel produced under `gs_panel/`. It exists so that a downstream consumer —
today, potentially [`genomic-prediction-resnet-hybrid`](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid)
— can read the panel without guessing at file formats, encodings, or what a
given value does and does not mean.

## Scope and what this is not

GS is a conditional, higher-level service on top of the variant-calling
pipeline: it only produces a usable panel once per-individual genotype data
has been called, filtered, and normalized. It is also conditional in a
second, literal sense (Issue #20): the entire GS lineage — normalization
through the reproducibility manifest — only runs at all when
`params.enable_gs_panel` is `true` (the default). When it is `false`, no
process in this lineage starts, and none of `variants/gs_normalized/`,
`variants/gs_classified/`, `variants/gs_filtered/`, `variants/gs_pass/`,
or `gs_panel/` is created; the primary `raw`/`filtered`/`pass`/QC lineage
this document does not otherwise describe is unaffected either way. This
exists specifically so that non-diploid variant calling — which this
schema's diploid-only constraint (see "Diploid-only constraint" below)
would otherwise always fail on — can be exercised end to end without
reaching the GS panel at all. This document, and the panel
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
| ALT is the literal `.` (no alternate allele observed at this site) | `no_alt` | no — excluded, counted |
| `len(REF) == 1` and `len(ALT) == 1` | `snp` | yes, subject to the duplicate-key check below |
| `len(REF) == len(ALT) > 1` | `mnp` | no — excluded, counted |
| `len(REF) != len(ALT)` (and not symbolic) | `indel` | no — excluded, counted |
| ALT contains `<`, `[`, `]`, or is `*` | `symbolic_or_star` | no — excluded, counted |

`no_alt` is checked before every other shape: `.` and a single-base REF
are the same string length, so a no-ALT record would otherwise be
silently miscounted as a `snp` rather than excluded and reported under
its own reason.

A row whose ALT still contains a comma after `bcftools norm -m-` is a
hard error (`MalformedVcfError`): this should never happen post-split, so
silently reclassifying it as some shape would hide a real anomaly.

**Duplicate `(CHROM, POS, REF, ALT)` keys** among `snp`-classified rows
are excluded **in full** — every occurrence of a colliding key, not just
the extras — because there is no automatic way to decide which of two
identical records is "correct"; keeping one at random would be a silent,
unreviewable choice. The count and the actual colliding keys are reported
in `cohort_gs.classification_accounting.summary.txt`.

**The VCF `ID` column is never used for variant identity.** Variant
identity throughout this panel — the duplicate-key check above, the
`variant_key` used in the matrix and variant metadata, and every record-
accounting figure — is `CHROM:POS:REF:ALT` only. `ID` is carried through
unmodified into the classified VCF's output rows, but it is not retained
in any GS panel metadata file and never inspected for uniqueness.
Duplicate `ID` values across otherwise-distinct records are therefore
permitted and silently ignored, not an error.

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
| `0/0` or `0\|0` (homozygous reference) | `-1` |
| `0/1`, `1/0`, `0\|1`, or `1\|0` (heterozygous) | `0` |
| `1/1` or `1\|1` (homozygous alternate) | `+1` |
| missing or non-standard (see below) | `nan` |

Phasing (`/` vs `\|`) never affects dosage. Per the VCF specification, the
separator in a `GT` field records only whether the call is phased —
whether the alleles' parental origin is known — not which alleles or how
many are present: `0\|1` carries exactly the same allele content as
`0/1` and must resolve to the same dosage. A genotype's phasedness is
tracked purely as an informational count (`phased_genotype_count`, see
below); it is never a reason to treat a call as missing. **An earlier
revision of this classifier got this wrong** — it nulled out every
phased call regardless of its allele content — and has been corrected;
`tests/bin/test_build_gs_panel.py`'s `ClassifyGenotypeTests` now pins
phased hom-ref/het/hom-alt calls to the same dosage as their unphased
equivalents.

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

A genotype is encoded as a dosage as long as it is **diploid** with a
biallelic-index call — `0/0`, `0/1`, `1/0`, or `1/1`, in either phasing
(phasing itself never disqualifies a call; see above). Every other shape
is treated as missing in the matrix (`nan`), but counted under its own
specific reason — never folded into one undifferentiated "missing"
bucket:

| Reason | Example | Counted as |
| --- | --- | --- |
| Missing | `.`, `./.`, `0/.` (any allele position is `.`) | `missing_calls` |
| Non-diploid | `0` (haploid), `0/0/1` (triploid) | `non_diploid_calls_treated_as_missing` |
| Non-biallelic allele index | `0/2` (defensive — should not occur given the input is already biallelic-only) | `non_biallelic_index_calls_treated_as_missing` |

`cohort.gs_panel.genotype_encoding_accounting.tsv` reports the cohort-wide
total for each reason (plus `total_treated_as_missing`, the sum of all
three); the sample and variant metadata files (below) report the same
breakdown at finer granularity. `phased_genotype_count` is reported
alongside these but is **not** part of `total_treated_as_missing` — it is
a purely informational count of how many calls (standard or otherwise)
used `|`, independent of whether that call was ultimately encoded as a
dosage or as missing.

### Diploid-only constraint (genotype encoding `diploid_additive_dosage_v1`)

This constraint belongs to the *genotype encoding* schema
(`diploid_additive_dosage_v1`), which is versioned independently of the
manifest's own `schema_version` (2 since Issue #52) — the encoding did
not change when the manifest's `containers` field did.

This encoding is diploid-only by design, not merely by convention: the
classification rule above has no way to assign a meaningful dosage to a
haploid or polyploid call, so a cohort with `sample_ploidy != 2` would
have *every* genotype fall into `non_diploid_calls_treated_as_missing`,
producing a matrix that is entirely `nan` while still exiting `0` and
looking, from the pipeline's own machine-readable status, like a normal
completed run. `bin/build_gs_panel.py` and `bin/build_gs_panel_manifest.py` **each
independently** fail fast — exit `1`, before writing any output —
whenever `--sample-ploidy` (wired from `params.sample_ploidy`) is not
exactly `2`. This check is deliberately duplicated rather than trusted
to run only once: in the normal Nextflow pipeline `build_gs_panel.py`
always runs first and would already have failed on a non-diploid input,
but `build_gs_panel_manifest.py` has no way to know that when invoked on
its own (e.g. directly, outside the pipeline), and a manifest recording
`parameters.sample_ploidy` alongside `genotype_encoding.ploidy ==
"diploid_only"` for a non-2 ploidy would itself be a self-contradictory
provenance record. `cohort.gs_panel.manifest.json` records the
configured `sample_ploidy` under `parameters` so a reader never has to
guess which ploidy a given panel was built under. A generalized,
ploidy-aware encoding is tracked as future work (see Issue #1's
roadmap); until it exists, this schema simply cannot be used for a
non-diploid cohort.

This per-script fail-fast is deliberately kept even though the pipeline
also has a coarser, earlier gate: `main.nf` refuses to start *any*
process at all when `params.enable_gs_panel` is `true` (the default)
and `params.sample_ploidy` is not `2` (Issue #20; see "Scope and what
this is not" above). That pipeline-level check exists to fail fast and
cheaply, before wasting a full variant-calling run; it is not a
replacement for the checks inside `bin/build_gs_panel.py` and
`bin/build_gs_panel_manifest.py` themselves, which remain the only
defense when either script is invoked directly, outside the pipeline.

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
not a binary array format. Panel scale was originally estimated here at
thousands of markers; Issue #33's real 20-sample cohort in fact produced
9,252,873 variant rows, which is what motivated the bounded-memory build
described below. TSV parsing itself has not become the bottleneck at that
scale -- retaining the parsed result was.

The matrix header row is always written, even when there are zero
variant rows (see "Empty panel contract" below) — an empty panel must
never lose the sample list.

### Compression: what is guaranteed, and why it is not `gzip.GzipFile`

The matrix is gzip-compressed with `mtime=0`, so the same logical
content always produces byte-identical compressed output. Without that,
two runs over the same input would write different checksums into the
manifest with nothing about the data having changed, undermining the
"same input reproduces the same panel" guarantee this contract rests on.

Since Issue #44 the matrix is compressed *incrementally* rather than in
one shot, and the choice of streaming primitive is part of the contract
rather than an implementation detail. `gzip.compress(payload, mtime=0)`
delegates to zlib's own gzip wrapper, which writes `OS=3` into the
header. `gzip.GzipFile` -- the obvious streaming replacement -- writes
its own header with `OS=255`. Those two produce **identical decompressed
bytes and different files**: swapping one for the other would have
silently invalidated every previously published matrix checksum while
appearing to change nothing.

The builder therefore streams through `zlib.compressobj` with
`wbits=31`, which is that same zlib wrapper at the same compression
level. Deflate's output does not depend on how the input was chunked as
long as no intermediate flush is forced, so the result is byte-identical
to the historical one-shot form. This was verified against a live
`gzip.compress` oracle rather than a stored digest, both on the
development host (Python 3.12.3, zlib 1.3) and inside the pinned
`python:3.12` production image (Python 3.12.13, zlib 1.3.1), across
chunk sizes from one byte to the whole payload.

**The compressed representation therefore did not change, and historical
matrix checksums remain valid.** The contract is pinned by
`tests/bin/test_build_gs_panel.py` (against the live oracle) and by
`tests/modules/build_gs_panel.nf.test` (asserting the published header's
flag bits, zeroed mtime and `OS=3` inside the real container), so a
future Python or zlib that broke the equivalence fails the suite rather
than quietly changing published artifacts.

The published gzip member also carries no `FNAME` and no `FCOMMENT`: the
builder writes through a hidden staging file, and its name -- an
absolute host path -- must not travel inside a published artifact.

### Build-time memory (Issue #44)

Nothing in the builder retains a quantity that scales with the variant
count. It makes one pass over the GS-eligible PASS VCF and, for each
data row, classifies that row's genotypes once and immediately writes
the matrix row into the streaming compressor, writes the variant
metadata row to an open handle, and folds the row into per-sample and
cohort-wide counters. The row is then dropped.

What is retained is the sample names, two per-sample counter lists, a
fixed-size set of cohort counters, the single row in hand, and the I/O
and compression buffers -- all `O(sample_count)` or constant. Sample
metadata, the accounting TSV and the summary are derived from the
accumulated counters after the scan, and are themselves
`O(sample_count)` or fixed size.

Each genotype is classified exactly once. The previous implementation
classified every cell five times over: once per output, plus a second
accounting pass to build the summary.

`bin/build_gs_panel.py` still defines `parse_gs_pass_vcf` and the
`build_*_rows` family, which do materialize whole documents. They are
**not** on the production path; they state each output's content
declaratively and serve as the oracle the streaming implementation is
tested against. A test replaces each of them with a raising stub and
runs the CLI, so a change that routed the production path back through
one of them fails immediately rather than at real-cohort scale.

Measured figures are in `docs/gs_panel_streaming_benchmark.md`.

## Output publication and failure semantics

The five panel artifacts -- matrix, sample metadata, variant metadata,
genotype accounting, and the accounting summary -- are built into
staging files beside their final paths and moved into place only after
the whole VCF has been read and all five documents have been produced.

What that guarantees:

- A malformed row **anywhere** in the input, including its last line,
  leaves no final output at all. A streaming builder has necessarily
  already written many good rows by then; none of them are published.
- A run that fails does not disturb the artifacts of a previous
  successful run.
- Each individual `os.replace` is atomic, and a failure partway through
  the publish sequence restores every already-replaced file from the
  copy moved aside moments earlier.
- Staging files are swept on every exit path, including the `return 1`
  paths that no exception handler would see.

What it does **not** guarantee: publication of the five files is not a
filesystem transaction. `os.replace` is atomic per file and there is no
primitive that commits five of them together. A process kill or power
loss in the middle of the replace sequence, or a failure during the
rollback itself, can leave a mix of new and old files on disk. The
rollback narrows that window to the sequence itself rather than to the
whole build; it does not eliminate it. Downstream readers should
continue to rely on `reconcile_gs_panel_accounting.py`, which re-reads
the artifacts and cross-checks them against the source VCF, rather than
assuming the set is internally consistent because it exists.

### Input validation

A data row must carry exactly `9 + sample_count` tab-separated fields,
checked against the sample list the `#CHROM` header declared. Both a
shortage and an excess are hard errors, at any point in the file. Before
Issue #44 the reader only required "at least 10" fields, so a row
carrying fewer genotypes than the header declared was accepted and
produced a matrix row narrower than its own header -- a column
misalignment the builder should never emit.

Also rejected, each as a diagnosable error naming the line: a FORMAT
column with no `GT`, a sample field with fewer subfields than FORMAT's
`GT` index requires (previously a bare `IndexError` traceback), a data
row before `#CHROM`, a second `#CHROM` header, a `#CHROM` header with no
sample columns, and a file with no `#CHROM` header at all.

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
| `non_standard_genotype_count` | Subset of `missing_genotype_count` caused by non-diploid/non-biallelic-index calls specifically (excludes plain missing; a phased call is not "non-standard" by itself — see "Dosage encoding" above) |

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
  -> normalized_records                 (post bcftools norm -m- splitting)
  -> classified_biallelic_snp_records   (post shape reclassification and duplicate-key exclusion)
  -> gs_pass_records                    (post GS-specific hard filtering)
  -> matrix_variant_records / variant_metadata_records   (read directly; must equal gs_pass_records and each other)
```

with sample counts reconciled the same way: `gs_pass_sample_count`,
`matrix_sample_count`, and `sample_metadata_records` must all agree.
`raw_all_records` and `normalized_records` are independently re-counted
directly from their VCFs (not trusted from any other script's output),
matching `bin/reconcile_variant_type_counts.py`'s established
"independently re-verify" convention. `classified_biallelic_snp_records`
is read from `classify_normalized_variants.py`'s own accounting rather
than re-derived, since re-implementing its shape-classification and
duplicate-key logic a second time would duplicate that logic rather than
verify it. Critically, `matrix_variant_records` and `matrix_sample_count`
are read directly from the genotype matrix file itself (its own header
and row count), not derived from the metadata files that were supposed
to describe it.

**Counting alone is not enough.** Four artifacts can each report the
same variant count while still disagreeing about *which* variant is in
which row, or listing the same rows in a different order: a shuffled
matrix sample column, a `variant_key` swapped for a different value, or
`variant_metadata.tsv`'s rows written in a different order than the
matrix. None of these change any count, so this tool compares the
actual, ordered sequence of values, not just their lengths --
`gs_pass_vcf.sample_ids == matrix.sample_ids == sample_metadata`'s own
`sample_id` column, and `gs_pass_vcf.variant_keys ==
matrix.variant_keys == variant_metadata`'s own `variant_key` column,
both compared as ordered tuples. A tuple-equality check on identity
subsumes the plain count check it replaced (two sequences cannot be
equal without also being the same length). The genotype matrix's own
row shape is checked too: every data row must have exactly
`1 + sample_count` columns, which catches a row with a silently dropped
or duplicated dosage column that a row-count-only check (the matrix
summarizer's previous behavior) would never see.

**Every one of these cross-file checks is a hard error
(`InconsistentGsPanelError`), not a warning: on any disagreement, the
tool exits `1` and writes no output at all.** This is a deliberate
departure from `bin/reconcile_variant_type_counts.py`'s convention of
reporting a negative `records_not_selected` as a `WARNING` rather than
failing the run. That warning is appropriate there because a negative
value in the *primary* lineage has a real, understood, non-buggy cause
(GATK's own multiallelic-record handling, documented in this repository's
Issue #15 work) — it is surprising but scientifically real. Here, there
is no equivalent legitimate scenario: hard filtering can only remove
records, never add them, so a negative
`gs_hard_filter_excluded_records` can only indicate a bug; and there is
no valid reason the matrix, the PASS VCF it was built from, and the
metadata files describing it should ever disagree on how many variants
or samples they contain. Once one of the GS panel's own artifacts
disagrees with its own source, none of the remaining numbers can be
trusted either, so reporting a partial, possibly-wrong accounting
alongside a "warning" would be worse than refusing to produce one.

This replaces two earlier, narrower revisions of this tool: one derived
`final_matrix_variant_records` from the variant metadata file's row
count without ever opening the matrix itself, so a bug that corrupted
only the matrix (a dropped row, a misaligned column, a truncated write)
would have gone completely undetected as long as the metadata files
still looked correct; a later revision opened the matrix but still only
compared counts, so a shuffled sample column or a reordered metadata
file would still have passed silently. The current tool also checks
that every data row in the raw/all, normalized, and GS-eligible PASS
VCFs has exactly as many sample fields as its own `#CHROM` header
declares, which catches a truncated or otherwise malformed row that a
plain field-count-only check would miss.

## Reproducibility manifest (concern 6)

`cohort.gs_panel.manifest.json` (`bin/build_gs_panel_manifest.py`) is
schema-versioned (`schema_version: 2`) and mirrors the *shape* of the
sibling repository's own `run_manifest.py` — a sortable `run_id`
(`<UTC-timestamp>-<uuid4-hex8>`), deterministic canonical JSON, a
self-referential `manifest_hash`, and filename-only checksums (never an
absolute path) — without depending on that repository's code: a
Nextflow bioinformatics pipeline has no reason to import a model-training
repository's Python package, so the pattern is reimplemented in
stdlib-only Python instead.

```json
{
  "schema_version": 2,
  "run_id": "20260814T074225Z-c22d6e72",
  "generated_at": "2026-08-14T07:42:25Z",
  "cohort_id": "cohort",
  "pipeline_version": "0.2.0-dev",
  "git_commit": null,
  "containers": {
    "gs_normalize_variants": "...",
    "classify_normalized_variants": "...",
    "gs_index_classified_variants": "...",
    "gatk_variantfiltration_gs": "...",
    "gatk_selectpassvariants_gs": "...",
    "build_gs_panel": "...",
    "reconcile_gs_panel_accounting": "...",
    "build_gs_panel_manifest": "..."
  },
  "parameters": { "sample_ploidy": 2, "snp_filter_qd_min": 2.0, "...": "..." },
  "genotype_encoding": {
    "schema": "diploid_additive_dosage_v1",
    "dosage_by_genotype": { "0/0": -1, "0/1_or_1/0": 0, "1/1": 1 },
    "phasing": "ignored for dosage; 0|1 encodes identically to 0/1",
    "missing_token": "nan",
    "matrix_orientation": "variant_rows_by_sample_columns",
    "ploidy": "diploid_only"
  },
  "panel_status": "empty",
  "checksums": { "cohort.gs_panel.genotype_matrix.tsv.gz": "sha256:...", "...": "..." },
  "manifest_hash": "sha256:..."
}
```

`genotype_encoding` is a fixed, static description of this schema (see
"Genotype matrix contract" above) recorded verbatim in every manifest, so
a reader never has to cross-reference `bin/build_gs_panel.py`'s source or
this document to know how to interpret the matrix. `checksums` covers not
only the panel's own deliverables (matrix, sample/variant metadata,
genotype-encoding accounting, record accounting) but also the inputs used
to build them — the raw/all cohort VCF and the reference FASTA/FAI used
for normalization — so a reader can reconstruct *which* inputs and
reference a given panel came from, not only verify the panel's own
outputs against themselves.

**`containers` (schema v2, Issue #52): one entry per GS-lineage process,
never merged by shared tool.** Each key is the exact GS process name
(lowercased) of every process the GS lineage actually runs, including
`build_gs_panel_manifest` — the process that writes this very file —
`gs_normalize_variants` and `gs_index_classified_variants` both use
bcftools by default; `build_gs_panel`, `classify_normalized_variants`,
`reconcile_gs_panel_accounting`, and `build_gs_panel_manifest` all use
Python by default; `gatk_variantfiltration_gs` and
`gatk_selectpassvariants_gs` are the GS lineage's own aliases of the
shared `GATK_VARIANTFILTRATION`/`GATK_SELECTPASSVARIANTS` modules — but
each value is that specific process's own *effective* container, not a
value shared across an entire tool category. Two processes that happen to
use the same tool by default can still diverge (a
`withName`/alias/fully-qualified-selector/profile override on only one of
them) and schema v2 is designed so that can never be silently collapsed
into a single ambiguous value.

The GS lineage runs exactly these eight processes, and all eight are
recorded: this field is the whole GS lineage, with no "except the one
that generated the record" exclusion. Processes outside the GS lineage
that happen to share these same images by default
(`SUMMARIZE_VARIANT_QC`, `BCFTOOLS_STATS`, the primary lineage's own
GATK stages) are deliberately *not* in scope here — a GS panel manifest
describes how the GS panel was produced; run-level provenance covering
every process in the pipeline is tracked separately as Issue #42.

Each value is Nextflow's own `task.container`, captured from inside that
exact task via a `container_id` output added to every upstream GS-lineage
process (see `modules/local/gs_normalize_variants.nf` and its sibling GS
modules), and — for `build_gs_panel_manifest` itself — read directly as
`task.container` in `BUILD_GS_PANEL_MANIFEST`'s own script, since
Nextflow resolves a task's container before rendering that task's script
(verified on Nextflow 26.04.6 against both a default and a `withName`
override). Either way it is the container **actually used**, already
resolved after any `withName`/alias/fully-qualified-selector/profile
override on top of that process's `container` directive default. The default itself now lives in
one place, `conf/containers.config` (`params.containers.bcftools/gatk/
python`), referenced by each module's `container` directive instead of
repeating the literal — but the manifest never trusts that default value
directly, precisely because a default can be overridden per-task without
this file changing at all. `workflows/adzuki_snp_pipeline.nf` wires each
upstream process's `container_id` output straight into
`BUILD_GS_PANEL_MANIFEST`; no literal container digest is defined in the
workflow file, nor in any GS module file, any more.

`tests/pipeline/adzuki_snp_pipeline.nf.test` holds the end-to-end proof
of that claim rather than restating it: one test changes
`params.containers` from a test config alone and asserts the changed
value appears both in the tasks' real `docker run` invocation
(`.command.run`) and in this manifest, while a process-specific
`withName` override on one process still wins over the changed default
and does not leak onto its siblings. The expected default values in those
tests are parsed out of `conf/containers.config` at test time rather than
copied into the test, so bumping a pinned image stays a one-file change.

**Schema v1 → v2: why this was a version bump, not an additive change.**
Schema v1's `containers.bcftools/gatk/python` was keyed by *tool
category*, one shared value for every process using that tool. That
meaning becomes ambiguous the moment two processes sharing a tool are
overridden differently (e.g. only `GS_NORMALIZE_VARIANTS`'s bcftools
container, not `GS_INDEX_CLASSIFIED_VARIANTS`'s) — there would be no
correct single value to put in `containers.bcftools`, and keeping the
field while its meaning became unclear was rejected as worse than
versioning: the field's own keys had to change from tool names to process
names, which is not a shape a v1 reader could parse as an extension of
the same field. Older, already-published `schema_version: 1` manifests
(including historical run/scale-validation manifests referenced elsewhere
in this repository's docs) are untouched by this change and remain valid
schema v1 documents; only newly generated manifests use schema v2.

Software versions are therefore recorded as each GS process's own
effective container identity (see `containers` above) rather than by
shelling out to e.g. `bcftools --version` at run time.

`git_commit` is Nextflow's own `workflow.commitId`, resolved by the
workflow and passed in as a plain argument rather than computed inside a
task container (a task's work directory holds only staged input files,
never the pipeline's `.git` checkout, so a container-side `git
rev-parse` would find no repository at all). `workflow.commitId` is
documented as populated only when Nextflow itself pulls a git-hosted
pipeline (`nextflow run owner/repo`); it is legitimately `null` for this
repository's own documented `nextflow run .` local-directory invocation,
and is reported as such honestly rather than worked around.

**Checksum reproducibility, precisely scoped.** The manifest document
itself (`cohort.gs_panel.manifest.json`) is **never** byte-for-byte
reproducible across independent runs, by design: it embeds a fresh
`run_id` and `generated_at` timestamp on every invocation, and — since it
now also checksums the raw/all and GS-eligible PASS VCFs for provenance
(see above) — it transitively embeds two checksums that are themselves
run-dependent. "Reproducible" therefore never means "the manifest file
matches"; it means a specific, named subset of the *checksums recorded
inside it* match. That subset is exactly the seven files this pipeline
writes purely from its own inputs and counts, with no embedded
timestamp:

- `cohort.gs_panel.genotype_matrix.tsv.gz`
- `cohort.gs_panel.sample_metadata.tsv`
- `cohort.gs_panel.variant_metadata.tsv`
- `cohort.gs_panel.genotype_encoding_accounting.tsv` and its `.summary.txt`
- `cohort.gs_panel.record_accounting.tsv` and its `.summary.txt`

These are byte-for-byte reproducible across independent runs over
identical input (`bin/build_gs_panel.py`'s matrix writer explicitly sets
`mtime=0` on its gzip output for exactly this reason), verified directly
by running the generated- and prebuilt-reference-index paths and diffing
each of the seven files above.

Every other checksum the manifest records —
`cohort_gs.snp.pass.vcf.gz` (the GS-eligible PASS VCF) and
`cohort.raw.vcf.gz` (raw/all) — is **not** reproducible run-to-run: GATK
embeds the actual wall-clock run time in its `##GATKCommandLine` VCF
header on every invocation that touches a VCF, a property of every VCF
this pipeline produces (not something introduced by, or fixable within,
the GS panel work), so these checksums will legitimately differ between
two runs even when every underlying record is identical. The reference
FASTA/FAI checksums are reproducible in practice (the reference file
itself does not change between runs against the same input), but that is
a property of the input, not something this pipeline guarantees or
verifies.

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
