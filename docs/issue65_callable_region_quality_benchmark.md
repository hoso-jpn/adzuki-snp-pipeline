# Variant-quality evidence by region and the delivery-support matrix (Issue #65)

This document records **which evidence exists, of which strength, for which
region, depth and variant type** on the pipeline's reference
GCF_016808095.1. It is the human-readable side of the machine-readable
records in [`evidence/issue65/`](evidence/issue65/).

**No accuracy is claimed anywhere in this document.** The audit found no
independent truth set for *Vigna angularis*, so no precision or recall exists
for real data. Every real-data number below is a concordance, a stability or a
descriptive count, and is named that way. No caller, threshold, sample rule,
MAF rule or imputation was adopted or changed: **production defaults are
unchanged**, and the only repository changes are the benchmark code, its
tests, and these records.

## 1. Evidence-class model

The contract lives in [`benchmarks/issue65/evidence_model.py`](../benchmarks/issue65/evidence_model.py).
Each class names the metrics it may carry, and an evaluation record that
reports a metric outside its class is rejected. The accuracy vocabulary
(`true_positive`, `precision`, `recall`, `f1`) belongs only to independent
truth and to the synthetic fixture that tests the harness.

| evidence class | metrics it may carry | strongest delivery status | what it can claim |
|---|---|---|---|
| `independent_truth` | TP/FP/FN, precision, recall, F1, GT concordance | `supported` | accuracy against that truth, in its region |
| `technical_replicate_concordance` | shared / only-in-one / no-call, `site_agreement_rate`, GT concordance | `supported_with_caveat` | reproducibility; not accuracy |
| `cross_platform_concordance` | concordance and descriptive metrics | `supported_with_caveat` | agreement with another platform on the same sample; not accuracy |
| `caller_concordance` | concordance metrics | `supported_with_caveat` | agreement between callers on the same reads; not accuracy |
| `downsampling_stability` | shared / only-in-full / only-in-downsampled, `call_retention`, `new_call_fraction`, GT concordance | `supported_with_caveat` | how calls change as depth falls; not accuracy |
| `descriptive_stratification` | counts, rates, missing / masked / heterozygous fractions | `supported_with_caveat` | counts by stratum; no quality claim |
| `synthetic_truth_fixture` | truth metrics | `validated_test_harness` | the harness counts correctly |

Statuses are `supported`, `supported_with_caveat`, `not_evaluated`,
`unsupported` and `validated_test_harness`. `delivery_row` refuses a status
stronger than its class allows, a `not_evaluated` cell that cites evidence,
and any other status that cites none.

## 2. Evaluation contract

Every record carries `evaluation_schema_version`, `evaluation_id`,
`reference`, `region_set`, `region_definition_hash`, `evidence_class`,
`query_dataset`, `comparator_dataset`, `comparison_unit`,
`eligible_denominator`, `excluded_denominator`, `exclusion_reasons` (which must
sum to the excluded denominator), `not_evaluated_reason`, `metrics` and
`limitations`.

- **`not_evaluated` is a record, not a missing row.** It has a non-empty
  reason and no metrics.
- **Comparison unit.** The unit is one ALT allele per sample, after
  multi-allelic splitting, trimming and left-alignment. That is the rule
  `bcftools norm -f ref -m -any` implements, and the unit tests pin the Python
  implementation against output from the pinned bcftools 1.24 container.
- **Exclusions are counted, never folded into a metric.** The reasons are
  `outside_region`, `region_boundary` (the span only partly overlaps the
  region), and, per side, `coordinate_mismatch`, `ref_mismatch` and
  `symbolic_allele`.
- **True negatives are not counted.** A VCF does not enumerate confident
  reference positions. A rate whose denominator is zero is `null`.
- **Coordinates.** BED is 0-based half-open and VCF is 1-based. The only
  conversion happens in `RegionSet.contains_span`. A region set is bound to the
  reference's ordered contigs and lengths, and an interval outside them is an
  error, not clipped.
- **Tags vs partitions.**
  - A *partition* group (GC bins, depth bins, mappability, callable,
    core/difficult) is checked to be disjoint and to cover its universe
    exactly, so its members' counts add up.
  - A *tag* (WindowMasker repeat, homopolymer, N) may overlap anything, and
    tag counts are not additive.
  - A record that straddles a stratum edge is counted in neither side and
    reported separately.
- **Splits.** Every comparison stratum is also split into `<stratum>:snp`,
  `:indel`, `:het` and `:hom_alt`. The dosage comes from the comparator's call,
  or the query's where the comparator is not positive.
- **Identity.** `region_definition_hash` is a SHA256 over the reference
  SHA256, contigs, merged intervals and provenance. The same membership always
  gives the same hash, and any membership change alters it (unit-tested).

## 3. Public asset audit (Phase 1)

The full inventory is [`benchmarks/issue65/public_asset_inventory.json`](../benchmarks/issue65/public_asset_inventory.json).
Each entry records source, accession/URL, organism, identity, platform,
assembly, coordinate system, version, checksum, licence, evidence class,
usable / unusable / not_evaluated, and the reason. Identity is established only
by BioSample accession, never by cultivar name.

| asset | status | reason |
|---|---|---|
| GCF_016808095.1 RefSeq assembly (SAMN14776547, PacBio) | usable (reference) | the pipeline's reference |
| RefSeq `genomic_gaps.txt` (32 gaps) | usable | every gap lies inside an N run of the FASTA (checked) |
| RefSeq lower-case bases (WindowMasker) | usable | NCBI README documents lower case as WindowMasker repeats |
| RefSeq RepeatMasker output | unusable | not published for this assembly |
| SRR11787767, Illumina NovaSeq WGS of SAMN14776547 | usable, cross-platform | same BioSample as the reference; not truth |
| SRR11787766, PacBio Sequel II of SAMN14776547 | not_evaluated | not independent of the reference; no pinned long-read caller |
| SRR11787765, Hi-C of SAMN14776547 | unusable | proximity-ligation library |
| PRJNA1138464 (327 WGS + 357 RAD, 684 BioSamples) | unusable for replicate or cross-method | one run per BioSample; no BioSample has both WGS and RAD |
| SAMN03488483 (68 GAII WGS runs, 58 libraries) | unusable | a "Resequencing" series under one BioSample; its same-library run pairs are mate-pair libraries or 4x-vs-20x runs |
| PRJNA1092869 and PRJNA1117856 long reads | unusable | no BioSample in common with the reference or the cohort |
| *V. angularis* small-variant truth | not found | ENA analyses for taxon 3914 are all REFERENCE_ALIGNMENT; EVA lists no *V. angularis* assembly |
| Human GIAB truth, stratifications, human-trained models | unusable by design | another organism |

**Truth candidates: none.**

- **Technical replicates:** none that are identity-confirmed and usable.
- **Cross-platform:** only the reference BioSample itself, which is
  self-consistency evidence and not truth.

## 4. Synthetic truth fixture (Phase 3)

[`make_synthetic_fixture.py`](../benchmarks/issue65/make_synthetic_fixture.py)
writes a 400 bp two-contig reference, a truth VCF, a query VCF and region BEDs.
Each record exercises one rule. The expected counts in
[`tests/bin/test_issue65_quality_evidence.py`](../tests/bin/test_issue65_quality_evidence.py)
were counted by hand from the case comments, not by running the engine.

| case | rule | expected |
|---|---|---|
| chrT:10, chrU:10 | both 0/1 | TP, GT concordant |
| chrT:20 | truth 1/1, query 0/1 | TP, GT discordant |
| chrT:30 / chrT:70 | no truth record / truth 0/0 | FP |
| chrT:40 / chrT:80 | no query record / query 0/0 | FN |
| chrT:50, chrT:90 | query `./.`, `0/.` | query no-call at truth variant |
| chrT:60 | truth `./.` | truth no-call (not in precision) |
| `CA>C`@99 vs `AA>A`@103 | same deletion in a homopolymer | TP after normalization |
| `G>GTT`@121 | insertion | TP, indel |
| truth `A>C,T 1/2` vs split query | multi-allelic | 2 TP |
| query `C>A,G 1/1` vs truth `C>A` | multi-allelic query | 1 TP |
| REF `T` where reference is `G` | REF mismatch | excluded |
| `<DEL>` | symbolic | excluded |
| `chrZ` | contig not in reference | excluded (coordinate mismatch) |
| 4 bp deletion across region end | boundary | excluded |
| chrT:250 | outside region | excluded |

The harness gives TP 8, FP 2, FN 2, query no-call 2 and truth no-call 1.
Precision is 8/10 and recall is 8/12, where recall's denominator includes
query no-calls. There are 15 eligible units and 5 exclusions (one per reason).
Partition strata add up, tags do not, and the output is deterministic. With the
same records, `caller_concordance` reports `site_agreement_rate` 8/15 and
`downsampling_stability` reports `call_retention` 8/12, with no accuracy words.

## 5. Region assets on GCF_016808095.1 (Phase 4)

All assets were computed from the exact FASTA (SHA256 `e9838db1…`), with fixed
parameters, by [`build_region_assets.py`](../benchmarks/issue65/build_region_assets.py)
and [`mappability.py`](../benchmarks/issue65/mappability.py). No human or other
mask was used. Hashes are in `evidence/issue65/region_assets.json` and
`mappability.json`.

| asset | kind | parameters | bases | fraction |
|---|---|---|---|---|
| `n_regions` | tag | runs of N | 3,200 | 0.0007% |
| `repeat_windowmasker` | tag | lower-case runs | 162,429,013 | 36.2% |
| `homopolymer_ge10` | tag | single-base run >= 10 | 1,089,560 | 0.24% |
| `gc_00_25` … `gc_45_100`, `gc_undefined` | partition | 1 kb windows, edges 25/30/35/40/45%, undefined if > half N | 58.4 / 68.8 / 124.2 / 122.8 / 49.2 / 24.9 / 0 Mb | |
| `high_mappability` | partition | 150 bp error-free tiles every 50 bp, BWA-MEM2 (pipeline image), primary alignment at origin with MAPQ >= 30 | 373,737,150 | 83.4% |
| `low_mappability` | partition | same, failing that rule | 74,613,350 | 16.6% |
| `mappability_undefined` | partition | no tile starts there (N, contig end) | 12,142 | 0.003% |

Resources: region assets 16.9 s at 289 MiB RSS. Mappability took 39.8 s wall
for 8,967,010 tiles aligned with 30 threads. Two mappability runs produced
byte-identical BEDs.

## 6. Callable and coverage (Phase 5)

This reuses Issue #45's existing 51 markdup BAMs on the Issue #45 scaling
window `NC_068975.1:1-20000000`. No FASTQ was reprocessed.

- `samtools depth -a -Q 20 -q 10` (samtools' default flag filter drops
  unmapped, secondary, QC-fail and duplicate reads).
- **Sample-level:** a sample is callable at a position when `5 <= depth <= 2.5 x`
  that sample's median depth over the window.
- **Cohort-level:** callable where at least 80% of the 51 samples are callable.
- **Reference-level:** N bases have depth 0 and are never callable.
- Cohort median-depth bins: 0-3, 3-6, 6-10, 10-15, 15-25, 25+.

These are benchmark stratification parameters, **not pipeline defaults**, and
"callable" describes alignment depth and quality. **It does not mean the
genotypes there are accurate.**

| result | value |
|---|---|
| cohort callable / non-callable | 15,142,446 / 4,857,554 bp (75.7% / 24.3%) |
| per-sample median depth over the window | 3 to 22 (median 11) |
| per-sample callable fraction | 0.20 to 0.95 (median 0.90); 7 samples below 0.5 |
| depth pass / labelling pass | 43.8 s / 297 s wall |

## 7. #64 genotype mask by region (Phase 6)

These are descriptive counts of the Issue #64 51-sample GS SNP PASS panel
(SHA256 `f1ab1c3a…`, unchanged) against the same panel after the illustrative
Issue #64 mask (`min_dp=10`, `min_gq=20`, SHA256 `c2fa4567…`). The two sides
are **site filter only** and **site filter + genotype mask** over the identical
records and samples (the tool fails if the record lists differ). No sample QC
is mixed in.

The mask policy is illustrative, not calibrated, and not a default.
Differences between strata mix sequence context, depth and cohort composition,
and cannot be separated here.

Genome-wide (9,746,661 records, 497,079,711 genotype cells):

| stratum | records/Mb | missing (site filter only) | masked by DP/GQ | het share of non-ref calls |
|---|---|---|---|---|
| all | 21,738 | 6.4% | 42.4% | 44.1% |
| high mappability | 24,120 | 6.1% | 43.0% | 40.9% |
| low mappability | 9,809 | 10.0% | 35.8% | 78.7% |
| WindowMasker repeat | 20,267 | 7.5% | 43.8% | 47.0% |
| homopolymer >= 10 | 25,768 | 7.7% | 49.8% | 20.7% |
| GC < 25% | 18,688 | 4.9% | 56.2% | 10.3% |
| GC >= 45% | 28,513 | 10.6% | 37.1% | 64.3% |

Window `NC_068975.1:1-20000000` (580,258 records):

| stratum | bases | records/Mb | missing | masked | het share |
|---|---|---|---|---|---|
| core | 8,262,745 | 34,423 | 5.8% | 39.5% | 54.1% |
| difficult | 11,737,255 | 25,205 | 8.6% | 37.6% | 72.8% |
| cohort callable | 15,142,446 | 32,123 | 6.6% | 38.8% | 59.9% |
| cohort non-callable | 4,857,554 | 19,318 | 10.4% | 37.0% | 79.2% |
| median depth 0-3 | 880,930 | 2,586 | 55.0% | 38.2% | 13.9% |
| median depth 6-10 | 3,321,020 | 19,681 | 6.2% | 55.9% | 16.7% |
| median depth 10-15 | 13,847,094 | 29,619 | 6.7% | 38.8% | 47.0% |
| median depth 15-25 | 1,081,320 | 70,906 | 8.9% | 26.2% | 95.3% |
| median depth 25+ | 244,862 | 86,375 | 7.8% | 16.5% | 97.2% |

What this shows is descriptive only:

- The DP/GQ mask removes the most cells where depth is low (3-10x, 56-74%) and GC is
  low.
- Strata with 1.5-2.5x the typical cohort depth carry 2.4-3.0x the record
  density, and 95-97% of their non-reference calls are heterozygous. Low
  mappability shows the same pattern. In a self-pollinating crop this is the
  signature expected of collapsed paralogs or mis-mapping, but that is an
  interpretation and was not verified.
- The mask removes *fewer* calls there, because depth is high. A DP/GQ mask
  therefore does not address this class of call.

Resources: genome-wide pass 229 s and window pass 22 s, two bcftools query
streams into one Python process.

## 8. Downsampling stability (Phase 7)

- **Input:** the existing markdup BAM of SRR29908806, the most deeply sequenced
  of the 51 (window mean depth 23.1x at MAPQ>=20/BQ>=10). No FASTQ was
  reprocessed.
- **Subsampling:** window reads with `samtools view --subsample F
  --subsample-seed S`. This hashes read names, so mates stay together.
- **Design:** fractions 0.5 and 0.25, each with seeds 65 and 66, plus the
  full-depth BAM. Every BAM went through the production HaplotypeCaller
  arguments (`--intervals` the window), then single-sample GenotypeGVCFs.
- **Comparator:** the full-depth calls. They are not truth, so retention
  measures agreement with them.
- The two seeds at one fraction draw overlapping reads. Their difference shows
  sampling variation, not replication.

| subsample | reads | mean depth | SNP retention | SNP new-call fraction | SNP GT concordance | indel retention |
|---|---|---|---|---|---|---|
| 0.5, seed 65 | 2,546,657 | 11.54 | 0.710 | 0.145 | 0.975 | 0.716 |
| 0.5, seed 66 | 2,548,214 | 11.54 | 0.711 | 0.146 | 0.975 | 0.724 |
| 0.25, seed 65 | 1,273,146 | 5.76 | 0.488 | 0.131 | 0.961 | 0.499 |
| 0.25, seed 66 | 1,273,115 | 5.77 | 0.485 | 0.133 | 0.961 | 0.505 |

(the full-depth set has 5,095,926 reads, mean depth 23.09)

The by-stratum, het/hom-alt and SNP/indel figures are in
`delivery_support_matrix.tsv` and `evaluations.json`.

This sample's full-depth calls are 82% heterozygous, and most of those have an
alternate-allele fraction of 10-20% (§10). Such calls are depth-sensitive, which
is why overall retention at half depth is only about 0.71. Retention is not a
statement about accuracy.

Resources: the whole run took 6 min 37 s wall, with the five HaplotypeCaller
jobs in parallel at 4 threads.

| step | wall | memory peak |
|---|---|---|
| HaplotypeCaller, full / 0.5 / 0.25 | 308 / 196 / 117 s | 1.87 / 2.13 / 2.03 GiB |
| GenotypeGVCFs | 10-16 s | 0.96 GiB |
| window BAM extraction | 10-22 s | 94-309 MiB |

Storage for all BAMs and VCFs was 1.06 GB. Memory is the container cgroup
`memory.peak`, polled each second.

## 9. Caller concordance

bcftools 1.24 (`mpileup -q 20 -Q 10 | call -m -v`, pinned container) was run
against the full-depth HaplotypeCaller calls on the same BAM and window. Both
callers read the same alignments, so a shared mapping error agrees with itself:
this is agreement, not accuracy.

| split | shared | only bcftools | only HaplotypeCaller | site agreement | GT concordance |
|---|---|---|---|---|---|
| all | 45,023 | 13,896 | 22,032 | 0.556 | 0.978 |
| SNP | 42,338 | 13,527 | 17,840 | 0.574 | 0.985 |

bcftools took 42.9 s at 92 MiB.

## 10. Cross-platform self-consistency

CROSSPLATFORM_SECTION

## 11. Delivery-support matrix (Phase 8)

Files:

- [`evidence/issue65/delivery_support_matrix.tsv`](evidence/issue65/delivery_support_matrix.tsv)
- [`delivery_support_matrix.json`](evidence/issue65/delivery_support_matrix.json)
  (rows and every per-class cell)
- [`evaluations.json`](evidence/issue65/evaluations.json) (all evaluation and
  not-evaluated records)

The rule is in [`assemble_evidence.py`](../benchmarks/issue65/assemble_evidence.py)
and introduces no numeric threshold:

- **Truth and replicate cells** are `not_evaluated` everywhere, because no
  usable asset exists.
- **A comparison or descriptive cell** is `supported_with_caveat` when its
  evaluation has at least one unit in the scope, and `not_evaluated` otherwise.
- **A scope** is `unsupported` when it lies outside the callable definition
  (`cohort_non_callable`), `supported_with_caveat` when any of its cells is,
  and `not_evaluated` otherwise.

MATRIX_SECTION

## 12. Not evaluated

| item | reason |
|---|---|
| independent truth, every scope | no *V. angularis* truth set found; the reference BioSample's long reads are not independent of the reference |
| technical replicate concordance | no identity-confirmed replicate (PRJNA1138464 one run per BioSample; SAMN03488483 unsuitable) |
| RAD vs WGS cross-method in the cohort | no BioSample with both; no name matching |
| long-read calls of SAMN14776547 | not independent of the reference; no pinned long-read caller |
| GS panel indel stratum | the GS panel is SNP-only |
| callable / depth outside the window | assets built from existing BAMs for the scaling window only |
| downsampling on other samples / whole genome | one BAM, one window |
| repeat-family strata | no RepeatMasker family output for this assembly |

## 13. Reproducibility

The run scripts used on seedcore-01 live under [`benchmarks/issue65/run/`](../benchmarks/issue65/run/).
Each script pins its container images by digest and records input and output
SHA256. Re-running from the recorded inputs is:

1. `run_mappability.sh`, then `build_region_assets.py`.
2. `run_callable.sh`.
3. `assemble_evidence.py derive-strata`.
4. `run_stratify.sh`, `run_downsampling.sh` and `run_crossplatform.sh`.
5. `assemble_evidence.py assemble --git-sha <sha>`.

Determinism checks:

- The mappability rerun gave byte-identical BEDs.
- `assemble_evidence.py assemble` was run twice on the final inputs and gave
  byte-identical `evaluations.json`, `delivery_support_matrix.json/.tsv` (see
  `evidence_manifest.json`).
- Unit tests assert identical records, denominators, membership and region
  hashes on re-runs.

Large intermediates (BAMs, VCFs, BEDs, depth stream) stay on seedcore-01 and
are identified by SHA256. Issue #45 and #64 evidence was read, never modified.

## 14. Limitations

- There is no truth, so no statement about accuracy is possible for any
  region.
- The real-data comparisons use one sample, one 20 Mb window of one
  chromosome, and single-sample raw calls without the cohort hard filter. They
  do not describe the delivered GS panel's accuracy.
- The callable and depth assets are limited to that window.
- Mappability models error-free 150 bp single-end reads. Real reads differ in
  errors, length and pairing.
- The cross-platform reads are the leading records of the run, not a random
  sample.
- The high heterozygous share in high-depth and low-mappability strata is
  consistent with paralog collapse, but that was not verified.
