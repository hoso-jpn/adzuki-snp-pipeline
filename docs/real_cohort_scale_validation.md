# Real cohort scale validation: 10 → 20 → 30 samples (Issue #33)

This document records Issue #33: a Gated, staged expansion of the real public-data cohort
from Issue #26's 5 samples toward the 20-30 sample scale Issue #26/#30 both left as an
unattempted, conditional next step. Each stage (10, then 20, then 30 samples) only proceeds
after the previous stage's Gate criteria are explicitly confirmed against real measurements;
a 327-sample full cohort remains out of scope regardless of this Issue's outcome.

This document is written incrementally as each stage completes -- sections below are filled
in as real evidence becomes available, not written in advance of the run they describe.

## Purpose

Determine, from real execution evidence (not extrapolation from Issue #26's 5-sample
baseline alone), whether 20-30 real samples is an operationally viable next step for
commissioned analysis on `seedcore-01`, and update Issue #30's **CONDITIONAL GO** decision
accordingly.

## Prerequisites confirmed before any real-data work (Stage 0)

- PR #29 (Issue #26) and PR #31 (Issue #30) merged into `main`; PR #32 (a follow-up
  comment-only fix) also merged. Issues #26 and #30 both closed.
- Working tree re-synced to `main` at `a3fe83e` on both the local checkout and `seedcore-01`.
- Baseline regression re-measured fresh on this exact commit (not reused from Issue #30's
  own baseline numbers):

  ```
  $ python3 -m unittest discover -s tests/bin
  Ran 252 tests in 0.220s — OK

  $ nextflow lint .
  ✅ 35 files had no errors

  $ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test tests/modules/gatk_genomicsdbimport.nf.test tests/modules/gatk_gathervcfs.nf.test --ci
  SUCCESS: Executed 26 tests in 148.846s

  $ git diff --check
  (exit 0)
  ```

## Environment (measured fresh, not reused from Issue #26/#30's own measurements)

Measured 2026-08-25, on host `seedcore-01`:

- Hardware Vendor **ASUS**, Hardware Model **ROG STRIX X870-F GAMING WIFI** -- identical to
  Issue #26/#30's own independently measured identity, confirming this is genuinely the same
  physical machine.
- CPU: AMD Ryzen 9 9950X3D (32 logical CPUs / 16 cores, SMT2).
- RAM: 123 GiB total, 15 GiB swap; 110 GiB available at measurement time.
- OS: Ubuntu 24.04.4 LTS, kernel 6.8.0-136-generic.
- GPU: an NVIDIA GeForce RTX 5090 is physically present (`nvidia-smi -L`), but is **not**
  used anywhere in this run -- every tool in this pipeline runs CPU-only.
- Docker: Server Version 29.7.2.
- Nextflow: 26.04.6.12646. Java: OpenJDK 21.0.11.
- Storage: `/home` on `/dev/nvme0n1p3`, 3.5 TB total, 3.1 TB available before this Issue's
  downloads began.

## Current resource contract (confirmed from source, post-Issue #30)

Read directly from `nextflow.config` on this checkout (not from docs), since docs can go
stale relative to code:

| Label | cpus | memory (attempt 1) | time | Consumers |
| --- | ---: | ---: | --- | --- |
| `process_haplotypecaller` | 4 | 16 GiB | 24h | `GATK_HAPLOTYPECALLER` only |
| `process_high` | 8 | 16 GiB | 24h | `GATK_GENOTYPEGVCFS` only (post Issue #30 split) |
| `process_mapping` | 8 | 16 GiB | 24h | `BWA_MEM2_MEM_SORT` |
| `process_bwa_index` | 8 | 16 GiB | 24h | `BWA_MEM2_INDEX` |
| `process_genomicsdb` | 8 | 16 GiB | 24h | `GATK_GENOMICSDBIMPORT` |
| `process_variant_classification` | 2 | 12 GiB | 2h | `CLASSIFY_NORMALIZED_VARIANTS` |
| `process_variant_qc_summary` | 2 | 12 GiB | 2h | `SUMMARIZE_FILTER_QC` |
| `process_gs_panel` | 2 | 8 GiB | 2h | `BUILD_GS_PANEL` |

All `memory` values are `{ N.GB * task.attempt }` closures; `maxRetries = 1`;
`errorStrategy` retries only on exit 137/140/143 (OOM-like), otherwise terminates. At 32
logical CPUs, `process_haplotypecaller` (4 cpus/task) permits up to `floor(32/4) = 8`
concurrent `GATK_HAPLOTYPECALLER` tasks in theory -- whether that is actually reached under
real scheduling/contention is exactly what Stage 1 measures below, not assumed.

## Existing 5-sample provenance (re-confirmed, not re-generated)

Re-verified against Issue #26's own artifacts on `seedcore-01`
(`issue26-cohort-2026-08-23/input/`), read-only (nothing moved, renamed, or deleted):

| Run | BioSample | Bytes (R1+R2) | Published MD5 | Result |
| --- | --- | ---: | --- | --- |
| `SRR29909135` (symlink to Issue #8's copy) | `SAMN42721563` | -- | verified | OK |
| `SRR29909069` | `SAMN42721579` | 1,513,649,942 | verified | OK |
| `SRR29909072` | `SAMN42721577` | 1,893,194,940 | verified | OK |
| `SRR29909067` | `SAMN42721581` | 1,915,571,965 | verified | OK |
| `SRR29909073` | `SAMN42721576` | 1,998,577,713 | verified | OK |

All 10 checksums in `checksums_expected.md5` verified via `md5sum -c` with `OK` for every
entry -- reproduced here as a fact check, not a re-download.

Reference: `GCF_016808095.1_ASM1680809v1_genomic.fna` (Longxiaodou 4 assembly), symlinked
from Issue #8 Phase 5's copy, unchanged.

## Stage 1: 10-sample cohort

### Accession selection

Additional 5 accessions were selected from BioProject `PRJNA1138464`'s remaining 322
WGS/paired-end runs (327 total WGS/paired-end runs minus the 5 already in use; the
BioProject's other 357 runs are RAD-Seq/single-end and out of scope, as established in
Issue #26), sorted by compressed FASTQ size ascending to keep this stage's per-sample
compute/storage/wall-time profile close to the existing cohort's. All candidates have a
unique BioSample (confirmed: no BioSample appears twice among the 327 WGS/paired-end runs)
and ENA-published MD5 checksums.

Approved (human-in-the-loop, before any download):

| Run | BioSample | Library | Bases | Compressed FASTQ | Est. coverage* |
| --- | --- | --- | ---: | ---: | ---: |
| `SRR29909074` | `SAMN42721575` | `AZ633` | 3.94 Gb | 2003 MB | ~8.8x |
| `SRR29909070` | `SAMN42721578` | `AZ639` | 3.89 Gb | 2031 MB | ~8.7x |
| `SRR29909068` | `SAMN42721580` | `AZ659` | 4.29 Gb | 2171 MB | ~9.6x |
| `SRR29909421` | `SAMN42721532` | `AZ662` | 4.33 Gb | 2237 MB | ~9.6x |
| `SRR29909418` | `SAMN42721535` | `AZ665` | 4.52 Gb | 2303 MB | ~10.1x |

\* Same raw-depth estimate methodology as Issue #26 (base_count / 448,362,642 bp).

All five report `library_strategy=WGS`, `library_source=GENOMIC`, `instrument_model=Illumina
HiSeq X` per ENA/SRA metadata -- the same blanket metadata value Issue #26 found to disagree
with real FASTQ read-header instrument-serial patterns for some of its own 5 samples. As in
Issue #26, this document does not assume ENA's `instrument_model` is correct; the real
header-derived identity for each new sample is recorded separately below once downloaded
(section pending download completion).

### Download and checksum verification

All 5 approved accessions' FASTQ (R1+R2) were downloaded from ENA's public FTP mirror and
verified against ENA's own published MD5 checksums. Two files (`SRR29909074_2.fastq.gz`,
`SRR29909068_2.fastq.gz`) hit transient ENA connectivity instability requiring several
retries with resume (`curl -C -`); `SRR29909074_2.fastq.gz` additionally reached its exact
expected byte count once but failed its checksum -- evidence of silent corruption during a
chaotic resume/restart history -- and was re-downloaded clean (no resume, full re-fetch via
`wget`) before it verified correctly. All 10 samples' FASTQ (5 existing + 5 new) verified
`OK` against published MD5 before the pipeline was run. Read-header inspection (whole-file,
not just the first read) confirmed all 5 new samples share the exact same real
instrument/flowcell/lane identity (`E00361:437:H3NC3CCX2:6`) as 3 of the existing 5 samples
-- consistent with the HiSeq 3000/4000-family serial pattern already documented in Issue #26,
still disagreeing with ENA's blanket `instrument_model=Illumina HiSeq X` metadata for every
run in this BioProject.

### Execution

Run via `nextflow run . -profile docker` against commit `f129a01` (the docs-only commit this
branch was on before the resource fix below; identical production code to `main`'s
`a3fe83e`), with the same 10-sample samplesheet, real (non-test-profile) reference config,
and `-with-trace`/`-with-timeline`/`-with-report` enabled, mirroring Issue #26's own
invocation exactly.

- Started: 2026-08-24T22:11:03Z. Completed: 2026-08-25T03:19:18Z (**5h 8m 12s**).
- `NEXTFLOW_RUN_EXIT=0`. 200 total task executions; 199 completed, 1 failed on its first
  attempt and completed on retry (see below) -- zero tasks failed permanently.
- CPU-hours (Nextflow's own report): 98.0.

### HaplotypeCaller concurrency

`process_haplotypecaller` (4 cpus/task, unchanged from Issue #30) permits a theoretical
maximum of `floor(32/4) = 8` concurrent tasks on this host. Computed directly from each of
the 10 real tasks' trace `submit` timestamp and `duration` (interval-overlap analysis, not
assumed from theory): **real max concurrency reached 7 of the theoretical 8** -- very close
to ideal, with only mild scheduling slack. Individual task wall times ranged from 33m58s to
3h36m5s (peak RSS 8.1-12.5 GB, all comfortably within the 16 GiB budget, 22-49% headroom);
as in Issue #26, wall time was not monotonic with sequencing depth, consistent with real
CPU contention among concurrently running tasks rather than a per-task regression.

Host-level `MemAvailable` never dropped below **62.1 GiB** (measured at 30s intervals via a
lightweight monitor for the run's full duration) during the HaplotypeCaller-heavy window --
comfortably healthy, no memory pressure from HaplotypeCaller concurrency itself.

### Downstream resource: a real regression found and fixed

Issue #30's 5-sample-derived memory budgets did **not** hold at 10 samples:

- `CLASSIFY_NORMALIZED_VARIANTS` was **OOM-killed (exit 137) on its 12 GiB first attempt**
  and only completed on its 24 GiB retry (`maxRetries = 1` -- no further retry margin left).
  The host-level monitor shows swap usage spiking to **10.5 GiB** in the exact window of
  that retry (2026-08-25T03:02-03:06Z), strongly suggesting the retry's own "success" was
  itself swap-assisted, not genuine headroom.
- `SUMMARIZE_FILTER_QC (cohort:snp)` and `BUILD_GS_PANEL` both "succeeded" on their first
  attempt, but the trace's own reported peak RSS landed **exactly at each one's memory
  ceiling** (12 GB and 8 GB respectively) -- the identical cgroup-capped/swap-assisted
  false-positive pattern Issue #30 itself identified in `BUILD_GS_PANEL` at 5 samples.

Following the exact same methodology Issue #30 established (reuse the real, already-published
10-sample artifacts -- `cohort_gs.normalized.vcf.gz`, `cohort.snp.filtered.vcf.gz`,
`cohort_gs.snp.pass.vcf.gz` -- replayed directly against the real production scripts inside
the real pinned container, at a deliberately generous 48 GiB memory ceiling so the reported
figure is genuine, unconstrained usage rather than another capped reading), true peak RSS at
10-sample scale was measured as:

| Process | 5-sample true peak (Issue #30) | 10-sample true peak | Previous budget | New budget | New headroom |
| --- | ---: | ---: | ---: | ---: | ---: |
| `CLASSIFY_NORMALIZED_VARIANTS` | 9.33 GiB | **28.47 GiB** | 12 GiB | 40 GiB | 28.8% |
| `SUMMARIZE_FILTER_QC` (snp) | 8.43 GiB | **16.28 GiB** | 12 GiB | 22 GiB | 26.0% |
| `BUILD_GS_PANEL` | 5.34 GiB | **12.65 GiB** | 8 GiB | 17 GiB | 25.6% |

(`SUMMARIZE_FILTER_QC`'s indel invocation remained lightweight at 10 samples too -- 2.6 GB
peak RSS per the real trace, well within budget -- consistent with the 5-sample finding that
indel records are far less numerous than SNP records in this cohort; not separately
re-benchmarked.) Two data points (5 and 10 samples) are not treated as evidence of a linear
scaling formula -- the growth from 5→10 samples was roughly 2.6-3.1x per process, well above
naive 2x sample-count scaling, and this is reported as an observed fact requiring
re-measurement at 20 and 30 samples, not extrapolated forward.

All three benchmark outputs were cross-checked against the real run's own production
outputs and found data-identical (classification/accounting record counts, GS panel matrix
row count, filter-breakdown contents all matched exactly; the only difference was a cosmetic
`cohort_id` label from the benchmark's own invocation choice).

**Fix applied**: `nextflow.config`'s `process_variant_classification`, `process_variant_qc_summary`,
and `process_gs_panel` first-attempt memory values were raised to 40/22/17 GiB respectively
(commit `556f38f`), restoring real positive headroom over the confirmed 10-sample true peaks.
`tests/bin/test_resource_label_contracts.py` was updated to assert the new values. Full
regression (252 unit tests, `nextflow lint .` 35/0, `git diff --check`, the same 26-test
nf-test suite) passed on this commit. The full 10-sample E2E pipeline was **not** re-executed
after this fix -- the change is a memory-budget parameter with no effect on scientific output,
and the benchmark itself (at a ceiling well above the new production values) already confirms
each process completes correctly and produces identical output at the new headroom.

### Joint Genotyping

`GATK_GENOMICSDBIMPORT` (36 intervals) and `GATK_GENOTYPEGVCFS` (36 intervals) both remained
on `process_genomicsdb`/`process_high` (8 cpus, 16 GiB, unchanged -- out of this Issue's
scope) with no OOM/retry on any of the 72 real tasks; peak RSS ranged up to ~1.1 GB
(GenomicsDBImport) and ~3.1 GB (GenotypeGVCFs), both far within budget. `GATK_GATHERVCFS`
completed in 17.5s at 742 MB peak RSS. `genomicsdb_batch_size` (default 50) remains
unexercised in its actual batching behavior at 10 samples, exactly as at 5 -- still below the
default batch size.

### Storage

| Location | Size |
| --- | ---: |
| `input/` (new downloads only; existing 5 are symlinks) | 10.5 GB |
| `runs/10sample/work/` (Nextflow work directory, incl. 36 GenomicsDB workspaces) | 111.4 GB |
| `runs/10sample/results/` (published outputs) | 58.4 GB |
| Total attributable to this run | ~180 GB |
| Host free space before / after | 3.1 TB / 3.1 TB |

Roughly 1.9x Issue #26's 5-sample ~95 GB, consistent with expected linear-ish scaling; no
storage concern at this or the next stage.

### Sample and variant accounting

- `bcftools query -l cohort.raw.vcf.gz`: exactly 10 unique samples
  (`SRR29909067,SRR29909068,SRR29909069,SRR29909070,SRR29909072,SRR29909073,SRR29909074,
  SRR29909135,SRR29909418,SRR29909421`), alphabetical order, matching the samplesheet's
  10 entries with no duplication or loss.
- Raw/all: 10,247,125 records (8,748,064 SNP + 1,422,224 indel per GATK's own type
  classification; 76,837 `records_not_selected` MIXED/MNP/other; 0 `snp_indel_duplicate_records`).
  Ti/Tv 1.85 (vs. 1.87 at 5 samples -- stable). Missing genotypes 12,627,888/102,471,250
  (12.3%, vs. 10.8% at 5 samples -- a modest, expected increase from more samples each
  contributing their own per-site missing calls).
- PASS SNP: 7,780,755. PASS indel: 1,419,193.

### GS panel accounting

`panel_status: populated`. Cross-checked internally by the pipeline's own record accounting
(never assumed to agree): genotype matrix has 7,939,188 variant rows x 10 sample columns;
variant metadata has 7,939,188 rows; sample metadata has 10 rows. All three agree, and the
manifest step itself would have hard-failed on any disagreement.

### 10-sample Gate decision: **GO**

| Criterion | Result |
| --- | --- |
| E2E exit 0 | Yes (with one documented, now-fixed OOM+retry) |
| No unexplained OOM/swap thrashing/host instability | The one OOM/swap event was fully explained (real resource gap) and fixed, not left unexplained |
| HaplotypeCaller concurrency operationally acceptable | Yes -- 7/8 real concurrency, `MemAvailable` never below 62.1 GiB |
| Downstream positive headroom, or resource fix completed via separate commit + regression | Fix completed (commit `556f38f`), full regression green |
| Sample count/order correct | Yes -- 10/10, correct order |
| Variant accounting non-negative and internally consistent | Yes |
| GS panel sample/variant dimensions match metadata | Yes -- cross-checked by the pipeline itself |
| Provenance manifest complete | Yes -- [`docs/real_cohort_scale_validation_10sample_manifest.json`](real_cohort_scale_validation_10sample_manifest.json) |

**Proceeding to Stage 2 (20 samples) requires separate human approval before any additional
accession is downloaded**, per this Issue's own standing rule -- this document's GO verdict
covers Stage 1 only and does not itself authorize Stage 2.

## Stage 2: 20-sample cohort

_(Pending Stage 1 Gate approval. Not executed.)_

## Stage 3: 30-sample cohort

_(Pending Stage 2 Gate approval and a decision that 30 samples adds sufficient information.
Not executed.)_

## Final operational decision

_(Pending all executed stages.)_

## Limitations and 327-sample unresolved risks

_(Pending final stage.)_
