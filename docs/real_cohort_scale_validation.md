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

_(Sections below -- download/checksum verification, execution, HaplotypeCaller concurrency,
downstream resource re-measurement, Joint Genotyping, storage, accounting, Gate decision --
are filled in as Stage 1 actually executes.)_

## Stage 2: 20-sample cohort

_(Pending Stage 1 Gate approval. Not executed.)_

## Stage 3: 30-sample cohort

_(Pending Stage 2 Gate approval and a decision that 30 samples adds sufficient information.
Not executed.)_

## Final operational decision

_(Pending all executed stages.)_

## Limitations and 327-sample unresolved risks

_(Pending final stage.)_
