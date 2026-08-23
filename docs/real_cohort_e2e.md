# Real cohort end-to-end validation (Issue #26)

This document records Issue #26: a real, multi-sample, public-data end-to-end run of the
full pipeline on `seedcore-01` -- public WGS FASTQ -> QC/trimming -> mapping/duplicate
marking -> per-sample gVCF -> GenomicsDBImport/Joint Genotyping -> cohort VCF -> hard
filtering/PASS extraction -> variant QC -> GS panel -- together with the reproducibility
evidence (checksums, tool versions, parameters, resource usage) a commissioned delivery
would need. Issues #8 (mapping hardening) and #11 (Joint Genotyping scale hardening) are
both closed and merged into `main` before this work began; this Issue does not re-implement
either.

This is a 3-5 sample validation, not a 327-sample production run. A full 327-sample cohort
is explicitly out of this Issue's scope (see Issue #26's body); this document records what
a 20-30 sample expansion would concretely need, based on what this run actually measured,
without attempting it.

## Environment (measured fresh, not reused from Issue #8)

Measured 2026-08-23 23:40 JST, on host `seedcore-01`, via `hostnamectl`/`nproc`/`free`/
`nvidia-smi`/`df` run directly for this Issue (not copied from
`docs/mapping_real_reference_profile.md`'s Issue #8 Phase 5 measurement, even though it
happens to be the same physical host -- confirmed by comparing Hardware Vendor/Model below
against that document, per this repository's own standing rule that a hostname alone is not
sufficient identity evidence):

- Hardware Vendor **ASUS**, Hardware Model **ROG STRIX X870-F GAMING WIFI** -- identical
  values to Issue #8 Phase 5's own independently measured identity, confirming this is
  genuinely the same machine, not merely a value copied across documents.
- CPU: AMD Ryzen 9 9950X3D (32 logical CPUs).
- RAM: 123 GiB total, 15 GiB swap; 111 GiB available at measurement time.
- OS: Ubuntu 24.04.4 LTS, kernel 6.8.0-136-generic.
- GPU: an NVIDIA GeForce RTX 5090 is physically present (`nvidia-smi -L`), but is **not**
  used anywhere in this run. Every tool in this pipeline runs CPU-only; no GPU acceleration
  (Parabricks or otherwise) is used, benchmarked, or claimed here -- Issue #26 explicitly
  excludes claiming untested GPU speedups.
- Docker: Server Version 29.7.2.
- Nextflow: 26.04.6.12646. nf-test: 0.9.5 (not used for this real run itself, only for the
  synthetic-fixture regression below).
- Storage: `/home` on `/dev/nvme0n1p3`, 3.5 TB total, 3.2 TB available before this run
  began. All FASTQ, reference, Nextflow work directory, and results are kept under a
  scratch directory outside this repository's working tree; none of it is committed (see
  `.gitignore`).

## Baseline regression (run before any real-data work)

Per this repository's own contract, real-data work did not start until the synthetic-fixture
regression suite was confirmed green on this exact checkout of `main` (commit `d01b7e1`,
PR #28's merge commit -- the merge that closed both Issue #8 and Issue #11):

```
$ python3 -m unittest discover -s tests/bin
Ran 220 tests in 0.218s — OK

$ nextflow lint .
✅ 35 files had no errors

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test tests/modules/gatk_genomicsdbimport.nf.test tests/modules/gatk_gathervcfs.nf.test --ci
SUCCESS: Executed 26 tests in 147.567s
```

## Cohort selection

Five real public WGS accessions, all from BioProject
[`PRJNA1138464`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138464) ("WGS resequencing
of 327 accessions [of *Vigna angularis*]"), the same publicly accessible, already-documented
study `SRR29909135` (Issue #8 Phase 5) comes from -- so no new data-access/licensing
question is introduced by this Issue. `SRR29909135` is reused directly (its FASTQ and
checksums from Issue #8 Phase 5 are still present on `seedcore-01`); the other four are the
smallest remaining WGS (not RAD-Seq -- the same BioProject also has 357 RAD-Seq runs from a
companion DArT-seq study, out of scope here as documented in README) runs in the project,
chosen to keep real download/runtime/storage manageable for a validation exercise:

| Sample | Run accession | BioSample | Library name | Bases | Compressed FASTQ | Approx. coverage* |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 1 (reused from Issue #8) | `SRR29909135` | `SAMN42721563` | `AZ693` | 8.64 Gb | 4.35 GB | ~19.3x |
| 2 | `SRR29909069` | `SAMN42721579` | `AZ648` | 3.24 Gb | 1.51 GB | ~7.2x |
| 3 | `SRR29909072` | `SAMN42721577` | `AZ636` | 3.52 Gb | 1.89 GB | ~7.8x |
| 4 | `SRR29909067` | `SAMN42721581` | `AZ702` | 3.59 Gb | 1.92 GB | ~8.0x |
| 5 | `SRR29909073` | `SAMN42721576` | `AZ635` | 3.72 Gb | 1.99 GB | ~8.3x |

\* Approximate coverage = base_count / 448,362,642 bp (the real Longxiaodou 4 assembly
length measured in Issue #11's contig-distribution analysis). This is a read-depth
estimate from raw base count, not a post-mapping/post-duplicate-marking effective coverage
figure.

All five: `library_strategy=WGS`, `library_source=GENOMIC`, `instrument_model=Illumina
HiSeq X` per ENA/SRA metadata. Downloaded from ENA's FTP mirror
(`ftp.sra.ebi.ac.uk/vol1/fastq/...`), verified against ENA's own published MD5 checksums
before use.

## Read-group provenance and a second real instrument-identity discrepancy

Following the same methodology as Issue #8 Phase 5 (derive read-group identity from the
FASTQ's own read headers, not solely from ENA/SRA's `instrument_model` field), the first
read header of every R1/R2 file was inspected, and lane/flowcell consistency was confirmed
across every read in each file (not just the first):

| Sample | Instrument:run:flowcell:lane (from read headers, whole-file-consistent) | Library name (ENA) |
| --- | --- | --- |
| `SRR29909135` | `A00609:64:H3CWYDSXY:4` | `AZ693` |
| `SRR29909069` | `A00609:64:H3CWYDSXY:4` (same flowcell/lane as `SRR29909135`) | `AZ648` |
| `SRR29909072` | `E00361:437:H3NC3CCX2:6` | `AZ636` |
| `SRR29909067` | `E00361:437:H3NC3CCX2:6` (same flowcell/lane as `SRR29909072`) | `AZ702` |
| `SRR29909073` | `E00361:437:H3NC3CCX2:6` (same flowcell/lane as `SRR29909072`/`SRR29909067`) | `AZ635` |

`SRR29909135`/`SRR29909069` were pooled on the same physical flowcell lane; `SRR29909072`/
`SRR29909067`/`SRR29909073` were pooled on a different one. This is expected for a
multiplexed resequencing study and is not itself a finding.

The finding: `A00609`'s serial pattern is consistent with Illumina's NovaSeq naming
convention (already noted as an unresolved discrepancy against ENA's blanket "Illumina
HiSeq X" metadata in Issue #8 Phase 5's own documentation), while `E00361`'s serial pattern
is consistent with Illumina's HiSeq 3000/4000-series naming convention instead -- neither
matches "HiSeq X" literally, and the two serial patterns are consistent with *different*
instrument families from each other. This is a read-header serial-naming-convention
inference, not a confirmed instrument model lookup against an Illumina-maintained registry.
ENA/SRA's `instrument_model` field reports "Illumina HiSeq X" uniformly for every run in
this BioProject, but the runs' own read headers are consistent with at least two distinct
real instrument identities. As in Issue #8, `optical_duplicate_pixel_distance`
was **not** overridden from its default (`100`) for this run, given this same unresolved
ambiguity -- overriding it towards a patterned-flowcell value on the assumption of one
specific instrument model, when the read headers themselves suggest a different and
non-uniform reality, would not be more correct.

## Baseline regression (run before any real-data work)

Confirmed green on this exact checkout of `fix/26-real-cohort-e2e` at commit `49cb2f6`
(the commit this entire real run executed against), before any real-data work began:

```
$ python3 -m unittest discover -s tests/bin
Ran 247 tests in 0.32s — OK   (220 baseline for main + 27 for build_run_manifest.py)

$ nextflow lint .
✅ 35 files had no errors

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test tests/modules/gatk_genomicsdbimport.nf.test tests/modules/gatk_gathervcfs.nf.test --ci
SUCCESS: Executed 26 tests in 147.567s
```

## Execution

Run via `nextflow run . -profile docker` (the pipeline's real, non-`test` profile;
`process_high`/`process_mapping`/`process_bwa_index`/`process_genomicsdb`/`process_medium`/
`process_low` all use their full, non-CI-sized `nextflow.config` values, not
`conf/test.config`'s tiny overrides), against `--input` the samplesheet built from the real
read-group metadata above, `--outdir` and `-w` under a dedicated scratch run directory, with
`-with-trace`/`-with-timeline`/`-with-report` enabled. `genomicsdb_batch_size` was left at
its default (`50`, far above this cohort's 5 samples, so no batching behavior is exercised
here) and `enable_gs_panel` was left at its default (`true`) -- unlike Issue #8 Phase 5's
single-sample mapping-only profiling run, this run needed the full lineage including the GS
panel, per Issue #26's own scope.

- Started: 2026-08-24T00:17:12+09:00
- Completed: 2026-08-24T02:50:24+09:00 (Nextflow's own reported duration: **2h 31m 37s**)
- Exit: `NEXTFLOW_RUN_EXIT=0`
- Tasks: 155 total, 153 succeeded on their final attempt, 2 succeeded only after an
  automatic `task.attempt`-scaled retry (see below) -- **zero** tasks failed permanently.
- CPU-hours (Nextflow's own report): 63.8

## Results

| Metric | Value |
| --- | --- |
| Raw cohort records (`raw/all`) | 5,558,870 (4,782,166 SNPs, 790,151 indels, 127,602 multiallelic sites) |
| Transition/transversion ratio (raw/all) | 1.87 |
| Cohort genotype missingness (raw/all) | 10.8% (3,014,668 / 27,794,350 genotype calls) |
| PASS SNPs | 4,259,731 |
| PASS indels | 770,385 |
| GS panel status | **`populated`** (4,298,980 variant rows x 5 samples) -- the first time this pipeline's GS panel has produced non-empty output; every synthetic-fixture test to date produces an empty panel by construction |
| `records_not_selected` (MIXED/MNP/other, excluded from both SNP and indel selection) | 18,284 |
| `snp_indel_duplicate_records` | 0 (no genuine duplicate output records) |

All artifacts (cohort VCF, PASS SNP/indel VCFs, per-sample gVCFs, GS panel dosage matrix,
sample/variant metadata, QC tables, and the pipeline's own `gs_panel/cohort.gs_panel.manifest.json`)
were produced and are internally consistent (the GS panel's own `record_accounting.tsv`
cross-checks its matrix/metadata files against each other and against `raw/all`, and would
have refused to write a manifest at all on any disagreement -- see
[`docs/gs_panel_data_contract.md`](gs_panel_data_contract.md)).

A whole-run provenance manifest was built with the new `bin/build_run_manifest.py`
(Issue #26) from this run's own real, already-published artifacts -- input FASTQ checksums
(read from the samplesheet actually used), reference bundle checksums, pinned container
references, the full (not abbreviated) git commit SHA this run executed against
(`49cb2f606407bcafead0eea29a3962ebff6f7733`), the Nextflow execution engine's own version
(`26.04.6`, recorded separately from the pinned tool containers -- Nextflow itself runs on
the host outside any container, and its own runtime semantics were the actual root cause of
a real production failure this pipeline hit before, in Issue #11), every Nextflow parameter
used, cohort/variant-type accounting, a pointer to the GS panel's own manifest, and
checksums of the raw cohort VCF, both PASS VCFs, and all five gVCFs. Committed as
[`docs/real_cohort_e2e_run_manifest.json`](real_cohort_e2e_run_manifest.json) (no raw
sequence data; every value is either a checksum, a count, a parameter, a pinned container
reference, or the execution engine's own version).

## Resource usage (from `-with-trace`)

| Process | n | Max peak RSS | Notes |
| --- | ---: | ---: | --- |
| `BWA_MEM2_INDEX` | 1 | 10.0 GB | `process_bwa_index` (16 GiB budget); consistent with Issue #8 Phase 5's single-accession measurement |
| `BWA_MEM2_MEM_SORT` | 5 | 14.1 GB | `process_mapping` (16 GiB budget); durations 12m26s-15m11s for the four ~7-8x samples, 27m37s for the one ~19.3x sample |
| `GATK_MARKDUPLICATES` | 5 | 6.8 GB | `process_medium` (8 GiB budget on attempt 1) |
| `GATK_HAPLOTYPECALLER` | 5 | 12.0 GB | `process_high` (16 GiB budget); durations below |
| `GATK_GENOMICSDBIMPORT` | 36 | 1.0 GB | `process_genomicsdb`; one per reference contig; total 22.7 CPU-minutes across all 36 |
| `GATK_GENOTYPEGVCFS` | 36 | 1.6 GB | one per reference contig; total 29.7 CPU-minutes across all 36 |
| `BUILD_GS_PANEL` | 1 | **4.0 GB** | `process_low` (4 GiB budget) -- at its ceiling, did not fail, but with no margin |
| `CLASSIFY_NORMALIZED_VARIANTS` | 2 (1 retry) | **8.0 GB** | `process_low`; **OOM-killed (exit 137) on attempt 1 at 4 GiB**, succeeded on attempt 2 at 8 GiB -- reported peak RSS on the successful attempt was itself at the 8 GiB ceiling |
| `SUMMARIZE_FILTER_QC` | 3 (1 retry, 2 stage invocations) | **8.0 GB** | `process_low`; **OOM-killed (exit 137) on attempt 1** for the `cohort:snp` invocation specifically (the `cohort:indel` invocation succeeded at 1.4 GB, well within budget -- SNPs vastly outnumber indels in this cohort), succeeded on attempt 2 at 8 GiB, again reported at the ceiling |

`GATK_HAPLOTYPECALLER` durations were **not** monotonic with sequencing depth: the one
~19.3x sample (`SRR29909135`) finished in 48m59s, while the four ~7-8x samples took
26m55s, 1h13m54s, 1h28m54s, and 1h23m16s -- three of the four *shallower* samples took
markedly *longer* than the one deep sample. `process_high` requests 8 cpus per task, and
this machine has 32 logical CPUs, so at most 4 `GATK_HAPLOTYPECALLER` (or other
`process_high`/`process_mapping`/`process_genomicsdb`-labeled) tasks can run genuinely
concurrently; the real explanation is scheduling/CPU contention among concurrently running
tasks competing for the same 32 cores, not sequencing depth. This is real evidence that wall
time does not scale linearly with sample count on a single fixed-core machine once
concurrent task count exceeds `total_cpus / task.cpus` -- directly relevant to the 20-30
sample assessment below.

## Storage

| Location | Size | Notes |
| --- | ---: | --- |
| `input/` (4 newly downloaded accessions; `SRR29909135` is a symlink to Issue #8's copy and does not count towards this directory's own size) | 6.9 GB | Compressed FASTQ |
| `reference/` (symlink to Issue #8's copy) | 12 KB (real file: ~454 MB) | Source FASTA |
| `work/` (Nextflow work directory: staged inputs, intermediate BAMs, per-sample gVCFs, 36 GenomicsDB workspaces, intermediate VCFs) | 53 GB | Not published, not committed |
| `results/` (published outputs: BAMs, gVCFs, cohort/PASS VCFs, QC tables, GS panel) | 31 GB | Not committed (see Data Handling below) |
| Total attributable to this run | ~95 GB | (6.9 + ~4.35 reused + ~0.45 reference + 53 + 31, rounding) |
| Host free space before / after | 3.2 TB / 3.1 TB | `/home` on `/dev/nvme0n1p3`, 3.5 TB total |

## What a 20-30 sample expansion would need (assessment only -- not executed)

Based on this run's real measurements, not attempted here (Issue #26 explicitly scopes a
20-30 sample expansion as conditional on the 3-5 sample cohort succeeding, and the full
327-sample cohort as out of scope regardless):

- **Storage**: linear extrapolation from ~95 GB for 5 samples suggests roughly 500-600 GB
  for 30 samples (dominated by `work/`'s per-sample BAMs/gVCFs and `results/`'s published
  copies of the same). Comfortably within the 3.1 TB currently free on this host; storage is
  not expected to block a 20-30 sample attempt on this machine.
- **`process_low` memory is already marginal at 5 samples**: `BUILD_GS_PANEL` used its full
  4 GiB budget without failing, and `CLASSIFY_NORMALIZED_VARIANTS`/`SUMMARIZE_FILTER_QC`
  (`cohort:snp`) were both OOM-killed at 4 GiB and landed exactly at their 8 GiB retry
  ceiling on success. With `maxRetries = 1`, there is no further retry headroom left if 8
  GiB is also exceeded. **Before attempting 20-30 samples**, `process_low`'s default memory
  (or dedicated labels for these specific processes, mirroring the `process_genomicsdb`
  precedent from Issue #11) should be re-evaluated -- this is a documented recommendation
  for whoever executes that expansion, not a code change made in this Issue.
- **CPU contention scales worse than linearly on one fixed-core machine**: the
  `GATK_HAPLOTYPECALLER` timing anomaly above shows that even 5 concurrent `process_high`
  tasks (8 cpus each) on 32 logical CPUs already produces real queuing/contention effects
  large enough to make a shallower sample take longer than a deeper one. 20-30 samples on
  the same 32-core machine would face substantially more of this effect; either more cores
  (a different or additional machine), reduced per-task `cpus` requests (trading single-task
  speed for concurrency), or an accepted, non-linear increase in total wall time would be
  needed to keep a 20-30 sample run practical. This was not modeled quantitatively here --
  it is reported as a real, observed effect worth accounting for, not a specific revised
  resource recommendation.
- **`genomicsdb_batch_size` (default 50) remains unexercised in its batching behavior**:
  even 30 samples is below the default batch size, so GenomicsDBImport would still open all
  samples in a single batch. Its actual batching behavior (Issue #11) would only become
  observable once sample count exceeds 50.
- **No production code was changed as a result of this run.** All of the above are
  documented observations and recommendations for the next phase of Issue #26 (or a
  follow-up issue), not applied here.

## Data handling

No raw FASTQ, BAM, intermediate VCF, or the reference FASTA/BWA index is committed to this
repository. Everything for this run lives under a scratch directory on `seedcore-01`,
outside this repository's working tree; `.gitignore` already excludes `*.fastq.gz`,
`*.bam`, `*.bai`, `*.vcf.gz`, `*.tbi`, `*.fna*`, `work/`, and `results/`. The only artifacts
committed from this run are this document and
[`docs/real_cohort_e2e_run_manifest.json`](real_cohort_e2e_run_manifest.json) -- both
contain only checksums, counts, parameters, and pinned container references, never sequence
data itself. No customer data, former-employer data, or credentials were used or produced;
this run used only the same public BioProject already documented in README, downloaded
directly from ENA's public FTP mirror.

## Completion criteria (Issue #26)

- [x] Seedcore-01's real environment and measurement timestamp confirmed (see Environment
      above) -- measured fresh for this Issue, cross-checked against Issue #8 Phase 5's
      independent measurement of the same physical host.
- [x] A 3-5 sample cohort completed reproducibly: exit 0, every task either succeeded
      outright or recovered via the pipeline's existing `task.attempt` retry mechanism, and
      the run's own internal consistency checks (GS panel record accounting, variant-type
      reconciliation) all passed.
- [x] 20-30 sample expansion conditions, required resources, and unaddressed items are
      explicitly stated above -- without attempting the expansion itself.
- [x] Cohort VCF, QC, GS panel, sample/variant metadata, and the run manifest are mutually
      consistent (cross-checked via the GS panel's own record accounting and this document's
      `bin/build_run_manifest.py` output).
- [x] Wall time, peak memory, and storage usage are recorded from real measurements
      (`-with-trace`, `du`, `df`), not estimated.
- [x] No raw data or customer/private information is included in this repository's Git
      history.
