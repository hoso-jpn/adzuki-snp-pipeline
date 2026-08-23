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
before use (see Provenance below).
