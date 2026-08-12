# adzuki-snp-pipeline

Research and development repository for building a reproducible SNP-calling pipeline for adzuki bean (*Vigna angularis*) from publicly available whole-genome sequencing data.

The repository currently documents a manually executed, single-sample pilot analysis. It does **not yet contain an executable workflow or a fully reproducible production pipeline**. A Nextflow DSL2 implementation, automated tests, and multi-sample cohort validation are planned in [Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1).

This work represents plant-genetics and bioinformatics research that informs the longer-term agricultural AI activities of Florigen AI. It does not imply a direct genomic-prediction-to-Physical-AI development path.

---

## Current Maturity

| Capability | Status | Evidence or limitation |
| --- | --- | --- |
| Manual single-sample SNP-calling pilot | Executed once | SRR29909135 was processed manually; the result has not yet been reproduced by an automated test |
| Documented command sequence | Available | Commands are recorded below, but software versions and execution parameters are not yet fully locked |
| Automated workflow | Not implemented | Nextflow DSL2 implementation is planned |
| Configurable reference bundle | Not implemented | The current instructions use GCF_016808095.1 |
| Multi-sample Joint Genotyping | Not validated | GenomicsDBImport-based cohort processing is planned |
| Trimming and integrated QC report | Not implemented | fastp and MultiQC are planned |
| Pipeline-level tests | Not implemented | nf-test and a synthetic test dataset are planned |
| Functional CI | Not implemented | Current CI checks repository structure only |
| Production use | Not supported | This is an experimental plant-research repository |

The figures and variant counts in this README are historical results from the single-sample pilot. They are not yet backed by an automated clean-environment reproduction test.

---

## Data and Reference Sources

### Sequencing data

| Item | Detail |
| --- | --- |
| BioProject | [PRJNA1138464](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138464) |
| Demonstration accession | SRR29909135 |
| Associated publication | [Chien et al. 2025, *Science* 388: eads2871](https://doi.org/10.1126/science.ads2871) |
| Public data described by the study | WGS resequencing of 327 accessions and DArT-seq data from 357 accessions |
| Scope validated in this repository | One WGS accession only |

The associated study also reports DArT-seq data from 357 accessions. The current and planned variant-calling scope of this repository is WGS; no RAD-seq or DArT-seq workflow is implemented or evaluated here.

### Reference genome

The manual pilot used the following independent reference assembly:

| Item | Detail |
| --- | --- |
| Accession | [GCF_016808095.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016808095.1/) |
| Assembly | ASM1680809v1 |
| Cultivar | Longxiaodou 4 |
| Assembly span | 447.8 Mb |
| Publication | [Li et al. 2024, *Scientific Data* 11:1074](https://doi.org/10.1038/s41597-024-03911-y) |

The approximately 540 Mb figure often cited for adzuki bean is a k-mer-based genome-size estimate for the cultivar Shumari ([Sakai et al. 2015, *Scientific Reports* 5:16780](https://doi.org/10.1038/srep16780)). For Longxiaodou 4, Li et al. estimated a genome size of 464.9 Mb by 21-mer analysis; the 447.8 Mb assembly represents 96.32% of that estimate. These values differ by cultivar and estimation method and should not be treated as interchangeable.

The sequencing data and reference assembly originate from different studies and genetic backgrounds. Future pipeline versions will treat the reference genome as an explicit, configurable analysis input rather than assuming that results are interchangeable across references.

---

## Pilot Workflow

![Manual pilot workflow](docs/workflow.png)

The following commands document the historical single-sample pilot. They are retained as technical evidence and as input to the planned Nextflow implementation.

They are **not yet a clean-environment reproduction procedure** because dependency versions, checksums, resource requirements, and all intermediate validation steps have not been fixed.

For readability and shell correctness, the recorded commands have been lightly edited to add output-directory creation, quote the reference path, and correct the long-option spelling for `--native-pair-hmm-threads`. These edits do not indicate that the complete procedure was rerun for this documentation update.

---

## Historical Environment Setup

```bash
conda create -n bioinfo -c conda-forge -c bioconda \
  fastqc sra-tools bwa samtools gatk4 bcftools ncbi-datasets-cli -y
conda activate bioinfo
```

---

## Historical Single-Sample Procedure

### 1. Download the reference genome

```bash
datasets download genome accession GCF_016808095.1 \
  --include genome \
  --filename adzuki_reference.zip

unzip adzuki_reference.zip -d adzuki_reference

REF=adzuki_reference/ncbi_dataset/data/GCF_016808095.1/GCF_016808095.1_ASM1680809v1_genomic.fna
```

### 2. Index the reference genome

```bash
bwa index "$REF"
gatk CreateSequenceDictionary -R "$REF"
samtools faidx "$REF"
```

### 3. Download the demonstration WGS data

```bash
prefetch SRR29909135 --output-directory ./raw_data

fasterq-dump ./raw_data/SRR29909135/SRR29909135.sra \
  --outdir ./raw_data \
  --threads 4 \
  --progress
```

### 4. Perform initial quality control

```bash
mkdir -p fastqc_results

fastqc \
  ./raw_data/SRR29909135_1.fastq \
  ./raw_data/SRR29909135_2.fastq \
  --outdir ./fastqc_results \
  --threads 4
```

This historical procedure performs inspection with FastQC but does not perform read trimming. fastp-based trimming and reporting are planned for the executable pipeline.

### 5. Map reads

```bash
bwa mem -t 4 \
  -R "@RG\tID:SRR29909135\tSM:SRR29909135\tPL:ILLUMINA" \
  "$REF" \
  ./raw_data/SRR29909135_1.fastq \
  ./raw_data/SRR29909135_2.fastq \
  | samtools sort -@ 4 -o SRR29909135.bam

samtools index SRR29909135.bam
```

### 6. Remove duplicates in the historical pilot

The historical command used `samtools markdup -r -d 2500`. The `-r` option removed duplicate reads rather than only setting the duplicate flag. This irreversible behavior is retained here solely as a record of the pilot procedure. The planned pipeline will omit `-r`, preserve duplicate records with the DUP flag, and allow downstream tools to exclude them.

```bash
samtools sort -n -@ 4 \
  SRR29909135.bam \
  -o SRR29909135.namesort.bam

samtools fixmate -m \
  SRR29909135.namesort.bam \
  SRR29909135.fixmate.bam

samtools sort -@ 4 \
  SRR29909135.fixmate.bam \
  -o SRR29909135.fixmate.sort.bam

samtools markdup -r -d 2500 -@ 4 \
  SRR29909135.fixmate.sort.bam \
  SRR29909135.markdup.bam

samtools index SRR29909135.markdup.bam
```

### 7. Call variants in GVCF mode

```bash
gatk HaplotypeCaller \
  -R "$REF" \
  -I SRR29909135.markdup.bam \
  -O SRR29909135.g.vcf.gz \
  -ERC GVCF \
  --tmp-dir /tmp \
  --native-pair-hmm-threads 4
```

### 8. Genotype the demonstration sample

```bash
gatk GenotypeGVCFs \
  -R "$REF" \
  -V SRR29909135.g.vcf.gz \
  -O SRR29909135.vcf.gz
```

This command was used only for the single-sample pilot. The planned multi-sample pipeline will use HaplotypeCaller GVCFs followed by GenomicsDBImport and GenotypeGVCFs for Joint Genotyping. Joint Genotyping will be the default for a multi-sample cohort rather than being conditional on a 30-sample threshold.

### BQSR policy

BQSR is intentionally excluded from the current procedure and from the planned default pipeline. GATK BaseRecalibrator uses a known-sites VCF to distinguish known polymorphisms from mismatches used to model sequencing errors. This repository has not identified or validated an appropriate known-sites resource for the Longxiaodou 4 reference bundle.

Bootstrapped BQSR is also outside the current scope because its known-sites construction and effect on the resulting calls have not been validated here. This is a deliberate design decision rather than an omitted implementation step. See the [GATK BaseRecalibrator documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator).

### 9. Select and filter SNPs

```bash
gatk SelectVariants \
  -V SRR29909135.vcf.gz \
  --select-type-to-include SNP \
  -O SRR29909135.snp.vcf.gz

gatk VariantFiltration \
  -V SRR29909135.snp.vcf.gz \
  --filter-expression "QD < 2.0" \
  --filter-name "QD2" \
  --filter-expression "FS > 60.0" \
  --filter-name "FS60" \
  --filter-expression "MQ < 40.0" \
  --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" \
  --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "ReadPosRankSum-8" \
  -O SRR29909135.snp.filtered.vcf.gz

gatk SelectVariants \
  -V SRR29909135.snp.filtered.vcf.gz \
  --exclude-filtered \
  -O SRR29909135.snp.pass.vcf.gz
```

The hard-filter thresholds above document the historical pilot. Their suitability across references, accessions, and cohort sizes has not yet been validated.

---

## Single-Sample Pilot Results

![Mapping statistics from the single-sample pilot](docs/mapping_stats.png)

| Step | Count | Note |
| --- | ---: | --- |
| Total reads | 57,597,756 | Paired-end, 150 bp |
| Mapping rate | 99.35% | BWA-MEM against GCF_016808095.1 |
| Properly paired | 93.59% | Single-sample pilot |
| Duplicate rate | ~10% | Historical run |
| Total variants | 783,836 | SNPs and indels |
| SNPs extracted | 678,212 | After SelectVariants |
| PASS SNPs | 610,790 | After the documented hard filters |

![Variant counts from the single-sample pilot](docs/variant_counts.png)

![SNP summary from the single-sample pilot](docs/snp_summary.png)

These values describe one historical execution for SRR29909135. They are not estimates for the complete 327-accession cohort and do not demonstrate SNP-calling accuracy or superiority over another workflow.

---

## Planned Pipeline

The implementation roadmap is tracked in [Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1).

The intended pipeline includes:

- Nextflow DSL2 with strict-syntax-compatible code
- samplesheet and schema validation
- configurable reference bundles
- fastp, FastQC, and MultiQC
- BWA-MEM2 mapping and duplicate marking without removing duplicate records
- GATK HaplotypeCaller GVCF generation
- multi-sample Joint Genotyping
- variant and cohort QC
- documented GS-panel output contracts
- nf-test coverage and a synthetic CI dataset
- software-version, parameter, and checksum manifests

These capabilities remain planned until their corresponding issues are implemented and validated.

---

## Data Handling Policy

This public repository is limited to:

- workflow and pipeline source code
- commands and configuration required for reproducibility
- analyses based on publicly available datasets
- synthetic or redistributable test data
- technical validation records that contain no confidential information

The following are not included:

- former-employer or other non-public research data
- customer data
- proprietary SNP panels or marker selections
- confidential business logic
- credentials, access tokens, or private infrastructure details
- material intended to remain protected as future intellectual property

Raw WGS data and full reference bundles are not committed to this repository. They must be obtained from their authoritative public sources.

---

## Related Repositories

The following repositories represent distinct stages of a longer-term research stack. They are not yet connected by an automated end-to-end workflow.

- [adzuki-snp-pipeline](https://github.com/hoso-jpn/adzuki-snp-pipeline) — public WGS data to cohort variants and SNP matrices; currently under reconstruction
- [adzuki-gwas-analysis](https://github.com/hoso-jpn/adzuki-gwas-analysis) — analysis of publicly available GWAS summary statistics
- [genomic-prediction-resnet-hybrid](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid) — auditable comparison of GBLUP and neural genomic-prediction models using compatible individual-level data

The public Dryad dataset currently used by `adzuki-gwas-analysis` does not contain individual-level genotypes and phenotypes. Therefore, these three repositories should not be read as an already-operational SNP-to-GWAS-to-genomic-prediction service.

---

## Author

**Yusuke Hosokawa**<br>
Independent researcher and AI engineer<br>
Building [Florigen AI](https://florigen.ai), a long-term agricultural AI initiative

Plant Genetics × Edge AI × Physical AI

- [GitHub](https://github.com/hoso-jpn)
- [researchmap](https://researchmap.jp/hosokawa-yusuke)
- [Breeding Science (2025) — QTL mapping in rice chromosome segment substitution lines](https://doi.org/10.1270/jsbbs.24058)

The genomics work in this repository documents plant-domain and bioinformatics expertise that informs longer-term agricultural AI and robotics research.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
