# adzuki-snp-pipeline


This repository provides a reproducible SNP calling workflow for adzuki bean (*Vigna angularis*) using publicly available whole-genome sequencing datasets.

The project originated from research activities in plant genetics and genomic prediction, and is maintained as part of the broader Florigen AI research ecosystem.

---

## Data Source

| Item | Detail |
|---|---|
| BioProject | [PRJNA1138464](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138464) |
| Publication | Chien et al. 2025, *Science* 388: eads2871 |
| Data type | WGS resequencing (327 accessions) + DArT-seq (357 accessions) |
| Reference genome | *Vigna angularis* GCF_016808095.1 (ASM1680809v1) |

---

## Pipeline Workflow

![Workflow](docs/workflow.png)

---

## Environment

```bash
conda create -n bioinfo -c conda-forge -c bioconda \
  fastqc sra-tools bwa samtools gatk4 bcftools ncbi-datasets-cli -y
conda activate bioinfo
```

---

## Usage

### 1. Download reference genome

```bash
datasets download genome accession GCF_016808095.1 \
  --include genome --filename adzuki_reference.zip
unzip adzuki_reference.zip -d adzuki_reference
REF=adzuki_reference/ncbi_dataset/data/GCF_016808095.1/GCF_016808095.1_ASM1680809v1_genomic.fna
```

### 2. Index reference genome

```bash
bwa index $REF
gatk CreateSequenceDictionary -R $REF
samtools faidx $REF
```

### 3. Download WGS data

```bash
prefetch SRR29909135 --output-directory ./raw_data
fasterq-dump ./raw_data/SRR29909135/SRR29909135.sra \
  --outdir ./raw_data --threads 4 --progress
```

### 4. Quality control

```bash
fastqc ./raw_data/SRR29909135_1.fastq \
       ./raw_data/SRR29909135_2.fastq \
       --outdir ./fastqc_results --threads 4
```

### 5. Read mapping

```bash
bwa mem -t 4 \
  -R "@RG\tID:SRR29909135\tSM:SRR29909135\tPL:ILLUMINA" \
  $REF \
  ./raw_data/SRR29909135_1.fastq \
  ./raw_data/SRR29909135_2.fastq \
  | samtools sort -@ 4 -o SRR29909135.bam
samtools index SRR29909135.bam
```

### 6. Duplicate marking

> **Note:** For patterned flow cells (HiSeq X, NovaSeq), use `-d 2500` for optical duplicate handling.

```bash
samtools sort -n -@ 4 SRR29909135.bam -o SRR29909135.namesort.bam
samtools fixmate -m SRR29909135.namesort.bam SRR29909135.fixmate.bam
samtools sort -@ 4 SRR29909135.fixmate.bam -o SRR29909135.fixmate.sort.bam
samtools markdup -r -d 2500 -@ 4 \
  SRR29909135.fixmate.sort.bam \
  SRR29909135.markdup.bam
samtools index SRR29909135.markdup.bam
```

### 7. SNP calling (GVCF mode)

```bash
gatk HaplotypeCaller \
  -R $REF \
  -I SRR29909135.markdup.bam \
  -O SRR29909135.g.vcf.gz \
  -ERC GVCF \
  --tmp-dir /tmp \
  -native-pair-hmm-threads 4
```

### 8. Genotyping

```bash
gatk GenotypeGVCFs \
  -R $REF \
  -V SRR29909135.g.vcf.gz \
  -O SRR29909135.vcf.gz
```

> **Note for multi-sample analysis:** Use `GenomicsDBImport` + `GenotypeGVCFs` (Joint Genotyping) for 30+ accessions.

### 9. Variant filtration

```bash
gatk SelectVariants \
  -V SRR29909135.vcf.gz \
  --select-type-to-include SNP \
  -O SRR29909135.snp.vcf.gz

gatk VariantFiltration \
  -V SRR29909135.snp.vcf.gz \
  --filter-expression "QD < 2.0"    --filter-name "QD2" \
  --filter-expression "FS > 60.0"   --filter-name "FS60" \
  --filter-expression "MQ < 40.0"   --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O SRR29909135.snp.filtered.vcf.gz

gatk SelectVariants \
  -V SRR29909135.snp.filtered.vcf.gz \
  --exclude-filtered \
  -O SRR29909135.snp.pass.vcf.gz
```

---

## Example Results

![Mapping Stats](docs/mapping_stats.png)

| Step | Count | Note |
|---|---|---|
| Total reads | 57,597,756 | Paired-end, 150 bp |
| Mapping rate | **99.35%** | BWA-MEM to GCF_016808095.1 |
| Properly paired | 93.59% | |
| Duplicate rate | ~10% | HiSeq X patterned flow cell |
| Total variants | 783,836 | SNP + indel |
| SNPs extracted | 678,212 | After SelectVariants |
| High-quality SNPs | **610,790** | After hard filtering (PASS) |

![Variant Counts](docs/variant_counts.png)

---

## Notes on Adzuki Bean Genomics

- **RAD-seq is not recommended** for adzuki bean. The genome (~540 Mb) contains abundant repetitive sequences. Low-coverage WGS is strongly preferred.
- For multi-sample GS panel construction, Joint Genotyping across all accessions is strongly recommended.
- The reference genome GCF_016808095.1 is a chromosome-level assembly suitable for population-scale SNP calling.

---

## Related Repositories

- [adzuki-gwas-analysis](https://github.com/hoso-jpn/adzuki-gwas-analysis) — GWAS for water absorption traits and FT gene exploration *(in progress)*
- [genomic-prediction-benchmark](https://github.com/hoso-jpn/genomic-prediction-benchmark) — GS model benchmarking with SoyNAM *(in progress)*
- [genomic-prediction-resnet-hybrid](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid) — ResNet + linear hybrid GS model

---

## Author

**Hoso**
Founder, [Florigen AI](https://florigen.ai)
Plant Genetics × Bioinformatics × Physical AI

- GitHub: [@hoso-jpn](https://github.com/hoso-jpn)
- Researchmap: [hosokawa-yusuke](https://researchmap.jp/hosokawa-yusuke)

---

## License

MIT License
