# adzuki-snp-pipeline

公開WGSデータからアズキ（*Vigna angularis*）のコホートvariantとgenomic selection（GS）向けSNP dosage matrixを生成する、Nextflow DSL2ベースの研究開発パイプラインです。

現在の`main`は、paired-end FASTQのQC・trimmingからBWA-MEM2 mapping、duplicate marking、GATK HaplotypeCaller gVCF、GenomicsDBImport / GenotypeGVCFsによるJoint Genotyping、hard filtering、variant QC、GS panel、MultiQCまでを一つの実行可能workflowとして実装しています。synthetic fixtureでCIを自動化し、公開WGSでは5→10→20検体までend-to-end実行しています。

> **位置づけ**: 実験的な植物ゲノミクス研究リポジトリです。20検体までの実データ動作は確認していますが、327検体scale、本番SLA、variant call精度の優位性、hard-filter閾値の生物学的最適性を保証するものではありません。

## 現在の到達点

| 項目 | 状態 | 根拠・制約 |
| --- | --- | --- |
| Nextflow DSL2 E2E | 実装済み | strict syntax、parameter / samplesheet validation、digest-pinned containers |
| QC / trimming / mapping | 実装済み | raw/trimmed FastQC、fastp、BWA-MEM2→samtools sort pipe、MarkDuplicates、SAMtools QC |
| Joint Genotyping | 実装済み | HaplotypeCaller gVCF → GenomicsDBImport → GenotypeGVCFs → GatherVcfs |
| Variant filtering / QC | 実装済み | SNP/indel hard filtering、PASS抽出、bcftools stats、FILTER/annotation coverage、record accounting |
| GS panel | 実装済み | normalization、split後再分類、独立再filter、`-1/0/+1/nan` dosage matrix、metadata、accounting、manifest |
| MultiQC | 実装済み | MultiQC 1.35。raw/trimmed FastQC ZIP、fastp JSON、MarkDuplicates、SAMtools QCを明示allowlistで集約 |
| Synthetic CI | 実装済み | Python unit tests、Nextflow lint、pipeline/module-level nf-test、GitHub Actions |
| Real cohort | 20検体まで検証済み | Issue #26で5検体、Issue #33で10→20検体。30検体は追加情報量とコストを比較し明示的にskip |
| Classifier memory | streaming化済み | Issue #35。10/20検体targeted replayでPython RSS約20.8/20.9 MiB、swap delta 0 |
| 327検体 | 未検証 | 50検体超のJoint Genotyping設計・benchmarkを先に実施する |
| Hard-filter biological calibration | 未確立 | 現在値は設定可能な運用上の出発点。truth setなしに最適性を主張しない |
| Production use | 非対応 | research / engineering validation repository |

## Pipeline

```text
paired-end FASTQ
  ├─ FastQC (raw)
  ├─ fastp
  ├─ FastQC (trimmed)
  └─ BWA-MEM2 | samtools sort
       └─ sample-level merge
            └─ GATK MarkDuplicates
                 ├─ SAMtools flagstat / stats / idxstats
                 └─ GATK HaplotypeCaller (gVCF)
                      └─ GenomicsDBImport
                           └─ GenotypeGVCFs
                                └─ GatherVcfs
                                     ├─ SNP / indel selection
                                     │    └─ VariantFiltration
                                     │         ├─ PASS VCF
                                     │         └─ variant / FILTER / annotation QC
                                     └─ GS lineage
                                          └─ bcftools norm -m-
                                               └─ split-record classification
                                                    └─ GS-specific re-filter
                                                         └─ dosage matrix + metadata + accounting + manifest

QC artifacts ───────────────────────────────────────────────→ MultiQC
```

## Quick start

必要要件はJava 17以降、Nextflow 26.04.6、Dockerです。現在pinしているproduction containersはLinux AMD64を主対象としています。

```bash
NXF_VER=26.04.6 nextflow lint .

NXF_VER=26.04.6 \
  nextflow run . \
  -profile test,docker
```

Apple Siliconでsynthetic functional testを行う場合はDocker emulation profileを使用します。これはperformance benchmark用ではありません。

```bash
NXF_VER=26.04.6 \
  nextflow run . \
  -profile test,docker_amd64
```

## 入力契約

### Samplesheet

CSVで、以下の列を使用します。

| 列 | 必須 | 意味 |
| --- | --- | --- |
| `sample_id` | Yes | biological sample ID。同一sampleの複数read groupを許容 |
| `read_group_id` | Yes | samplesheet全体で一意 |
| `fastq_1` | Yes | paired-end R1 |
| `fastq_2` | Yes | paired-end R2 |
| `library_id` | Yes | library ID |
| `platform` | Yes | 現行contractは`ILLUMINA` |
| `platform_unit` | No | flowcell / lane / barcode等 |

必須値欠損、存在しないFASTQ、R1/R2同一path、重複`read_group_id`、FASTQのread-group間再利用、想定外列は解析process開始前に拒否します。

### Reference bundle

`reference_id`、`reference_name`、`reference_fasta`を必須とし、`reference_accession`、`reference_species`、`reference_cultivar`、prebuilt FAI / sequence dictionary / BWA-MEM2 indexを任意指定できます。index未指定時はworkflow内で生成します。

Longxiaodou 4の設定例は[`conf/references/longxiaodou4.config.example`](conf/references/longxiaodou4.config.example)です。実データ検証ではGCF_016808095.1 / ASM1680809v1を使用しました。参照assemblyの異なる結果を直接互換とみなしてはいけません。

## 主要パラメータ

| Parameter | Default | Contract |
| --- | ---: | --- |
| `sample_ploidy` | `2` | run全体のglobal ploidy。HaplotypeCaller / GenotypeGVCFsへ明示伝播 |
| `enable_gs_panel` | `true` | GS lineageのON/OFF。GS schema v1はdiploid-only |
| `optical_duplicate_pixel_distance` | `100` | GATK MarkDuplicatesへ伝播。flowcell特性に応じて設定 |
| `genomicsdb_batch_size` | `50` | GenomicsDBImportへ伝播。50+ sampleでの実batching性能は未検証 |

SNP / indel hard-filter thresholdもparameter化されています。現在値はGATK系の一般的なhard-filter値を出発点としたもので、Longxiaodou 4 / PRJNA1138464 / アズキ集団に対する最適値ではありません。

## 主要出力

```text
<outdir>/
├── qc/
│   ├── fastqc/
│   ├── fastp/
│   ├── alignment/
│   ├── variants/
│   └── multiqc/
│       ├── multiqc_report.html
│       ├── multiqc_config.yaml
│       ├── multiqc_version.txt
│       └── multiqc_data/
├── alignment/
├── variants/
│   ├── gvcf/
│   ├── raw/
│   ├── by_type/
│   ├── filtered/
│   ├── pass/
│   ├── gs_normalized/
│   ├── gs_classified/
│   ├── gs_filtered/
│   └── gs_pass/
└── gs_panel/
    ├── cohort.gs_panel.genotype_matrix.tsv.gz
    ├── cohort.gs_panel.sample_metadata.tsv
    ├── cohort.gs_panel.variant_metadata.tsv
    ├── cohort.gs_panel.genotype_encoding_accounting.tsv
    ├── cohort.gs_panel.record_accounting.tsv
    └── cohort.gs_panel.manifest.json
```

詳細なGS contractは[`docs/gs_panel_data_contract.md`](docs/gs_panel_data_contract.md)、MultiQC contractは[`docs/multiqc.md`](docs/multiqc.md)を参照してください。

### GS dosage contract

variant rows × sample columnsのTSVで、diploid biallelic genotypeを以下へ変換します。

| GT | dosage |
| --- | ---: |
| `0/0` | `-1` |
| `0/1`, `1/0` | `0` |
| `1/1` | `+1` |
| missing / unsupported | `nan` |

phased genotypeはallele countが同じなら同じdosageです。variant keyは`CHROM:POS:REF:ALT`相当の正規化済みidentityを使用し、sample順はVCF header順を維持します。

## MultiQC

Issue #38 / PR #40でMultiQC 1.35を統合しました。containerはtagとdigestで固定し、version check、AI機能、MegaQC uploadを無効化しています。

明示的に取り込む7 input categoriesは以下です。

- raw FastQC ZIP
- trimmed FastQC ZIP
- fastp JSON
- GATK MarkDuplicates metrics
- SAMtools flagstat
- SAMtools stats
- SAMtools idxstats

FastQC / fastp HTML、BWA logs、FASTQ/BAM/VCF、variant QC、GS artifactsはMultiQC inputへ渡しません。synthetic main fixtureでは23 source filesを6 standard module/section source classesとして取り込みます。single-sample / single-read-group境界もreal workflow lineageを通るnf-testで検証しています。

## Tests and CI

```bash
python3 -m unittest discover -s tests/bin -v

NXF_VER=26.04.6 nextflow lint .

NXF_VER=26.04.6 \
  ./nf-test test \
  tests/pipeline/adzuki_snp_pipeline.nf.test \
  tests/modules/*.nf.test \
  --profile "test,docker"
```

現在のtest suiteは、synthetic 2 sample / 3 read-group E2E、diploid/haploid genotype contract、mapping + sort、duplicate marking、reference FAI/dict整合、single-input GenomicsDBImport/GatherVcfs、MIXED variant contract、GS normalization、MultiQC input/output semanticsなどを自動検証します。

GitHub Actionsはmainへのpush / pull requestでNextflow lint、Python unit tests、nf-testを実行します。real WGS cohortはCIでは実行しません。

## Releases

Versioned research releasesは[GitHub Releases](https://github.com/hoso-jpn/adzuki-snp-pipeline/releases)で公開します。`v0.2.0`のrelease scopeとvalidation notesのreview可能な正本は[`docs/releases/v0.2.0.md`](docs/releases/v0.2.0.md)です。

- [`v0.2.0`](https://github.com/hoso-jpn/adzuki-snp-pipeline/releases/tag/v0.2.0) — 最初の明示的なversioned research release。release対象commitは`dc87eb8fdba0294482fc5fdba991c2a34f1569e2`です。

`v0.2.0`が固定するcontractは、synthetic profileとCIで検証されるNextflow DSL2 workflow、FastQCからGS panel生成までのmodule contract、pin済みcontainer digest、および公開WGS 20検体までのreal E2E evidenceです。reference FAI/dictionaryのname、length、order整合性は、sample依存のmappingと全variant calling processに対するhard dataflow gateです。

以下は`v0.2.0`時点で未検証です。

- MultiQC統合はsynthetic fixtureのみで検証しています。real 20検体でのMultiQC再実行は未実施です。
- hard-filter thresholdの生物学的妥当性は未確立です。truth setがないため`FILTER=PASS`は「設定された規則を通過した」以上の意味を持たず、accuracy / precision / recallは主張しません。
- 30検体stageはコストと情報量の比較により明示的にskipしており、失敗ではなく未測定です。
- 50検体超のGenomicsDB batchingと327検体full cohortは未検証です。`genomicsdb_batch_size=50`は初期運用値であり最適化値ではありません。
- BQSRは検証済みknown-sites未確立のため意図的に除外しています。
- production SLAとvariant call精度の優位性は保証しません。

## Real-cohort evidence

### 5 samples — Issue #26

公開WGS 5検体をSeedcore-01でFASTQからGS panelまでE2E実行しました。wall time、resource、storage、variant accounting、input/reference/container checksumを記録しています。

- [`docs/real_cohort_e2e.md`](docs/real_cohort_e2e.md)
- [`docs/real_cohort_e2e_run_manifest.json`](docs/real_cohort_e2e_run_manifest.json)

### 10 → 20 samples — Issue #33

公開WGSを10検体、20検体へ段階拡張しました。20検体でHaplotypeCaller実concurrency 8/8を確認し、GS panelは9,252,873 variants × 20 samplesまで生成しました。

30検体は「失敗したから未実施」ではなく、10/20検体でdownstream memory budget再調整が繰り返し必要になったこと自体をscale findingとし、追加情報量と計算コストを比較して明示的にskipしました。20–30 sample operational decisionは**Conditional Go**です。

- [`docs/real_cohort_scale_validation.md`](docs/real_cohort_scale_validation.md)
- `docs/real_cohort_scale_validation_10sample_manifest.json`
- `docs/real_cohort_scale_validation_20sample_manifest.json`

### Streaming classifier — Issue #35

`CLASSIFY_NORMALIZED_VARIANTS`の全件materializationをlocus-local streamingへ変更しました。formal 10/20-sample targeted replayではclassifier Python peak RSSが約20.8 / 20.9 MiB、host swap delta 0となり、旧実装の28.47 / 53.36 GiBから構造的に改善しました。

- [`docs/classify_normalized_variants_streaming_benchmark.md`](docs/classify_normalized_variants_streaming_benchmark.md)

## 科学的な境界

### BQSR

実施しません。Longxiaodou 4 referenceに対して検証済みのknown-sites resourceを確立していないためです。bootstrap BQSRも現在のscope外です。

### Hard filtering

`FILTER=PASS`は「設定されたhard-filter規則を通過した」ことだけを意味します。truth setがないため、variant callの正しさや最適なthresholdを意味しません。annotation欠損はthreshold通過と区別してQCします。

### GS panel

MAF / call-rate filtering、LD pruning、imputation、GS model trainingはこのrepositoryのscope外です。現在のGS panelはvariant-calling outputから下流解析へ渡すための監査可能なdata packageです。

## 次の開発課題

現在までのIssueはすべて完了しました。リポジトリ全体の再監査後、次のfollow-upを起票しています。

- [#42](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/42) — **P0**: run-level provenance manifestをNextflow DAGへ自動統合する
- [#44](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/44) — **P0**: `BUILD_GS_PANEL`をbounded-memory化する
- [#43](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/43) — **P1**: `SUMMARIZE_FILTER_QC`をstreaming化する
- [#45](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/45) — **P1**: 50検体超のJoint Genotyping scale strategyをtargeted検証する
- [#46](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/46) — **P1**: real cohortでhard-filter annotation分布とthreshold sensitivityを評価する

327検体full E2E、Parabricks/GPU化、全module網羅のnf-test、下流`genomic-prediction-resnet-hybrid`のアズキingestionは、上記より後または別repositoryで扱います。

## データ取扱い

この公開repositoryには公開データ由来のsanitized evidence、synthetic/re-distributable fixtures、code/config/docsのみを置きます。顧客データ、前職由来の非公開データ、生FASTQ/BAM/gVCF/VCF、credential、private infrastructure情報はcommitしません。

実際の受託案件におけるtransfer / retention / deletion / external AI service non-disclosureの前提は[`docs/customer_data_handling.md`](docs/customer_data_handling.md)に記載しています。

## 関連リポジトリ

- [`adzuki-gwas-analysis`](https://github.com/hoso-jpn/adzuki-gwas-analysis) — 公開GWAS summary statisticsの解析
- [`genomic-prediction-resnet-hybrid`](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid) — GBLUP / neural genomic prediction benchmark

これらは現時点で自動化されたSNP→GWAS→GS end-to-end systemとして接続されていません。

## Author

**Yusuke Hosokawa** — independent researcher / AI engineer  
Plant Genetics × Edge AI × Physical AI

- [GitHub](https://github.com/hoso-jpn)
- [researchmap](https://researchmap.jp/hosokawa-yusuke)
- [Florigen AI](https://florigen.ai)
- [Breeding Science (2025)](https://doi.org/10.1270/jsbbs.24058)

## License

MIT License. See [LICENSE](LICENSE).
