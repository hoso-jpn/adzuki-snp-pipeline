# MultiQC統合レポート

## 目的と出力

MultiQC統合は、親Issue #1からP1 follow-upのIssue #38として切り出したQC review usability機能です。
variant calling、hard filtering、variant accounting、GS panelのscientific semanticsは変更しません。

実行ごとに`${outdir}/qc/multiqc/`へ次をpublishします。

```text
multiqc_report.html
multiqc_data/
multiqc_config.yaml
multiqc_version.txt
```

`multiqc_report.html`をブラウザで開くと統合レポートを閲覧できます。machine-readableな主契約は
`multiqc_data/multiqc_data.json`です。実行に使ったtracked configはbyte-for-byteそのまま
`multiqc_config.yaml`として保存し、runtime versionは`multiqc_version.txt`に保存します。

## Versionとoffline設定

- MultiQC: 1.35
- BioContainers tag: `1.35--pyhdfd78af_1`
- immutable image: `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1@sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac`

`conf/multiqc_config.yaml`は`no_version_check: true`、`no_ai: true`、
`megaqc_upload: false`を明示します。AI summaryを生成せず、MegaQCその他のexternal serviceへ
reportをuploadしません。MultiQCはcontainer内の固定path `/multiqc/input/`に置いたQC reportの
symlinkだけをscanします。これによりNextflow host workdirの絶対pathをpublished provenanceへ
記録しません。

## 明示的input allowlist

| Pipeline artifact | MultiQC 1.35 standard module / section | Granularity | synthetic count |
| --- | --- | --- | ---: |
| raw FastQC ZIP | FastQC (`fastqc`) | read group × mate × raw | 6 |
| trimmed FastQC ZIP | FastQC (`fastqc`) | read group × mate × trimmed | 6 |
| fastp JSON | fastp (`fastp`) | read group | 3 |
| GATK MarkDuplicates metrics | Picard / DuplicationMetrics (`picard`) | biological sample | 2 |
| SAMtools flagstat | Samtools / flagstat (`samtools`) | biological sample | 2 |
| SAMtools stats | Samtools / stats (`samtools`) | biological sample | 2 |
| SAMtools idxstats | Samtools / idxstats (`samtools`) | biological sample | 2 |

synthetic fixtureの合計input provenanceは23 filesです。標準moduleだけを使用し、custom parserや
pluginはありません。各categoryは必須かつnon-emptyで、`require_logs: true`により明示した
moduleがparseできない場合もfailureになります。

FastQC HTMLとfastp HTMLは、同じrunを二重取り込みしないため除外します。BWA mapping logs、
BCFtools variant stats、variant/filter QC TSV、VCF、BAM、FASTQ、GS normalized VCF、GS accounting、
GS panelもIssue #38のscope外であり、scanしません。`${outdir}`全体をrecursive scanすることは
ありません。

## Sample naming

FastQC inputはcopyせず、task workdir内のsymlinkを次の決定的な名前にします。

```text
<read_group_id>.raw.R1.fastq.gz
<read_group_id>.raw.R2.fastq.gz
<read_group_id>.trimmed.R1.fastq.gz
<read_group_id>.trimmed.R2.fastq.gz
```

MultiQC sample IDは末尾のFASTQ/FastQC拡張子をcleanした
`<read_group_id>.<raw|trimmed>.<R1|R2>`です。元FASTQが異なるdirectoryで同じbasenameでも衝突
せず、raw/trimmed、R1/R2、read groupを常に区別します。MultiQC既定cleanerの`.trim` tokenは
`trimmed.R1/R2`を縮約するため、tracked configの明示的な`fn_clean_exts`から除外しています。

fastpは`<read_group_id>`のread-group granularityです。MarkDuplicatesと全SAMtools reportは
`<sample_id>`のbiological-sample granularityです。この粒度差は意図的です。

## Machine-readable test contractと決定性

nf-testは`multiqc_data/multiqc_data.json`の公開された安定subsetをparseし、次を検証します。

- `report_saved_raw_data`のFastQC 12、fastp 3、Picard duplication 2、SAMtools各section 2 samples
- `report_data_sources`のstandard module/section、7 source class、合計23 source files
- raw/trimmed、R1/R2、全read-group sample IDの集合
- sourceが`/multiqc/input/`配下のallowlisted QC reportだけであること
- project/output/work/tmpのabsolute pathがJSON、source provenance、logへ漏れないこと

byte-stable contractはtracked config、source naming rule、module allowlist、上記semantic subsetです。
HTML、生成日時、runtime metadata、MultiQC log、plot representationはbyte-stableではなく、golden
hashを保証しません。

## 制約

synthetic 3-read-group／2-sample fixtureでのみ検証しています。実データ20検体のMultiQC rerunは
実施しておらず、327検体でも未検証です。MultiQC resource scalingのbenchmarkもIssue #38の
scope外です。variant QCやGS artifactは統合対象ではありません。
