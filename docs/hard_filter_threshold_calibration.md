# Hard-filter threshold calibration protocol

Issue #46のプロトコルと2026-09-10の20-sample public cohort実測です。
Issue #33のfiltered SNP/indel VCFを再利用し、事前定義scenarioを変更せずに解析しました。

## 目的と解釈境界

SNP/indel別にQD、QUAL、SOR、FS、MQ、MQRankSum、ReadPosRankSumの分布、欠損率、
現行および少数の代替thresholdによるvariant-setの変化を定量化します。truth setがないため、
accuracy、precision、recall、生物学的妥当性、「最適threshold」は評価しません。とくにQUALは
sample/cohort sizeへの依存が強く、20検体で得た値をcohort非依存の品質基準と解釈しません。

## 事前定義scenario

正本は[`../conf/hard_filter_sensitivity_scenarios.json`](../conf/hard_filter_sensitivity_scenarios.json)
です。

- `current`: 現在の`nextflow.config` default。既存FILTER tagとの一致確認にも使う。
- `lenient`: currentより除外が少なくなる方向の探索条件。
- `stringent`: currentより除外が多くなる方向の探索条件。

代替値は感度幅を見るためのexploratory conditionであり、候補値自体を推奨値とはみなしません。
実データを見た後に都合よくscenarioを追加・変更する場合は、config SHAと変更理由を別記します。

## 実行

20-sample run manifestが指す同一pipeline/reference lineageのfiltered SNP/indel VCFを使います。
異なるreference、GATK version、pipeline SHAのartifactを混在させません。

```bash
python3 bin/analyze_hard_filter_sensitivity.py \
  --filtered-vcf cohort.snp.filtered.vcf.gz \
  --cohort-id PRJNA1138464_20sample \
  --variant-type snp \
  --scenario-config conf/hard_filter_sensitivity_scenarios.json \
  --distribution-output snp.annotation_distribution.tsv \
  --sensitivity-output snp.threshold_sensitivity.tsv \
  --summary-output snp.threshold_sensitivity.summary.txt
```

indelも`--variant-type indel`で同様に実行します。summaryにはinputとscenario configのSHA256を
記録します。TSVの`current`行では、annotation値から再計算したhit数とVCFに記録済みのFILTER
tag数の差を`predicted_minus_observed`として出します。0でない場合はthreshold判断を中止し、
input/config lineage、GATK expression semantics、annotation parsingを調査します。

件数差が0でも、異なるrecordの過剰・不足が相殺される可能性があります。
`discordant_records`はcurrentの予測と観測をrecordごとに照合した不一致件数です。
annotation別と`ANY_FILTER`の両方で0を要求します。SNP/indelの想定外tagや未評価のFILTER `.`は
入力エラーとし、出力は全tableの書き込み成功後に公開します。`ANY_FILTER`のpresentは
対象annotationが一つ以上あるrecord数であり、全annotationが揃うrecord数ではありません。

## Decision gate

次をすべて確認した後にだけdefault policyを判断します。

1. SNP/indel両方でinput SHA、record count、pipeline/reference/GATK identityが揃う。
2. distributionの各annotationで`present + missing = total`、histogram合計がpresentに一致する。
3. `current`のannotation別および`ANY_FILTER`の`predicted_minus_observed = 0`かつ`discordant_records = 0`。
4. lenient/current/stringentのhit率をannotation別・unionで比較し、QUALを過剰解釈しない。
5. truth set不在と20-sample固有の限界をdecision recordへ残す。

判断は次の三択です。

- **KEEP**: evidenceは得たが、truth setなしにdefaultを変える根拠がない。
- **CHANGE IN SEPARATE SCIENTIFIC PR**: 明示的な研究仮説とbefore/after accountingを伴って変更する。
- **DEFER TO DOWNSTREAM DECISION PACK**: 一律defaultにせず、目的別selection policyとして扱う。

## 20-sample実測と判断: KEEP

公開BioProject **PRJNA1138464**、reference **GCF_016808095.1 / Longxiaodou 4**、
source pipeline `556f38fd93008e0e4093a2fbd836e8d141afd8a4`、Nextflow 26.04.6、GATK 4.6.2.0。
source manifestはrepository内の20-sample manifestとbyte-identicalです。
実行前にreference FASTA/FAI/dictのSHA256、VCFのordered contigs、20検体のsample order、
record数、raw VCFのmanifest checksumを照合しました。input SHA256と全accessionは
[`measurement.json`](evidence/issue46/measurement.json)に記録しています。

| Variant type | Records | lenient ANY_FILTER | current ANY_FILTER | stringent ANY_FILTER |
| --- | ---: | ---: | ---: | ---: |
| SNP | 10,296,980 | 300,698 (2.9203%) | 1,278,951 (12.4206%) | 3,543,704 (34.4150%) |
| indel | 1,649,635 | 1,388 (0.0841%) | 6,802 (0.4123%) | 97,608 (5.9169%) |

分母は各variant typeの全recordです。各scenarioは複数thresholdを同時に変えるため、union差を
単一annotationの効果とは解釈しません。per-annotationのpresent分母によるhit率は別に出力します。

currentのSNP tag件数はQD 67,992、QUAL 0、SOR 410,137、FS 10,371、MQ 890,542、
MQRankSum 0、ReadPosRankSum 6です。indelはQD 6,800、QUAL 0、FS 5、ReadPosRankSum 0です。
tagは重複し得るため単純和はunionと一致しません。全current行で
`predicted_minus_observed = 0`かつ`discordant_records = 0`でした。

QD欠損はSNP 18、indel 12 recordsです。MQRankSumのevaluable率はSNP 40.4789%、indel 30.9910%、
ReadPosRankSumはSNP 40.3749%、indel 30.5634%です。欠損recordをthreshold通過と扱いません。
QUALの最小値は両typeとも30で、このartifact上の`QUAL < 30`は0件です。これは上流selectionを経た
当該cohortの観測であり、QUAL閾値の一般的有効性やsample数非依存性を示しません。

全7 annotationのmin/max/mean・固定bin分布・欠損率は
[SNP distribution](evidence/issue46/snp.annotation_distribution.tsv) /
[indel distribution](evidence/issue46/indel.annotation_distribution.tsv)、scenarioごとの件数と率は
[SNP sensitivity](evidence/issue46/snp.threshold_sensitivity.tsv) /
[indel sensitivity](evidence/issue46/indel.threshold_sensitivity.tsv)が正本です。
各histogram合計がpresent、present + missingがtotalに一致し、lenient ≤ current ≤ stringentを確認しました。

解析はseedcore-01 / Python 3.12.3でSNP 99.19秒・22.25 MiB RSS、indel 15.08秒・21.97 MiB RSS。
host swap-used増分とswap-out増分は両方0でした。SNP中のhost swap-inは12 KiBであり、
host全体のカウンタのためこのprocessへの帰属はできません。

**KEEP**: 現行filterは保存済みannotationから再現でき、代替条件の影響も定量化できました。
truth setがなく、variant集合の増減から誤陽性・誤陰性は判定できないため、defaultを変える根拠は
ありません。生物学的妥当性・accuracy・precision/recall・最適性は未検証のままです。
判断はこの20検体cohortに限定し、327検体への外挿や下流GS/GWAS性能の改善を主張しません。
