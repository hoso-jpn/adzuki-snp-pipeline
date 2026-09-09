# Hard-filter threshold calibration protocol

Issue #46のdecision用プロトコルです。現時点では解析コードと事前定義scenarioを固定し、
20-sample public real cohortの結果値は記録していません。対象VCFがない環境で数値を推定・
転記せず、同一artifact上で実行したmachine-readable evidenceから判断します。

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

## Decision gate

次をすべて確認した後にだけdefault policyを判断します。

1. SNP/indel両方でinput SHA、record count、pipeline/reference/GATK identityが揃う。
2. distributionの各annotationで`present + missing = total`、histogram合計がpresentに一致する。
3. `current`のannotation別および`ANY_FILTER`の`predicted_minus_observed = 0`。
4. lenient/current/stringentのhit率をannotation別・unionで比較し、QUALを過剰解釈しない。
5. truth set不在と20-sample固有の限界をdecision recordへ残す。

判断は次の三択です。

- **KEEP**: evidenceは得たが、truth setなしにdefaultを変える根拠がない。
- **CHANGE IN SEPARATE SCIENTIFIC PR**: 明示的な研究仮説とbefore/after accountingを伴って変更する。
- **DEFER TO DOWNSTREAM DECISION PACK**: 一律defaultにせず、目的別selection policyとして扱う。

現時点の状態は**PENDING REAL-COHORT EVIDENCE**です。default thresholdは変更していません。
