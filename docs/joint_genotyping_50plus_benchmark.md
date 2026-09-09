# 50+ sample Joint Genotyping targeted benchmark protocol

Issue #45で、327-sample full FASTQ→GS E2Eへ進む前にJoint Genotypingだけを隔離して評価する
ためのプロトコルです。現時点では入力検証・candidate plan・evidence schemaだけを固定し、
50+ sample実測値や採否判断は記録していません。

## 前提

- 公開データ由来の51検体以上のgVCFを使い、`genomicsdb_batch_size=50`の実batchingを発生させる。
- 全gVCFは同じreference、GATK、pipeline SHA、ploidy、HaplotypeCaller条件で生成されたものに限る。
- 全327 FASTQの再処理は前提にせず、既存gVCFを再利用する。
- 異なるlineageのgVCF混在、customer/private data、327-sample full E2Eは対象外。

## Preparation gate

manifestは`sample_id,gvcf,gvcf_index`の3列CSVです。次のCLIは51検体未満、sample名不一致、
重複path、index欠損、reference contig name/length/order不一致を拒否し、全入力checksum、
GATK sample-name-map、baseline/candidate interval plan、未記入のevidence JSONを生成します。

```bash
python3 bin/prepare_joint_genotyping_benchmark.py \
  --gvcf-manifest gvcfs.csv \
  --reference-fai GCF_016808095.1_ASM1680809v1_genomic.fna.fai \
  --pipeline-commit <full-40-hex-sha> \
  --gatk-container 'broadinstitute/gatk:4.6.2.0@sha256:...' \
  --batch-size 50 \
  --window-size-bp 20000000 \
  --small-scaffold-max-bp 1000000 \
  --sample-name-map-output sample_name_map.tsv \
  --interval-plan-output interval_plan.tsv \
  --evidence-template-output benchmark_evidence.json
```

`sample_name_map.tsv`はGATK 4.6.2.0の契約どおりheaderなしの3列
（sample name、gVCF path、明示的index path）です。形式の正本は
[pinned GATK source](https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/genomicsdb/GenomicsDBImport.java#L126-L145)
を参照してください。このfileはbenchmark実行用で絶対pathを含むためcommitしません。公開repositoryへ
commitするのは、basename/checksumへsanitizeされたevidenceとmethod/result documentだけです。

## Targeted experiment order

一度に複数要因を変えず、evidence templateの順に比較します。

| ID | gVCF input | interval | Reblock | consolidate | 判定対象 |
| --- | --- | --- | --- | --- | --- |
| E0 | repeated `--variant` | per-contig baseline | No | No | 現行baseline |
| E1 | `--sample-name-map` | per-contig baseline | No | No | sample-name-map単独効果 |
| E2 | sample-name-map | chromosome split + small-scaffold group | No | No | interval戦略単独効果 |
| E3 | sample-name-map | candidate | Yes | No | Reblock単独効果 |
| E4 | sample-name-map | candidate | No | Yes | consolidate単独効果 |

E3のReblockGVCFsは、入力gVCFのchecksumを保持したまま別の派生input setとして作り、size、時間、
sample/header/order、Joint Genotyping後のvariant/accounting equivalenceを比較します。採用前に
scientific equivalenceを確認し、単なる容量削減だけでGOにしません。

## 必須計測

各experimentで次をevidence JSONへ記録します。

- GenomicsDBImport wall time / peak RSS / workspace bytes / file count / retry / swap delta
- GenotypeGVCFsとGatherVcfsのwall time
- output sample countとsample order checksum
- variant/accounting checksum
- batch size、interval plan、Reblock/consolidate状態
- reference identity、pipeline SHA、GATK container、全input checksum

## Decision gate

sample-name-map、interval戦略、Reblock、consolidateはそれぞれ`ADOPT`/`REJECT`/`DEFER`を根拠付きで
記録します。327-sample full runは、resource envelopeを実測から外挿し、storage/timeにheadroomが
あり、retry/swapが許容範囲で、sample/order/variant accountingに回帰がない場合のみ`GO`です。
不確実性が残る場合は`CONDITIONAL GO`、lineageまたはresource contractを満たせない場合は`NO-GO`
とします。

現在の状態は**PENDING REAL BENCHMARK**です。production workflow、resource label、batch-size default、
interval strategy、Reblock/consolidate policyは変更していません。
