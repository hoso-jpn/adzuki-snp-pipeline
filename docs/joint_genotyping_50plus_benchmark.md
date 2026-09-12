# 50+ sample Joint Genotyping targeted benchmark protocol

Issue #45で、327-sample full FASTQ→GS E2Eへ進む前にJoint Genotypingだけを隔離して評価する
ためのプロトコルです。入力検証・candidate plan・実行／監査helper・evidence schemaを用意しました。
現時点では50+ sample実測値や最終採否判断は記録していません。

## 前提

- 公開データ由来の51検体以上のgVCFを使い、`genomicsdb_batch_size=50`の実batchingを発生させる。
- 全gVCFは同じreference、GATK、pipeline SHA、ploidy、HaplotypeCaller条件で生成されたものに限る。
- 全327 FASTQの再処理は前提にしない。再利用可能なpublic artifactを優先し、同一lineageの51 gVCFを確保する。
- 異なるlineageのgVCF混在、customer/private data、327-sample full E2Eは対象外。

## Preparation gate

manifestは`sample_id,gvcf,gvcf_index`の3列CSVです。次のCLIは51検体未満、sample名不一致、
重複path、index欠損、reference contig name/length/order不一致を拒否し、全入力checksum、
GATK sample-name-map、baseline/candidate interval plan、未記入のevidence JSONを生成します。
相対pathはmanifestの親directoryを基準に解決します。sample名の空白は除去せず拒否します。
出力先がmanifest・FAI・gVCF・indexと重なる場合は拒否します。3成果物はすべて書き込み後に公開し、
通常のwrite/rename failureでは既存結果を保ちます。process/host crashをまたぐtransactionは保証しません。

このCLIは**実行許可gateを完了しません**。`pipeline_commit`と`gatk_container`はCLIからの宣言であり、
FAIの一致はreference配列の同一性を証明しません。`.tbi`の名前一致もindex内容の対応を証明しません。
evidenceの`input_validation.benchmark_ready` / `lineage_verified`は常に`false`です。
実行前にsource-run manifestと各gVCFのchecksumを結び、FASTA identity、GATK version、
ploidy・HaplotypeCaller条件、index内容の対応を別途検証する必要があります。

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

実行前の全件lineage/index検証、測定方法、資源制限、service復旧、失敗保存は
[`benchmarks/issue45/README.md`](../benchmarks/issue45/README.md)に定義します。
実行helperはproduction checkoutと分離し、production SHAとhelper SHAを別々に記録します。

## Targeted experiment order

一度に複数要因を変えず、evidence templateの順に比較します。

| ID | gVCF input | interval | Reblock | consolidate | 判定対象 |
| --- | --- | --- | --- | --- | --- |
| E0 | repeated `--variant` | per-contig baseline | No | No | 現行baseline |
| E1 | `--sample-name-map` | per-contig baseline | No | No | sample-name-map単独効果 |
| E2a | sample-name-map | chromosome split、scaffold個別 | No | No | E1に対するwindow分割単独効果 |
| E2b | sample-name-map | chromosome split + small-scaffold group | No | No | E2aに対するgrouping単独効果 |
| E3 | sample-name-map | candidate | Yes | No | Reblock単独効果 |
| E4 | sample-name-map | candidate | No | Yes | consolidate単独効果 |

E3のReblockGVCFは、入力gVCFのchecksumを保持したまま別の派生input setとして作り、size、時間、
sample/header/order、Joint Genotyping後のvariant/accounting equivalenceを比較します。採用前に
scientific equivalenceを確認し、単なる容量削減だけでGOにしません。
E3/E4はいずれもE2bと比較します。evidenceの`compare_to`が比較対象を明示します。
groupingはreference dictionary中で連続するsmall scaffoldだけに限定し、GatherVcfsの入力順を
保ちます。windowは1-based closed intervalで、gap/overlapのないpartitionとして検証します。
small scaffold groupの合計もwindow size以下に制限し、大量のscaffoldが一つの巨大taskに
集約されることを防ぎます。

実測前に、各interval taskは8 CPU / 16 GiB、最大3 task並列、windowは20 Mb、small scaffoldは
1 Mb以下と固定します。これは実測前の事前登録値です。16 GiBでE0が実際にOOM-killedとなったため、
GenotypeGVCFsのみ実測に基づく共通ceilingへ改めました。経緯と根拠は
「2026-09-11: E0 baselineは16 GiBで失敗」以降の各節を参照してください。
16 GiBでのbaseline failureは置き換えず保持し、比較用measurementとは別に記録します。GenotypeGVCFsには全experimentで
`--only-output-calls-starting-in-intervals true`を指定し、window境界のvariant開始位置の所有を
一意にします。E3はGQ bands 20/100、`keep-all-alts=true`、`floor-blocks=false`、
`drop-low-quals=false`とし、追加のallele droppingを導入しません。

full-record checksumを保持したうえで、sample/contig order、GT、AC/AN、各INFO/FORMAT差を監査します。
pinned GATKの高QD補正には乱数が使われるため、task分割によるQD差を別に数えます。
QUAL/AD/GTから補正対象であることを説明でき、他の全fieldが同一で、current QD filter判定も
変わらない場合だけ「説明済み」と分類します。欠損・NaN・評価不能なannotationを同等扱いしません。
この分類は自動採用判断ではなく、実測後のscientific reviewの材料です。

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

production workflow、resource label、batch-size default、interval strategy、
Reblock/consolidate policyは変更していません。実測の進捗は以下の日付別節に追記します。

## 2026-09-10: 51-sample preparation in progress

Ownerによる追加public FASTQ取得の承認後、ENAの全684 runを再取得しました。
WGS / paired-endは327 run、327 unique BioSampleで、RAD-Seq / single-end 357 runを除外しました。
既存20検体と、残りのeligible runをaccession昇順で選んだ31検体を固定しました。
全51検体はPRJNA1138464、Vigna angularis、Illumina HiSeq Xです。
この選定はvariant結果に依存しませんが、集団の無作為標本でもありません。

- [全候補と除外理由](evidence/issue45/candidate_cohort.tsv)
- [固定51検体と公開checksum](evidence/issue45/cohort_selection.json)
- [既存20検体のchecksum / lineage audit](evidence/issue45/old20_lineage_audit.json)
- Public metadata source: [ENA filereport API](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1138464&result=read_run&format=tsv)

旧20検体の40 FASTQはENA published MD5とsource manifest SHA256に一致し、20 gVCFのSHA256も
source manifestと一致しました。reference FASTA/FAI/dictのSHA256と、全20 HaplotypeCallerの
canonical parametersも一致します。旧実行から現在のmainへの上流module差分には、scientific
commandの変更はありません。ただし旧20 gVCFには前段10検体runからのcached outputが含まれ、
生成元SHAは単一ではありません。今回のprotocolのstrict single-SHA contractを維持するため、
**FASTQは再利用し、gVCFは全51検体REGENERATE**とします。旧gVCFを新SHA由来とは表記しません。

固定production SHAは`3158ca50c2c13c31bdc80db302c7df4bbb5670bf`です。
production checkoutとPR #61 checkoutを分離し、productionのmoduleを変更せず、別途hashを
記録するgVCF専用wrapperから呼び出します。synthetic 2検体（うち1検体は2 read groups）で、
wrapperと固定mainのfull workflowのgVCF全recordがbyte-identicalであることを確認しました。
これは実行経路の回帰検証であり、51-sample real benchmarkの代替ではありません。

不足31検体のFASTQ（約150.6 GB）の取得と、再利用分を含む全102 filesのchecksum検証が完了しました。
全51検体の公開FASTQ容量は207.4 GBです。実gVCFの生成・全件validation・E0〜E4測定はまだ完了しておらず、
`benchmark_ready=false` / `lineage_verified=false`、PRはDraftを維持します。
327-sampleの最終Gateは実測後に決定します。

全102 FASTQのpublished MD5検証・SHA256固定後、最初のgenerationを開始しましたが、
Docker内から外部`bin` symlinkの参照先が見えず、reference validationがexit 127で停止しました。
gVCFは0件で、4 FASTP taskのみ完了しています。解析container停止とlocal trial復旧を確認し、
失敗runの約47.1 GBとログを保持しました。
[sanitized failure evidence](evidence/issue45/generation_failure_20260910_194131.json)に記録します。
production codeは変更せず、git-tracked scriptsをchecksum付きでlaunch directoryへコピーする
staging修正と、mount範囲を限定した回帰を追加しました。新run IDで全51検体を再実行し、
失敗runのpartial outputは入力に使いません。

## Earlier inventory, before download authorization

`seedcore-01`の既知の公開cohort保管領域（Issues #11/#26/#33/#35）、pipeline checkouts、
Issue #44 replay、data volumeのdirectory一覧を調査しました。既存run metadataから退避済みの
Issue #33 artifactを発見しましたが、公開gVCFのheaderで確認できたユニーク検体は20です。
5/10/20検体runの重複は追加sampleとして数えていません。詳細は
[`evidence/issue45_public_gvcf_inventory.json`](evidence/issue45_public_gvcf_inventory.json)に記録します。

51検体以上のhomogeneous inputがないためE0〜E4の実測と327検体resource envelopeは未完了です。
full 327-sample runの実行判断は**NO-GO（benchmark evidence不足）**とし、各手法の採否は保留します。
新規FASTQ download、sample名の書き換え、異なるlineageの混在は行っていません。
再開には、最低31の追加ユニークpublic sampleを含む51+ homogeneous gVCFと、全inputの
reference・pipeline/GATK・index provenanceが必要です。現在の実機はMemAvailable約3 GiBであり、
targeted GenomicsDB実行には併せて十分なRAM余裕を確保する必要があります。

コード上もsource-run provenance検証とindex内容検証を完了するまで実行gateを開きません。
PR #61はDraft、Issue #45はopenのまま維持します。

## 2026-09-11: 51検体gVCF生成完了とlineage検証

`issue45-51gvcf-20260910-200521`で全51検体のgVCFを生成し、独立validationを完了しました。
51 unique run accession、51 unique BioSample、51 unique gVCF SHA256で、gVCFの複製や
sample名の書き換えはありません。ordered contigs、shared INFO/FORMAT/FILTER/ALT定義、
expanded HaplotypeCaller parametersはいずれも51検体で単一値です。全gVCFのBGZF CRC/EOFと、
通常読込みとindex経由の全record一致を確認しました。入力gVCFは合計48.03 GB、3,566,258,766 records、
全taskがattempt 1で成功しています。この時点で`benchmark_ready=true` / `lineage_verified=true`
となりました。

## 2026-09-11: E0 baselineは16 GiBで失敗（negative evidenceとして保持）

固定した16 GiB / heap 15 GiBで最初のE0を実行し、最長contig `NC_068970.1`（65,407,200 bp）の
GenotypeGVCFsがexit 247 / `OOMKilled=true`で停止しました。peak RSSは15.98 GiBで、自身の
16 GiB上限の99.9%です。progress meterは約1,218,671 records/分から3 records/分へ崩壊しており、
これはGCのdeath spiralで、live setが本当にheap上限へ達したことを示します。
host側のMemAvailable最小値は83.23 GiB、swapは3.318→3.306 GiBでdeltaは負であり、
host資源の枯渇ではなくcontainer cgroup上限単独が原因です。

この結果は**production baselineに関する有効なnegative evidenceとして保持**します。
通るまでallocationを上げて「E0成功」に置き換えることはしていません。
[sanitized failure evidence](evidence/issue45/e0_oom_failure_20260911_160616.json)に
`allocation_regime = production_baseline_16gib`として記録しています。

完了した7 taskから、GenomicsDBImportのpeak RSSは1.27–1.55 GiBでinterval長にほぼ依存せず、
16 GiB割当は十分余裕があると分かりました。GenotypeGVCFsのpeak RSSはinterval長とともに増加します
（51検体で概ね0.28 GiB/Mb）。ただしこれは**51検体・本reference・本GATK version・本ploidy・
本HaplotypeCaller条件での観測**であり、sample数方向のscalingを測ったものではなく、327検体へ
線形外挿しません。さらに最長2 contigの読み値自体が検証対象の上限で切られているため下限値です。

`genomicsdb_batch_size=50`の実batchingは、開始した全7 intervalで
`Importing batch 1 with 50 samples` / `Done importing batch 1/2` /
`Importing batch 2 with 1 samples` / `Done importing batch 2/2` /
`Import of all batches to GenomicsDB completed!`として確認済みです。
51検体を入力しただけでなく、GATK側で2 batchとして実処理されました。

## 2026-09-12: 真のpeak RSS calibrationと共通ceilingの決定

上限で切られた読み値から比較用ceilingを決めないため、本repositoryの既存methodology
（`nextflow.config`のIssue #30/#33: 対象processを単独で十分大きなceilingで再測定する）に従い、
最長contigのGenotypeGVCFsを単独・96 GiB / heap 80 GiBで再測定しました。保持した失敗evidenceを
変更しないよう、GenomicsDB workspaceはコピーして使用しています。

| 指標 | 16 GiB試行（上限で切断） | calibration（真値） |
| --- | --- | --- |
| 結果 | exit 247 / OOMKilled | exit 0 / 完了 |
| peak RSS | 15.98 GiB（上限の99.9%） | **17.91 GiB** |
| wall time | 84.2分で88%地点 | 49.0分で完走 |
| major page faults | 10,565,494 | 80 |
| host swap | — | 0.00 GiB |

真の必要量17.91 GiBは16 GiB上限を超えており、OOMは測定上の人工物ではなく実際の資源不足です。
major page faultsの10,565,494対80と、84.2分で88%対49.0分で完走という2つの独立した読み値が、
progress meterに現れたGC death spiralを裏付けます。

したがって共通ceilingは**GenotypeGVCFs 32 GiB / 派生heap 26 GiB**（native reserve 6 GiB、
実測真値の1.79倍）とし、**E0〜E4の全experimentへ同一に適用**します。
allocation差が比較対象ペアの交絡要因にならないようにするためです。
Java heapはcontainer上限からnative reserveを引いて導出します。GenomicsDBのTileDB bufferと
JVMのmetaspace/stack/code cacheはheap外に存在するため、heapを独立した数値として書くと
上限と一致してしまう今回のbug classが再発します。GenomicsDBImportは実測1.55 GiBに対して
16 GiBのままで変更しません。

[calibration evidence](evidence/issue45/genotype_memory_calibration_20260912.json)に記録しています。
このcalibrationはE0〜E4の比較experimentではなく、16 GiB baseline failureとも比較用measurementとも
別に報告します。公開evidenceの`allocation_regime`は、比較ペアを交絡し得る値の一様性を明示的に検査します。

科学的parameter、threshold、ploidy、reference、interval semanticsはいずれも変更していません。
E0/E1のper-contig planも同一のままです。変更したのは資源割当のみです。
