# 50+ sample Joint Genotyping targeted benchmark protocol

Issue #45で、327-sample full FASTQ→GS E2Eへ進む前にJoint Genotypingだけを隔離して評価する
ためのプロトコルと、その実測結果・判定の記録です。51検体でのE0〜E4実測、各手法の採否、
sample数scaling補助測定、327検体resource envelopeとGateは末尾の日付別節に記録します。

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

## 2026-09-13: E0 completed under the common ceiling

共通ceiling（32 GiB / 派生heap 26 GiB）でE0を完走しました。

| 指標 | 実測値 |
| --- | --- |
| status | COMPLETED |
| interval数 | 36（chromosome 11 + scaffold 25） |
| retry / OOM | 0 / 0 |
| GenomicsDBImport peak RSS | 1.03–2.13 GiB（interval長にほぼ非依存） |
| GenotypeGVCFs peak RSS | 最大17.55 GiB |
| GatherVcfs | 3秒 / 0.34 GiB |
| IndexFeatureFile | 36秒 |
| workspace | 55.14 GB / 3,132 files |
| output | 51 samples / 13,219,170 variants / 1,247,857,386 called alleles |
| sample order / contig order / AC-AN accounting | すべてvalid |
| elapsed（validation込み） | 4.17時間（task時間合計11.42時間、3並列） |
| host MemAvailable最小 / swap delta | 88.4 GiB / 0.00 GiB |

GenotypeGVCFsのpeak RSS 17.55 GiBは、単独96 GiBで測ったcalibrationの17.91 GiBと2%以内で一致します。
これによりceilingが読み値を切っていないことと、calibrationが妥当だったことが相互に裏付けられ、
同時に16 GiBでの失敗が実際の資源不足だったことも確認されます。

`genomicsdb_batch_size=50`の実batchingは**全36 intervalで**`((1,50),(2,1))`の単一patternとなり、
全intervalが`2/2`まで完了しました。serial reader fallbackは0件です。

### interval長とpeak RSSの関係（unclipped実測に基づく修正）

chromosome-scale 11 intervalすべてがceilingに触れていない実測値なので、これを用いて回帰しました。

```
peak_rss_GiB = 3.172 + 0.2192 * interval_Mb     (R^2 = 0.9332, max residual 0.88 GiB)
```

関係は比例ではなく**affine**です。約3.2 GiBのinterval長に依存しない固定成分があるため、
peak/長さの単純比はintervalが長くなるほど0.3373→0.2683 GiB/Mbと低下します。
原点を通る単一slopeを使うと短いintervalを過小評価します（20 Mbで約5.6 GiB対本fitの約7.6 GiB）。
16 GiB runで得た約0.28 GiB/Mbは、最長2 contigの読み値が検証対象の上限で切られた比であり、
sizingやprojectionには使用せず、当該失敗runが示した値としてのみ保持します。

20 Mb windowでの予測値は**observedではなくassumption**です。E0には27.7 Mb未満の
chromosome-scale intervalが存在せず、20 Mb地点は内挿にすぎません。これはE2aが同一cohort・
同一ceilingで実測するため、仮定ではなく測定値に置き換わります。

両係数は51検体でのみfitしたものであり、sample数方向のscalingは測定していません。
327検体へ外挿していません。

- [E0以降のsanitized metrics](evidence/issue45/benchmark_results.json)
- [memory modelとscope](evidence/issue45/genotype_memory_model_51_samples.json)

## 2026-09-13: E1 sample-name-map（ADOPT、条件付き）

E0との差分を入力形式のみに限定して実行しました。gVCF 51件、per-contig plan、reference、
batch size、ceiling、並列数はすべて同一です。

### correctness

| 指標 | 結果 |
| --- | --- |
| shared variants | 13,219,170 |
| identical records | 13,219,170 |
| different records | 0 |
| left/right only variants | 0 / 0 |
| record / GT / accounting / header SHA256 | すべて一致 |
| sample order / contig order | 一致 |

出力callsetは全recordで一致します（record / GT / accounting / shared header checksumが一致）。
gathered fileのbyte列は、各runの引数と日時を記録する`##GATKCommandLine` header行の分だけ異なり得るため、
一致の判定はfile bytesではなくrecord単位のchecksumで行います。入力形式は結果を一切変えませんでした。

### performance

| 指標 | E0 | E1 | 差 |
| --- | --- | --- | --- |
| import wall合計 | 5.33 h | 5.34 h | +0.17% |
| genotype wall合計 | 6.09 h | 6.09 h | -0.06% |
| elapsed（validation込み） | 4.17 h | 4.14 h | -0.63% |
| workspace | 55.14 GB / 3,132 files | 同一 | 0 |
| import peak RSS | 2.13 GiB | 1.76 GiB | — |

51検体では**性能差は測定されません**。wall timeの差はrun間ノイズの範囲で、import peak RSSは
いずれも16 GiB割当を大きく下回るため意味のある効果ではありません。

### operational

| 指標 | E0 | E1 |
| --- | --- | --- |
| argv要素数 | 148 | 50 |
| command文字数 | 2,871 | 1,013 |
| `--variant`引数 | 51 | 0 |

**command長は採用理由になりません。** repeated-variant形式を327検体へ投影しても13,359文字で、
本hostの`ARG_MAX` 2,097,152の約0.6%にすぎません。よく挙げられるこの理由は本件では成立しない
ものとして明記します。

### limitation: sample orderは検証できていない

本cohortの51 run accessionは既にlexicographic順であり、repeated-variant順・map file順・
出力sample順のすべてがsorted順と一致します。GATKはmap指定時にsample順をmapから導くため、
**sorted順と入力順が一致する本cohortでは、両形式のordering差を検出できません**。
51検体でrecordが一致したことは、sample名がsorted順でないcohortで両形式のorderingが
一致する証拠にはなりません。327検体で採用する前に、production samplesheet順とsorted順の
一致を確認するか、mapのsorted順を契約として受け入れた上で出力sample順を明示的に検証してください。

### 判定

**ADOPT**（条件は下記のとおり解消済み）。根拠は性能ではなくdata integrityとauditabilityです。
mapはsampleごとに明示的なindex pathを持ち、GATK自身のmap validationを有効化でき、
どのgVCFがrunに入ったかを単一のchecksum可能なartifactとして記録します。
これは51検体より327検体で価値が高くなります。測定可能なコストはなく、callsetも変えません。
command長は根拠に含めません。採用は上記sample ordering確認を条件とします。

性能上の利点は主張しません。truth setがなく、そもそもcallsetのrecordが一致しているため、
精度改善も主張しません。

### 条件の解消（327検体のsample order確認）

`candidate_cohort.tsv`の全eligible setで確認しました。WGS / paired-endは327 run、
327 unique run accession、327 unique BioSampleです。本pipelineはsample IDに
**run accessionをそのまま**使用します。

- 327 run accessionは既にlexicographic順
- 全accessionが同一長（11文字）のためlexicographic順と数値順が一致（zero-padding問題なし）

したがって327検体でもmapのsorted順はsamplesheet順と一致し、`--sample-name-map`採用によって
output sample orderが変わることはありません。**条件は解消**です。

ただしBioSample accessionはrun accession順とは**一致しません**。BioSampleで並べ替えると
run順が変わります。将来sample IDをrun accession以外（BioSample等）へ変更する場合は、
この性質を再確認する必要があります。

- [sample-name-map decision evidence](evidence/issue45/sample_name_map_decision.json)

## 2026-09-13: E2a ADOPT / E2b REJECT（interval strategy確定）

data-contract gateとarchitecture gateを分離して評価しました。
data-contract gateの通過は「比較可能な候補への昇格」であり、採用ではありません。

### E2a: 20 Mb chromosome window splitting → **ADOPT**

data-contract gate（対E1）: **ADOPT_CANDIDATE**

| 分類 | 件数 | 差分中の割合 |
| --- | --- | --- |
| directly explained（raw QD ≥ 35.01） | 3,251,539 | 99.957% |
| reconstruction-limit consistent | 1,397 | 0.043% |
| **unresolved** | **0** | **0.0%** |

QD<2 FILTER membership flip 0件、非QD差分0件。variant set・genotype・AC/AN・QUAL・FILTER・
FORMAT・全sample列・sample order・contig orderはすべて一致し、差分は`INFO/QD`のみです。
1,397件は[35.000000, 35.009756]に収まり、explainedへは昇格させません。
詳細は[E2a QD difference audit](evidence/issue45/e2a_qd_difference_audit.json)。

architecture gate（対E1）:

| 指標 | E1 | E2a | 差 |
| --- | --- | --- | --- |
| **GenotypeGVCFs peak RSS** | 17.56 GiB | **7.43 GiB** | **-57.7%** |
| GenotypeGVCFs wall合計 | 6.09 h | 6.04 h | -0.8% |
| GenomicsDBImport peak RSS | 1.76 GiB | 1.76 GiB | 0.0% |
| workspace bytes | 55.14 GB | 55.12 GB | -0.0% |
| workspace files | 3,132 | 4,611 | +47.2% |
| elapsed | 4.14 h | 4.09 h | -1.4% |

20 Mb window 20本の実測peak RSSは5.40–7.43 GiBでした。
唯一ceilingに迫っていたprocessを57.7%削減し、wall timeもworkspace sizeも悪化させません。
memory要件を「最長contigの性質」から「選んだwindow長の性質」へ移す点が本質です。
代償はworkspace file数+47.2%とtask数+17であり、資源制約ではなくscheduling overheadです。

### E2b: small-scaffold grouping → **REJECT**

data-contract gate（対E2a）: **ADOPT_CANDIDATE**。差分975件はすべてhelper自身の
保守的な`>= 35.01`分岐に収まり、unresolved 0、filter flip 0、非QD差分0です。
E2aと同じgateで独立に再分類した結果も同一で、directly explained 975件（100%）、
reconstruction-limit consistent 0件です。E2aと異なりreconstruction-limitの議論自体が不要でした。
本比較975件での最小出力QDは20.29、QD<2境界からの実測marginは18.29です
（本比較での実測値であり、jitterが常にQD=2から安全という一般論ではありません）。
**科学的な問題による却下ではありません。**

architecture gate（対E2a）:

| 指標 | E2a | E2b | 差 |
| --- | --- | --- | --- |
| **GenomicsDBImport peak RSS** | 1.76 GiB | **6.57 GiB** | **+272%** |
| workspace files | 4,611 | 4,515 | **-2.1%** |
| interval count | 53 | 29 | -45.3% |
| GenotypeGVCFs peak RSS | 7.43 GiB | 7.35 GiB | -1.1% |
| import wall合計 | 5.34 h | 5.28 h | -1.1% |
| elapsed | 4.09 h | 4.07 h | -0.4% |

grouped task単体（25 contig、合計818,535 bp）:

| 指標 | E2a（25 task） | E2b（1 grouped task） |
| --- | --- | --- |
| import wall | 40.2 s | 13.3 s（-26.9 s） |
| import peak RSS | 1.14 GiB | **6.57 GiB** |
| workspace files | 2,175 | 2,079（-96） |
| serial reader fallback | なし | **発生** |

決定的な観測: grouped taskの6.57 GiBは**run全体のimport peak**となり、
20 Mb window（最大1.76 GiB）を上回ります。しかもこのtaskが担当するのは合計818,535 bpで、
20 Mb window 1本の1/24未満です。GenomicsDBImportのmemoryは総塩基数ではなく
**同時に書き込む個別contig array数**に追随するため、25個の小scaffoldをまとめる方が
1 chromosomeの20 Mbを取り込むより3.7倍のmemoryを要します。

groupingの動機であったfile/task増加の抑制は実現しませんでした。
file数は2.1%（4,611中96）、import wallは約19,000秒中26.9秒の削減にとどまる一方、
run全体のimport peak RSSは272%増加し、grouped taskがimport側の律速になります。
さらにGATKはこのtaskでserial reader initializationへfallbackし、要求した8 reader threadsを失います。
96 fileの削減のために3.7倍のimport memoryとreader並列性を失い、grouping logicと
順序制約の複雑性を恒常的に抱えるのは割に合いません。

将来、small scaffoldが十分多くtask/file増加が実際の運用制約になる場合、
またはGenomicsDBImportのper-array memory挙動が変わる場合は再評価します。

### 確定したinterval strategy

```
ADOPT : 20 Mb chromosome window splitting + small scaffoldは個別task（E2a）
REJECT: small-scaffold grouping（E2b）
```

本referenceでのinterval数は53です。E2bのbaselineはE2aなので、grouping却下によって
window splittingの採用根拠は影響を受けません。

20 Mbは実測前にprotocolで固定した値であり、結果を見て探索していません。最適値の主張はしません。
いずれの数値も51検体・本reference・GATK 4.6.2.0での実測であり、327検体の挙動は示しません。

- [interval strategy decision](evidence/issue45/interval_strategy_decision.json)

## 2026-09-13: E3 ReblockGVCF → **REJECT**（既存判定の記録）

判定は`3251645`で確定済みです。本節は既存evidenceを要約するだけで、E3の再実行や判定の再審はしていません。

E3は29 interval taskとReblockをすべて完了しましたが、gathered callsetのoutput validationで停止しました。
51 sample列が丸ごと`.`で出力され、`joint_integrity`がdiploid検査で拒否したためです。
validatorはfail-closedに正しく動作しました。当時のメッセージはploidy errorとだけ表示していましたが、
これは`0a8f9a0`で、全列欠損、ploidyを持たない`.`、allele数不一致を区別して報告するよう改善しました。
いずれも従来どおり拒否します。

| 区分 | 件数 |
| --- | --- |
| representational change（phasing collapse、no-callの全列欠損表記） | 4,073,660 |
| hom-ref call → missing | 86,612,301 |
| called non-reference genotype → missing | 19,566 |
| variant set change（E2bのみ 46,135 / E3のみ 130,776） | 176,911 |
| AN changed records | 12,831,569 |
| AC changed records | 28,037 |
| FILTER changes | 0 |

判定理由: **the storage-saving trade-off is incompatible with this repository's current exact
scientific/data-contract requirement**。ReblockGVCFが誤っているという判断ではありません。GQ bands 20/100の
設計どおり参照信頼度の解像度を捨てた結果です。gVCFサイズ76.9%削減（48.03 GB → 11.11 GB）と
好ましい資源傾向は記録しますが、判定には含めていません。

- [ReblockGVCF decision evidence](evidence/issue45/reblock_decision.json)

## 2026-09-14: E4 consolidate → **REJECT**

E2bとの差分は`--consolidate true`だけです。plan、sample-name-map、batch size、入力、reference、
GATK image、allocation、並列数は同一です。

### data-contract gate: ADOPT_CANDIDATE

| 指標 | 結果 |
| --- | --- |
| shared / identical records | 13,219,170 / 13,219,170 |
| different records / E2bのみ / E4のみ | 0 / 0 / 0 |
| record / GT / accounting / shared header checksum | すべて一致 |
| sample order / contig order | 一致 |

partitionが変わらないため、E2a/E2bで見られたQD jitterすら生じません。E2b・E4のgathered fileは
byte列が異なりますが、header diffとrecord本体hash（両者`1d367557…`）により、差は
`##GATKCommandLine` 2行（`--consolidate`値、map path、日時）だけであることを確認しました。

### architecture gate: REJECT

同一taskで比較できる28本のchromosome window（E2a・E2b・E4で同一task）での比較です。

| 指標 | E2b | E4 | 差 |
| --- | --- | --- | --- |
| **GenomicsDBImport wall合計** | 18,997 s | 27,278 s | **+43.6%** |
| import wall差（task別） | — | — | 中央値+43.6%（300 s超の全taskで+40.7〜46.4%） |
| GenomicsDBImport peak RSS | 1.73 GiB | 2.43 GiB | 16 GiB割当内 |
| GenotypeGVCFs wall合計 | 21,691 s | 21,320 s | -1.7%（task別 -3.4〜+0.1%） |
| GenotypeGVCFs peak RSS | 7.35 GiB | 7.31 GiB | -0.5% |
| workspace bytes | 55.10 GB | 55.35 GB | +0.45% |
| workspace files | 2,436 | 1,316 | -1,120（windowあたり87 → 47） |
| experiment elapsed（全29 task、validation込み） | 4.07 h | 4.84 h | **+19.1%** |

retry 0、OOM 0。E4実行中のhostはMemAvailable最小99.98 GiB、swapは0.75 MiBで変化なし、
storage空き最小2.03 TBでした。

コストの所在は明確です。`NC_068975.1:1-20000000`では、batch 1（15分12秒）とbatch 2の1検体取込み（約20秒）が
両runで同じでした。E4はその後`Consolidating GenomicsDB array`を記録し、batch 2完了まで6分55秒を
要しています。import overheadはすべて、読込み後に配列を1 fragmentへ書き直す工程に由来します。

判定: consolidationの代償は実測でimport wall +43.6%、elapsed +19.1%です。見返りはGenotypeGVCFs wall -1.7%、
ceilingが懸念されるprocessへのmemory改善なし、window側のfile数-1,120にとどまります。本構成でfile数は
scheduling overheadであり、制約として観測されていません。**採用しません。** productionのmoduleは
flagを渡しておらず、GATK既定のconsolidate=falseのままです。callsetはどちらでも変わらないため、
これは純粋な資源判断です。

pinned GATK 4.6.2.0のdocは「100 batchを超える場合にconsolidateを使う」とし、効果は
「top Java layerのoverheadで目立たない可能性がある」としています
（[GenomicsDBArgumentCollection](https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/genomicsdb/GenomicsDBArgumentCollection.java)）。

**限界**: 51検体では配列あたり2 fragmentしか測っていません。327検体では`ceil(327/50)=7` fragmentになりますが、
これは算術上の帰結で実測ではありません。7 fragmentでの読込み側の利得と、327検体規模での書き直しコストは
未測定であり、327 projectionではassumptionとして扱います。7はpinned guidanceの100 batchより大幅に少なく、
実測にもdocにもconsolidateを支持する材料はありません。E2bのgrouped task（25 scaffold）では
file 2,079 → 1,079、import peak 6.57 → 5.39 GiBでしたが、groupingはREJECT済みで採用planへ
転用できないため、判定に含めていません。

再評価の条件: 配列あたりのbatch数がguidanceへ近づく場合（incremental importや大幅に小さいbatch size等）、
または327検体でfragment数に起因するGenotypeGVCFs読込み問題が観測された場合です。

- [consolidate decision evidence](evidence/issue45/consolidate_decision.json)
- [E0〜E4 sanitized metrics（全6 experiment）](evidence/issue45/benchmark_results.json)
- [suite evidence freeze manifest](evidence/issue45/suite_freeze.json)

### 確定したJoint Genotyping architecture（51検体実測）

| 要素 | 判定 |
| --- | --- |
| `--sample-name-map` | ADOPT（E1、sample order条件は327検体で解消済み） |
| 20 Mb chromosome window splitting | ADOPT（E2a） |
| small-scaffold grouping | REJECT（E2b） |
| ReblockGVCF | REJECT（E3） |
| `--consolidate` | REJECT（E4） |
| `genomicsdb_batch_size` | 50のまま（全36 / 53 / 29 intervalで50 + 1のbatchingを実測） |
| GenotypeGVCFs allocation | 16 GiBでは不足（E0 OOM）。20 Mb windowでの実測最大は7.43 GiB（E2a） |

## 2026-09-14: sample数scaling補助測定（13 / 26 / 51）

E0〜E4はすべて51検体での測定で、sample数方向の情報を持ちません。そこで補助測定を行いました。
これはどのone-factor判定の入力でもありません。

| 固定条件 | 値 |
| --- | --- |
| interval | `NC_068975.1:1-20000000`（suite E2b/E4の`interval_0017`と同一） |
| subset | validated manifest順の決定的prefix、13 ⊂ 26 ⊂ 51（nesting検証済み） |
| reference / GATK / allocation | suiteと同一（`execute_experiments.py`はsuite lineageとbyte一致） |
| 並列 | 各level単独、suite controller終了後に実行（51検体も単独で再測定） |
| helper | `0a8f9a0` |

### observed

| samples | variants | genotype cells | GenotypeGVCFs RSS / wall | GenomicsDBImport RSS / wall | workspace | batches |
| --- | --- | --- | --- | --- | --- | --- |
| 13 | 181,326 | 2.36 M | 1.27 GiB / 249 s | 1.08 GiB / 275 s | 0.80 GB / 47 files | 1 |
| 26 | 273,562 | 7.11 M | 2.11 GiB / 453 s | 0.80 GiB / 527 s | 1.55 GB / 47 files | 1 |
| 51 | 792,503 | 40.42 M | 7.46 GiB / 1,009 s | 1.31 GiB / 901 s | 2.71 GB / 87 files | 2 |

全levelでretry 0、sample order・contig order・AC/AN accountingはvalid、swap増加は最大0.25 MiBでした。
51検体を単独で再測定した出力は、suite E2bの同じintervalとrecord本体SHA256が一致しました。
resource値は単独実行のほうがwallで3.5〜5.0%速く、peak RSSの差は5%以内です。suiteの3並列はこのintervalを実質的に歪めていません。

### cohort構成との交絡

26→51でsample数は1.96倍ですが、variant数は2.90倍、GenotypeGVCFs RSSは3.54倍、出力sizeは4.19倍に跳ねます。
51検体の出力を解析した結果は次のとおりです。

- manifest位置33〜50の7検体は、non-reference callが16万〜22万（他は4.3万〜8.4万）で、大半がhom-altでした。
  referenceから大きく離れた系統です。公開metadataには事前にこれを識別できる属性がありません。
- 792,503 siteのうち、最初のnon-reference carrierが1〜13番は213,235、14〜26番は87,400、27〜51番は491,868でした。

**このprefix系列ではsample数と構成が同時に変わるため、Nだけにfitした指数はsample数の性質ではありません。**
GenotypeGVCFs memoryをNだけで外挿することはしません。

### empirical trend

GenotypeGVCFs peak RSSは、genotype cell数（samples × 出力variants）に対して、独立な3つの推定値が一致しました。

| 推定 | 傾き（GiB / 百万cell） |
| --- | --- |
| E0の11 chromosome（N=51固定、contig間でvariant数が変化） | 0.188（R² 0.951。同じ点をinterval長で回帰するとR² 0.933） |
| scaling 13→26 | 0.177 |
| scaling 26→51 | 0.161 |

wall timeは構成の影響を受けにくい指標です。E0ではinterval長でR² 0.992、variant数でR² 0.870でした。
Nに対しては、genotype wallの327/51比が6.46〜7.08、import wallが5.12〜6.00、workspace bytesが5.35〜6.08です（linear / power / 上側pair傾きの3形式）。
import peak RSSにNとの傾向はありません（batch sizeで上限が決まる）。workspace fileは配列あたり`7 + 40 × batch数`でした。

- [sample scaling evidence](evidence/issue45/sample_scaling_20260914.json)

## 2026-09-14: 327検体resource envelope

`benchmarks/issue45/project_resource_envelope.py`が、commit済みevidenceだけから
[resource_envelope_327.json](evidence/issue45/resource_envelope_327.json)を再導出します。
unit testはcommit済みJSONと再導出結果の一致を検査します。scriptは判定を行いません。

### observed（51検体、本host）

- 採用plan（E2a、53 task）: import wall合計19,216 s、genotype wall合計21,752 s、elapsed 4.09 h、workspace 55.12 GB / 4,611 files、gathered callset 4.86 GB
- 入力gVCF: 48.03 GB（1検体0.40〜1.28 GB）。FASTQ: 選定51検体207.4 GB、全327 runの公開size合計1.55 TB
- 上流生成（Joint Genotypingの範囲外、参考値）: 17.1 h、272.9 CPU-h、中間ファイル全保持で正味780.7 GB
- host: 32 CPU / RAM 123.5 GiB / `/data` 3.94 TB（現在の空き2.03 TB）
- production契約: GenotypeGVCFsは16 GB、retry時32 GB、per-contig、cleanup設定なし

### assumptions

| ID | 内容 |
| --- | --- |
| A1 | GenotypeGVCFs RSSは、最大観測task（88.1 M cells）を超えてもcell数に線形 |
| A2 | 327検体でのvariant数（未測定）を4 scenarioで表す: S0 = site増加なし（下限）、S1 = Watterson a_nによる増加、S2 = 14〜26番のdiscovery rateで減衰せず増加、S3 = 27〜51番のrateで減衰せず増加（悲観側） |
| A3 | 測定windowのvariant増加率がgenome全体に当てはまる（51検体時点で平均より高密度なので、memory側は保守的） |
| A4 | 7 batchのimportは1〜2 batchの観測と同じ挙動 |
| A5 | 同一host・同一I/O・3並列 |
| A6 | 上流の生成量はFASTQ bytes比またはsample数比で増える |
| A7 | 中間ファイルは51検体生成時と同様にすべて保持される |

### projected range

| 項目 | S0 | S1 | S2 | S3 |
| --- | --- | --- | --- | --- |
| window variants（51検体比） | ×1.00 | ×1.42 | ×3.34 | ×7.85 |
| GenotypeGVCFs RSS / 20 Mb window | 42.6–48.6 GiB | 59.9–68.8 GiB | 140–163 GiB | 328–382 GiB |
| productionのper-contig最長contig | 94–107 GiB | 132–151 GiB | 307–356 GiB | 717–835 GiB |
| 最大観測cell数に収まるwindow長 | 6.8 Mb（72 task） | 4.8 Mb（96 task） | 2.0 Mb（224 task） | 0.87 Mb（523 task） |
| callset（per-interval copy込み） | 45–62 GB | 60–88 GB | 122–208 GB | 246–489 GB |
| full run保持storage | 6.90–7.82 TB | 6.91–7.85 TB | 6.98–7.97 TB | 7.10–8.25 TB |

scenarioに依存しない項目:

- Joint Genotyping task時間: import 27.3–32.0 h、genotype 39.0–42.8 h
- elapsed: 3並列で23.8–26.9 h、2並列で35.7–40.3 h
- workspace 295–353 GB、配列あたり287 files（E2a planで15,211 files）
- 入力gVCF 308–360 GB
- 上流生成: 保持storage 5.01–5.85 TB、1,750–2,046 CPU-h、同一allocationで110–128 h

### uncertainty

- 327検体でのvariant数は未測定です。S0〜S3は約8倍の幅があり、memoryとcallsetの幅の大部分を占めます。
- A1は観測範囲外への外挿です。S0でも20 Mb windowには最大観測の約3倍のcell数が必要です。
- peak RSSはheap 26 GiB下のresident setで、live setの実測ではありません（ただしE0のOOMは、ceiling付近でlive setが追随することを示しています）。
- 51検体はaccession順の選定で、327検体からの無作為標本ではありません。divergent群が327検体中に占める割合は不明です。
- 7 batchのimportは1〜2 batchからの外挿です。
- wall timeの幅は1 windowと3つのtrend形式に基づきます。構成がgenotype時間を押し上げる可能性は、memoryほど大きくはありませんが残ります。

## 2026-09-14: full 327-sample run Gate → **NO-GO**

| 基準 | 結果 |
| --- | --- |
| GenotypeGVCFs memory | **FAIL**: S0（site増加なし）でも20 Mb windowが32 GiB tierを超え、launch gate内に最大2 task。productionのper-contigはretry時32 GBの約3倍 |
| storage | **FAIL**: 中間ファイル全保持で6.9–8.25 TB。volume全体の3.94 TBを超える |
| time | 非決定的: Joint Genotypingは3並列で約1日、上流を含めても日単位。ただし並列数はmemoryが許す場合に限る |
| retry / swap | PASS（51検体でretry 0、OOM 0、swap増加は最大1 MiB） |
| scientific / data contract | PASS（採用した変更はrecord一致、または説明済みのQD jitterのみ。callを変えたReblockGVCFは不採用） |
| lineage | 51検体はPASS。327検体は全gVCFを単一SHAで生成する必要があり、未実施 |

CONDITIONAL GOではなくNO-GOとする理由です。CONDITIONAL GOは、envelopeが収まり不確実性だけが残る場合に使います。
今回は投影範囲の最も楽観的な端ですでにmemory tierとstorage volumeを満たさず、productionも採用architectureを実装していません。
これはprotocolの「resource contractを満たせない」に該当します。
cells modelの外挿（A1）がmemoryを過大に見積もっている可能性はあります。ただしS0で20 Mb windowが32 GiBに収まるには、1/4〜1/3の過大評価が必要です。
またstorageのFAILはmemory modelに依存しません。

本判定は、51検体で決めたarchitecture判定を再審するものではありません。327検体をgenotypeできないとも主張しません。
「このpipelineとhostで、現行の契約のまま327検体runを開始すべきではない」という判定です。

### 再評価の条件

1. **R1**: productionが`--sample-name-map`とchromosome window splittingを採用する（small scaffoldは個別task、consolidateはoff）
2. **R2**: window長を固定の20 Mbではなく、taskあたりのgenotype cell予算から決める。最大観測task（88.1 M cells）以内に収めるなら、S0で6.8 Mb以下、S1で4.8 Mb以下（72〜96 chromosome task）
3. **R3**: GenotypeGVCFs allocationをその予算と既定のnative reserveから決め、並列数 × tierを110 GiBのlaunch gate内に収める
4. **R4**: 6.9–8.25 TBの保持量に対応するretention/cleanup方針またはvolumeを用意し、生成中に1検体あたりの保持byteを実測で確認する
5. **R5**: 327 gVCFを単一SHAで生成した後、1 windowでGenotypeGVCFsを単独・十分なceilingで測定し、A1/A2を実測で置き換える。あわせて7 batchのimport（A4）を確認してから、full Joint Genotypingを見積もる

R1〜R4を満たしR5を実測した後、同じprotocol ruleで再評価します。自動的にGOにはなりません。

- [327 gate decision](evidence/issue45/gate_327_decision.json)
