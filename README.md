# adzuki-snp-pipeline

公開のホールゲノムショットガンシーケンシング(WGS)データから、アズキ(*Vigna angularis*)の
再現可能なSNPコーリングパイプラインを構築するための研究開発リポジトリです。

このリポジトリには、ペアエンドWGSの前処理からサンプル単位のGVCF生成、複数サンプルの
Joint Genotyping、設定可能なハードフィルタリング、PASS抽出、variant QC、そしてゲノミック
セレクション(GS)パネル生成までを行う実行可能なNextflow DSL2ワークフローが含まれています。
5検体の実データコホート(Issue #26)をSeedcore-01実機上で実行し、Joint Genotyping・
GSパネルまでの一連のend-to-end動作と再現性証跡を確認しました(詳細は
[実データコホートE2E検証](#実データコホートe2e検証issue-26)を参照)。あわせて、履歴として
手動実行した単一サンプルのSNPコーリング試行(pilot)も残しています。パイプラインレベルの
nf-testスイートがJoint Genotyping fixture契約・GenomicsDBImportのbatch-size/メモリ契約・
参照ゲノムFAI/dict契約をカバーしています。MultiQC集約、全モジュールを網羅したnf-test、
20〜30検体・327検体スケールでの検証は今後の課題として[Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1)
で追跡しています。

この研究はFlorigen AIの長期的な農業AI活動を支える植物遺伝学・バイオインフォマティクス研究の
一部です。ゲノム予測からPhysical AIへの直接的な開発パスを意味するものではありません。

---

## 現在の成熟度

| 能力 | 状態 | 根拠・制約 |
| --- | --- | --- |
| 手動単一サンプルSNPコーリング(pilot) | 1回実行済み | SRR29909135を手動で処理した結果であり、自動テストによる再現はまだ行っていない |
| 手順の文書化 | 記録済み | コマンドは本文書に記録済みだが、ソフトウェアバージョンと実行パラメータは完全には固定されていない |
| Nextflow DSL2ワークフロー | フィルタ済みコホートVCF・variant QC・GS SNPパネルまで実装済み | strict構文パーサ対応の前処理・マッピング・重複マーク・GVCF生成・Joint Genotyping・ハードフィルタリング・PASS抽出・QC処理・GSパネル正規化/行列構築が利用可能 |
| Variant calling / genotyping | syntheticデータで機能検証済み、実データ5検体コホートで動作確認済み | HaplotypeCaller、contig単位のGenomicsDBImportとGenotypeGVCFs、reference順のGatherVcfsが接続済み。Issue #26で5検体・36contigの実データを通した(下記参照) |
| Variant filtering / QC | syntheticデータで機能検証済み | SNPとindelを分離し、設定可能なハードフィルタを適用し、PASSレコードを抽出し、raw/filtered/PASSの各段階でQC成果物を生成する。FILTER値・annotation別のevaluable rate会計はannotation欠損としきい値未達を区別し、cohort全体のreconciliationは`raw/snp`+`raw/indel`が`raw/all`からどれだけ乖離しているかを報告する。しきい値の妥当性は実コホートでは未検証 |
| GSパネル(genomic selection) | syntheticデータで機能検証済み、実データ5検体で動作確認済み | `raw/all`を正規化(`bcftools norm -m-`)し、split後のREF/ALT形状で再分類し、独立したlineageで再フィルタし、sample/variant metadataとfull-lineage record会計、再現性manifestを備えたdosage matrixへ変換する([`docs/gs_panel_data_contract.md`](docs/gs_panel_data_contract.md)参照)。MAF/コールレートフィルタ、LD pruning、imputationはスコープ外。しきい値の妥当性は実コホートでは未検証であり、現時点でこのパネルを読み込む下流リポジトリは存在しない |
| 設定可能な参照ゲノムbundle | 実装済み | 互換性のあるprebuilt indexを受け付けるか、FASTA・sequence dictionary・BWA-MEM2 indexを生成する |
| 複数サンプルJoint Genotyping | syntheticデータで機能検証済み、実データ5検体・36contigで動作確認済み(Issue #26) | syntheticでは2サンプル・2contigでGenomicsDBImportベースのJoint Genotypingが完了し、両サンプルが両contigをカバーするため各deterministic SNP locusは`1/1`/`0/0`の非欠損ペア(`AC=2;AN=4;AF=0.5`)に解決される。実データでは5サンプル・36contigのLongxiaodou 4参照ゲノムで動作を確認したが、20〜30検体・327検体規模の検証はまだ行っていない(Issue #26/#11) |
| 読み取り前処理・QC | MultiQCなしで実装済み | raw/trimmed FastQC、ペアエンドfastp、mappingログ、重複metrics、SAMtools QCを生成する |
| パイプラインレベルテスト | 部分実装 | syntheticなDocker smoke testが3リードグループ・2サンプルGVCF・両locusで確信度の高い`1/1`/`0/0`遺伝子型を持つ2つのraw SNP・ハードフィルタannotation・indexed PASS出力・7種のvariant QCタスク・38個のQC成果物・GSパネルのempty-panel契約を検証する。[nf-test](https://www.nf-test.com/)によるpipeline-level testがこの契約を自動化し(`tests/pipeline/adzuki_snp_pipeline.nf.test`)、haploid(`sample_ploidy=1`、`enable_gs_panel=false`)のend-to-end遺伝子型契約とGSパネル有効時の`sample_ploidy=1`拒否(Issue #20)、mappingタスク数と重複マーク件数の厳密一致・負の`optical_duplicate_pixel_distance`拒否・prebuilt BWA-MEM2 indexパス(Issue #8)、`genomicsdb_batch_size`の伝播とXmx比率・prebuilt/生成両方の参照bundle・contig順序不一致の拒否(Issue #11)を検証する追加testを含む。module-level nf-testは実際にpinされたcontainerに対して、`GATK_SELECTVARIANTS`がMIXED型レコードをSNP/indel双方の選択から除外すること(`tests/modules/gatk_selectvariants.nf.test`)、`GS_NORMALIZE_VARIANTS`がMIXEDレコードを正しい遺伝子型再割当てで分割すること(`tests/modules/gs_normalize_variants.nf.test`)、`BWA_MEM2_MEM_SORT`が中間`.sam`を生成せず文書化されたCPU分割を適用すること(`tests/modules/bwa_mem2_mem_sort.nf.test`)、`VALIDATE_REFERENCE_CONTIGS`がFAI/dictの名前・長さ・順序不一致を検出すること(`tests/modules/validate_reference_contigs.nf.test`)、`GATK_GENOMICSDBIMPORT`/`GATK_GATHERVCFS`が単一sample入力(Nextflowのpath入力がList→scalarへ暗黙変換される境界条件)を正しく扱うこと(`tests/modules/gatk_genomicsdbimport.nf.test`、`tests/modules/gatk_gathervcfs.nf.test`)を確認する。全モジュールを網羅するnf-testは[Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1) #Gで計画中 |
| Functional CI | syntheticフィクスチャpathで実装済み | `.github/workflows/test.yml`が`nextflow lint .`と結合nf-testスイート(pipeline-level + module-level、`-profile test,docker`)をmain向けのpush/pull requestごとに実行する。`.github/workflows/lint.yml`はリポジトリ構造のみを個別に検査する。実データコホートのCIはスコープ外 |
| 実データコホート検証(Issue #26) | 5検体E2E検証を実施 | 公開WGSデータ5検体(同一BioProject)をSeedcore-01実機でQC→mapping→重複マーク→gVCF→Joint Genotyping→hard filtering→variant QC→GSパネルまで通した。詳細は[実データコホートE2E検証](#実データコホートe2e検証issue-26)を参照。20〜30検体・327検体規模の拡張は未実施(意図的にスコープ外) |
| Base Quality Score Recalibration (BQSR) | 意図的に除外 | 検証済みのknown-sitesリソースが存在しない。[設計判断](#設計判断)を参照 |
| 本番利用 | 非対応 | これは実験的な植物研究リポジトリである |

本文書中の図と大規模なvariant件数は、単一サンプルpilotに由来する歴史的な結果であり、実行可能な
ワークフローによって再現されたものではありません。synthetic fixtureを用いたワークフロー経路は、
決定論的なexpected variantsを伴うクリーンなDocker smoke testを別途完了しています。

---

## データと参照ゲノム

### シーケンシングデータ

| 項目 | 詳細 |
| --- | --- |
| BioProject | [PRJNA1138464](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1138464) |
| 単一サンプルpilotの検体 | SRR29909135 |
| 実データコホート検証(Issue #26)で使用した検体 | SRR29909135、SRR29909069、SRR29909072、SRR29909067、SRR29909073(5検体、詳細は[実データコホートE2E検証](#実データコホートe2e検証issue-26)参照) |
| 関連論文 | [Chien et al. 2025, *Science* 388: eads2871](https://doi.org/10.1126/science.ads2871) |
| 論文が報告する公開データ | 327検体のWGS再シーケンシングと357検体のDArT-seqデータ |
| このリポジトリで検証済みのスコープ | WGS 5検体(327検体中) |

同論文は357検体分のDArT-seqデータも報告していますが、このリポジトリの現在および計画中のvariant
callingスコープはWGSのみであり、RAD-seq/DArT-seqワークフローは実装・評価していません。

### 参照ゲノム

単一サンプルpilotおよびIssue #26の実データコホート検証は、いずれも以下の独立した参照アセンブリを
使用しています。

| 項目 | 詳細 |
| --- | --- |
| Accession | [GCF_016808095.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016808095.1/) |
| Assembly | ASM1680809v1 |
| Cultivar | Longxiaodou 4 |
| Assembly span | 447.8 Mb(実測: 448,362,642 bp、contig 36本。内訳は[Joint Genotyping scale hardening](#joint-genotyping-scale-hardeningissue-11)参照) |
| 論文 | [Li et al. 2024, *Scientific Data* 11:1074](https://doi.org/10.1038/s41597-024-03911-y) |

アズキでよく引用される約540 Mbという値は、cultivar Shumariのk-merベースのゲノムサイズ推定です
([Sakai et al. 2015, *Scientific Reports* 5:16780](https://doi.org/10.1038/srep16780))。
Longxiaodou 4についてLi et al.は21-mer解析により464.9 Mbと推定しており、447.8 Mbのアセンブリは
その推定値の96.32%に相当します。これらの値はcultivarと推定手法が異なるため、互換なものとして
扱ってはいけません。

シーケンシングデータと参照アセンブリは異なる研究・異なる遺伝的背景に由来します。実行可能な
ワークフローは参照ゲノムを明示的な設定可能入力として扱っており、異なる参照ゲノムに対する結果を
互換なものと仮定してはいけません。

---

## Nextflow WGS Variant-Calling ワークフロー

Issue [#4](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/4)、
[#6](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/6)、
[#9](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/9)、
[#13](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/13)が、ペアエンドWGSリードから
フィルタ済みコホートVCFとvariant QCまでの入力契約と実行可能ワークフローを確立しています。

現在のワークフローは以下を行います。

- `nf-schema` 2.8.0でパイプラインパラメータを検証する
- samplesheetの構造・必須値・パス・read group一意性を検証する
- read 1/read 2に同一パスが指定された場合、およびFASTQファイルが複数read groupで再利用された場合に拒否する
- 1つの生物学的サンプルに対して複数のread group・シーケンシングレーンを許可する
- ペアエンドfastpトリミングの前後でFastQCを実行する
- 互換性のあるFASTA・sequence dictionary・BWA-MEM2 indexを生成または受け付ける
- 各read groupを単一のpiped `BWA_MEM2_MEM_SORT`タスクで独立にmapping・座標ソートする(BWA-MEM2のアラインメントストリームを`samtools sort`へ直接pipeし、中間`.sam`ファイルは一切生成しない)。`ID`・`SM`・`LB`・`PL`・オプションの`PU`メタデータを保持する
- サンプル自身の期待read group数が揃った時点で、他サンプルの進捗を待たずにそのサンプルのread groupを最終BAMへmergeする([Read groupを意識した早期merge](#read-groupを意識した早期merge)参照)
- 重複レコードを削除せずGATK MarkDuplicatesでlibrary単位の重複をマークする。`optical_duplicate_pixel_distance`は設定可能
- BAM indexとSAMtools flagstat/stats/idxstatsレポートを生成する
- GATK HaplotypeCallerで生物学的サンプルごとに1つのindexed GVCFを生成する
- 参照ゲノムのcontigごとに1つのGenomicsDBワークスペースを作成し、全サンプルのGVCFを共同genotypingする
- contig単位のraw VCFを参照ゲノムのindex順にgatherし、indexed raw cohort VCFを作成する
- cohort VCFをindexedなSNP VCFとindel VCFに分離する
- GATK VariantFiltrationで設定可能なSNP/indel別のハードフィルタを適用する
- indexedなPASSのみのSNP/indel VCFを抽出する
- raw/filtered/PASSの各段階で`bcftools stats`を実行し、機械可読なcohortおよびサンプル単位のQCテーブルと人間可読なサマリーを生成する
- filtered SNP/indel VCFについてFILTER値・annotation別のevaluable rate会計を報告し、`raw/snp`と`raw/indel`のレコード数を`raw/all`に対してreconcileする
- 決定論的なexpected SNPを伴う再配布可能なsynthetic functional-testデータセットを提供する

まだ実装していないのは、MultiQC実行、全モジュールを網羅したnf-test、そして20〜30検体・327検体
規模でのコホート検証(Issue #26/#11で計画)です。5検体規模の実データE2E検証は完了しています。

### 必要要件

- Bash 3.2以降
- Java 17以降
- Nextflow 26.04.6
- コンテナ実行用のDocker

Nextflow 26.04以降はデフォルトでstrict構文パーサを使用します。グローバルにインストールされた
launcherを変更せずにテスト済みバージョンを選択できます。

```bash
NXF_VER=26.04.6 nextflow -version
```

初回実行時にpin済みの`nf-schema`プラグインをダウンロードする場合があります。synthetic FASTQと
参照fixture自体はこのリポジトリに保存されています。

各processのcontainerはimage tagとmanifest digestの両方でpinされています。BioConda BWA-MEM2
image(単独の`BWA_MEM2_INDEX`processでのみ使用)はversion 2.3としてtagされていますが、同梱の
実行ファイルはversion 2.2.1を報告します。これは上流の2.3リリースが以前のinternal version
stringを保持したままだったためです。
([BioConda recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/bwa-mem2/meta.yaml)、
[upstream issue #283](https://github.com/bwa-mem2/bwa-mem2/issues/283)参照)

`BWA_MEM2_MEM_SORT`(mapping + sort)は代わりに結合された
[Seqera Wave](https://seqera.io/wave/)のmulti-package image、
`community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113`をtagとdigest両方で
pinして使用しており、`docker run ... bwa-mem2 version` / `samtools --version`で直接検証済みです:
bwa-mem2 **2.2.1**(上記2.3-tag imageと同じ実バイナリ)とsamtools **1.22.1**。このパイプラインの
他所すべて(`SAMTOOLS_MERGE`、`SAMTOOLS_INDEX`、`SAMTOOLS_QC`、`SAMTOOLS_FAIDX`)でpinされている
samtools **1.24**とbwa-mem2を組み合わせたmulled/Wave imageは見つかりませんでした。この
1.22.1対1.24のギャップは意図的な、文書化されたこのcontainer選択の制約であり(Issue #8)、
見落としではありません。このパイプラインのsynthetic fixtureは、生成されるBAMが期待通りの
coordinate-sort headerを持つこと、そして生成/prebuiltの両indexパス(いずれも同じsamtools
1.22.1でソート)が1.24でpinされた下流stageを通した後にbyte-identicalなraw cohort VCFを
生成することを確認しています。これはsamtools 1.22.1自身の出力を1.24自身の出力と比較したもの
ではありません。samtoolsバージョン間のBGZF相互運用性は、同一のsort tie-break順序・圧縮・
バグ修正を保証するものではなく、SAM仕様は同一RNAME/POSのレコード順序を未規定のままにしている
ため、実データにおける1.22.1と1.24自身の間の挙動差は未検証です。Issue #8 Phase 5の実参照ゲノム
実行(詳細は[`docs/mapping_real_reference_profile.md`](docs/mapping_real_reference_profile.md))
は実WGSデータ上で同じ`@PG`バージョンと`samtools quickcheck`/`flagstat`がクリーンなBAMを
再確認しましたが、1.22.1と1.24自身のsort出力同士を比較したわけではありません。結合image
(2つの独立したcontainerではなく)が必要な理由は、bwa-mem2とsamtools sortが1つのshell pipe内で
*同時に*動作し、1タスク分のCPU/メモリ予算を共有する必要があるためです(下記参照)。

### Mapping + sortのリソース契約

`BWA_MEM2_MEM_SORT`のpipe内の両ツールは同時に動作し、1タスク分の`cpus`/`memory`割り当て
(`process_mapping`ラベル)を共有するため、どちらか一方にタスクの全予算を割り当てるとoversubscribe
になってしまいます。

- **CPU分割**: bwa-mem2が`task.cpus`の80%(floor 1)、samtools sortが残り20%(floor 1)を得ます。
  bwa-mem2のアラインメント/スコアリング処理がペア全体のCPUコストを支配し、スレッド数にほぼ線形に
  スケールする一方、samtools sortのスレッド追加による高速化は数スレッドを超えると小さくなります。
  **Issue #8 Phase 5**は、この分割を実際のLongxiaodou 4参照ゲノムと実際の約19.3xカバレッジの
  WGS検体1件に対して1台のマシン(Seedcore-01、[`docs/mapping_real_reference_profile.md`](docs/mapping_real_reference_profile.md)参照)
  上で測定しました: 変更を加えていない80/20分割はOOMやretryなく完走し、タスク全体のCPU使用率は
  bwa-mem2がペア全体のコストを支配しているという説明と整合していました。この単一のベースライン
  実行ではボトルネックや失敗は見つからなかったため、**変更せずに維持**しています。これは80/20が
  その1つの参照ゲノム/検体/マシンにおいて*安全*であることを確認したものであり、*最適*であることを
  証明したものではなく、比較対照となる別分割のベンチマークは行っていません。**Issue #26**では
  この同じ80/20分割を5検体・36contigの実データコホートに対して適用し、OOMやretryなく完走する
  ことを追加で確認しました(詳細は[実データコホートE2E検証](#実データコホートe2e検証issue-26)参照)。
- **メモリ分割**: `samtools sort -m`は*スレッドあたり*のメモリを指定します(省略時はタスクが実際に
  利用可能な量に関わらず768 MiB/スレッド)。この計算式は`task.memory`の50%をbwa-mem2用に確保し
  (そのフットプリントは参照ゲノム/index sizeにスケールするため、この式は実測に先立ってそれを
  予測しようとしません)、固定512 MiBをOS/ツールのオーバーヘッド用に確保し、残りをsort自身の
  スレッドに均等に分配します。以前のバージョンのこの計算式と異なり、`task.memory`/`task.cpus`の
  組み合わせがsortに最低1 MiB/スレッドを残せない(あるいは両ツールに各1スレッドを与えられない)
  ほど小さい場合、割り当てを黙ってclampしてタスクの実際のメモリ予算を超過させる代わりに、診断
  可能なエラーでタスクを即座に失敗させます(PR #25レビュー)。同じPhase 5実行はこの計算式を
  `task.memory=16 GiB`で測定し、peak RSSは`BWA_MEM2_MEM_SORT`で12.9 GB(約19.4%の余裕)、
  `BWA_MEM2_INDEX`で10 GB(37.5%の余裕)に達し、いずれもOOMは発生しませんでした。この計算式と
  `16.GB * task.attempt`のデフォルトは**変更せずに維持**しています。この1回の成功実行は、それを
  引き上げる・引き下げる根拠にはなりません。

両方の計算式は`modules/local/bwa_mem2_mem_sort.nf`にinline実装・文書化されています。
`BWA_MEM2_MEM_SORT`と`BWA_MEM2_INDEX`はそれぞれ`process_high`とも互いとも独立した専用の
resourceラベル(`process_mapping`、`process_bwa_index`)を持っており、より大規模あるいは異なる
形状の実データセット(Issue #11/#26)が今回の1回の実行では見つからなかったボトルネックを
明らかにした場合でも、個別にチューニングできます。

タスクスクリプトでは明示的に`set -o pipefail`を設定しています。bashのデフォルトの`-e`だけでは
shell pipe内の前段のexit codeをチェックしないため、これがなければbwa-mem2の失敗(segfault、
入力の破損)が、破損/切り詰められた入力に対して`samtools sort`がexit 0で終了することでマスク
されてしまう可能性があります。

OOM様のexit code(137/140/143)によるretryは、失敗した実行と全く同じメモリ制限で単純にretryする
代わりに、すべてのresourceラベルにわたって`task.attempt`でメモリをスケールします(例:
`{ 16.GB * task.attempt }`)。単純なretryでは、最初の試行が本当にOOM killされた場合に成功する
ことはできません。

### Read groupを意識した早期merge

Nextflowの`groupTuple()`は明示的な`size:`を指定しない場合、上流channel全体が閉じるまで
(=全サンプルの全read groupのmappingが終わるまで)、たとえあるサンプルの全read groupが
とっくにmapping済みであっても、完了した1つのサンプルグループすら出力せずにすべてのtupleを
バッファします。`main.nf`はmaterializeされたsamplesheetから直接各サンプルのread group数を
数え(`read_group_counts_by_sample`、channel処理が始まる前に一度だけ計算)、その数を
`groupKey(sample_id, expected_read_group_count)`(Nextflow自身が文書化しているこの
パターン)へ渡します。これにより`groupTuple()`は、他のサンプルの進捗に関わらず、そのサンプル
自身のread groupが揃い次第、そのサンプルのグループと`SAMTOOLS_MERGE`タスクを出力します。
`groupKey`自身の`.target`プロパティが後から元の`sample_id`を復元します。既存の決定論的な
BAM名順序でのmerge契約は変更されていません。

### Optical duplicate検出のためのread-name形式

synthetic FASTQのread nameは、以前の`<read_group_id>_<contig>_<number>`方式ではなく、CASAVA
1.8スタイルの7コロンフィールドIllumina形式(`SIM:1:FC1:<lane>:<tile>:<x>:<y>`)に従っています。
GATK MarkDuplicatesのデフォルト`READ_NAME_REGEX`は、optical対PCR重複の下位分類のためtile/x/yを
抽出するために正確にこの形式を期待しています。以前の名前は実際に確認されたMarkDuplicatesの
警告(`did not start with a parsable number`)を引き起こし、optical duplicate clusterを常に
0件と報告していました。これは基盤となるPCR/library重複FLAG自体(座標/CIGARベース)には影響
しませんが、その下位分類(library size推定)を黙って劣化させていました。laneはread group自身の
`_L<NNN>`接尾辞から解析され、tileとx/yのベースラインは`read_group_id`とcontig名の安定した
ハッシュから導出されます。これにより、laneの番号を共有するが異なるサンプルに属する2つの
read group(例: 現実的に1本の物理レーンへ多重化される`sample_a_L001`と`sample_b_L001`)が
異なる模擬tileへ配置され、2つのサンプルのDNAが同一の物理クラスタを占有することは決してない
という事実と整合します。正確な導出方法は`tests/scripts/generate_synthetic_data.py`の
`casava_read_name()`を参照してください。

### ワークフローの検証と実行

pin済みcontainerは現時点でLinux AMD64を対象としています。ネイティブLinux AMD64では`docker`を、
Apple SiliconでのDocker emulationによる機能テストでは`docker_amd64`を使用してください。emulation
profileはパフォーマンスベンチマーク用ではありません。

```bash
NXF_VER=26.04.6 nextflow lint .

NXF_VER=26.04.6 \
  nextflow run . \
  -profile test,docker
```

Apple Siliconでは:

```bash
NXF_VER=26.04.6 \
  nextflow run . \
  -profile test,docker_amd64
```

synthetic datasetには2本の5 kb contig、2つの生物学的サンプル、3つのread groupが含まれています。
`sample_a`と`sample_b`はいずれも両方のcontig上にリードを持ち、各サンプルは自身のdeterministic
SNP locusにalternate alleleを、もう一方のサンプルのlocusにreference支持リードを供給するため、
Joint Genotypingは欠損遺伝子型ではなく同一site上でnon-referenceコールとreferenceコールの両方を
生成します([Issue #12](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/12)参照)。
2つの`sample_a`のread groupは`library_id=lib_a`を共有し、`chrSynthetic1`上に意図的な
cross-lane重複フラグメントを含みます。成功した実行では`sample_a`の40リードすべてが保持され、
うち4リードが重複としてマークされます。`sample_b`の24リードには意図的な重複は含まれません。

fixtureはさらに`sample_a`用の`chrSynthetic1:1501 C>G`と`sample_b`用の
`chrSynthetic2:1601 A>C`にdeterministic SNPを符号化しています。各locusでキャリアサンプルは
`1/1`、もう一方のサンプルは(`./.`ではなく)`0/0`とコールされ、2サンプルcohortで
`AC=2;AN=4;AF=0.500`となります。期待されるref/alt allele、遺伝子型、正確なサンプル単位/site
depth、`AC`/`AN`/`AF`は`tests/data/variants/expected_variants.tsv`に記録されており、これが
下記のnf-testスイートが直接読み込んで(下限値としてではなく)厳密一致で検証する正規の契約です。
クリーンなDocker smoke実行は、`sample_a`、`sample_b`のサンプル列順で、これら2つのraw SNP
レコードのみを、想定外のvariantレコードなしに生成します。

デフォルトの`snp_filter_sor_max=3.0`では、両方のsynthetic SNPともSOR値が3.0を超えるため
(`chrSynthetic1:1501`で`SOR=4.407`、`chrSynthetic2:1601`で`SOR=3.912`)、`SNP_SOR_HIGH`
フィルタを受けます。そのためデフォルトのPASS SNP VCFは有効だが空であり、`--snp_filter_sor_max 10.0`
を指定した別途の許容的なsmoke実行では、このfixtureが強化される前と同様に両SNPがPASSとして
保持されます。もう一方のサンプルの各locusにreference支持リードを追加したことで、期待通り
`AC`、`AN`、`AF`、サンプル単位/site `DP`が変化し、raw cohortの遺伝子型欠損率は`0.5`から`0.0`へ
低下しましたが、`SOR`、`FS`、`MQ`、`QD`は数値上変化しませんでした。これらのannotationは
cohort全体のpileupではなくalternate保持サンプル自身のリードから計算されているように見え、
そのため同一site上にhomozygous referenceサンプルを追加してもこれらの値はシフトしませんでした。
この挙動はfilterラベリング・PASS抽出・空VCF処理を検証するものであり、どちらのしきい値が
生物学的に適切であることの証拠ではなく、上記のSOR/FS/MQ/QD不変性はこのGATKバージョンの
annotation挙動に関する経験的観察であり、他所で頼るべき文書化された保証ではありません。

これらのfixtureはワークフローの挙動、read-groupを意識した重複マーク、GVCF生成、Joint
Genotyping、hard filterの仕組み、variant-QC生成をテストするものであり、生物学的または本番
環境での検証を構成するものではありません。両サンプルとも、両方のdeterministic SNP locusで
確信度の高い非欠損遺伝子型を受け取るようになり、以前
[Issue #12](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/12)で追跡していた
reference対alternate遺伝子型のギャップを解消しました。

deterministic fixtureは以下で再生成できます。

```bash
python3 tests/scripts/generate_synthetic_data.py
```

### 自動化されたパイプラインレベルテスト

パイプラインレベルの[nf-test](https://www.nf-test.com/)は`tests/data/variants/expected_variants.tsv`
を直接読み込み、上記の遺伝子型・annotation契約をend-to-endで検証します: gVCFとraw cohort VCFの
存在、想定外のvariantなしの正確に2つのraw SNPレコード、サンプル列順、そして両locusでの
`GT`、`REF`/`ALT`、`AC`、`AN`、`AF`、サンプル単位/site `DP`の(下限値ではなく)厳密一致、加えて
PR #14のfiltering/QC回帰guard(デフォルトの`SNP_SOR_HIGH`フィルタリング、空のデフォルトPASS
SNP VCF、空のindel VCF)。またvariant-QC出力契約(36個の`qc/variants`成果物、raw cohortの
欠損率metrics、あるstageのレコード数が0の場合のサンプル単位欠損率`NA`、FILTER値の内訳と
annotationカバレッジ契約、type単位のレコード会計)も検証します。契約ファイルを文字列として
テスト内に複製せず直接読み込むことで、fixtureとそのtestが乖離することを防いでいます。

さらに2つのpipeline-level test(Issue #20)が同じsynthetic fixtureに対して`--sample_ploidy 1`で
再実行し、`sample_ploidy`/`enable_gs_panel`を検証します。1つは`enable_gs_panel=false`が
end-to-endで正しいhaploid遺伝子型(`GT`、`AC`、`AN`、`AF` -- diploidのものではなく専用の
`tests/data/variants/expected_variants_haploid.tsv`契約から読み込む)を解決し、
`variants/gs_*`/`gs_panel/`出力が一切生成されないことを検証します。もう1つは、無効な
デフォルトの組み合わせ(`sample_ploidy=1`で`enable_gs_panel`をデフォルトの`true`のまま)が
どのprocessも開始する前に(`workflow.trace.tasks()`が空)診断可能なエラーメッセージで失敗し、
出力ディレクトリが作成されないことを検証します。

メインの遺伝子型契約testはさらにIssue #8のmapping/重複マークhardeningを検証します: 正確に3つの
`BWA_MEM2_MEM_SORT`タスクと`SAMTOOLS_SORT`タスクが皆無であること(`workflow.trace.tasks()`)、
各最終BAMの`@HD`/`@RG`行(実際のBGZF圧縮BAMから`java.util.zip.GZIPInputStream`で直接decode、
テストランナーに`samtools`バイナリが存在することへの依存なし)、`sample_a`の40保持リードのうち
4件が重複フラグ付き、`sample_b`の24保持リードのうち0件、MarkDuplicates metricsファイルの存在、
そしてCLI相当のString `optical_duplicate_pixel_distance="2500"`で、その値が実際のGATKコマンド
ライン(metricsファイル自身に埋め込まれたコマンドラインheaderから読み込む)に到達すること。
さらに2つのtestが、負の`optical_duplicate_pixel_distance`がどのprocess実行前にも即座に失敗する
こと(Issue #20の負のtestを反映)と、完全にprebuilt(非生成)のBWA-MEM2 indexが`BWA_MEM2_INDEX`
を完全にスキップしながらパイプライン全体を完走することを検証します。専用のmodule-level
nf-test(`tests/modules/bwa_mem2_mem_sort.nf.test`)は、実際にpinされた結合containerに対して
`.sam`ファイルがタスク自身の作業ディレクトリに一切作成されないこと、そしてbwa-mem2/samtools
sortが(それぞれ自身の`@PG` headerから読み込んだ)文書化された非oversubscribeスレッド数を、
新規生成/prebuiltの両方のindexで受け取ることを追加で確認します。

Issue #11の`GATK_GENOMICSDBIMPORT`/`GATK_GATHERVCFS`のsingle-input回帰guardは
[Joint Genotyping scale hardening](#joint-genotyping-scale-hardeningissue-11)で詳述します。

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash

NXF_VER=26.04.6 ./nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test
```

Apple Siliconでは、`nf-test.config`で宣言されたprofileをDocker emulationへ上書きしてください。

```bash
NXF_VER=26.04.6 \
  ./nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test \
  --profile "test,docker_amd64"
```

GitHub Actionsは同じnf-testスイートに加え、下記のvariant QC summarizer自身のunit test、
[Filter and annotation-coverage QC](#filterとannotationカバレッジqc)で説明する
`GATK_SELECTVARIANTS`module-level nf-test、`GS_NORMALIZE_VARIANTS`module-level nf-test、
`BWA_MEM2_MEM_SORT`module-level nf-test(Issue #8)、`VALIDATE_REFERENCE_CONTIGS`・
`GATK_GENOMICSDBIMPORT`・`GATK_GATHERVCFS`のmodule-level nf-test(Issue #11)を、
mainへのすべてのpush/pull requestに対してネイティブの`docker` profileで実行します
(`.github/workflows/test.yml`)。`.github/workflows/lint.yml`はリポジトリ構造のみを個別に
検査します。このnf-testスイートはIssue #12で検証したJoint Genotyping fixture契約、Issue #16で
検証したvariant QC出力契約、Issue #15で検証したFILTER/annotationカバレッジ/type会計契約、
Issue #8で検証したmapping/重複マークhardening契約、Issue #11で検証したJoint Genotyping
scale hardening契約にスコープされています。全モジュールを網羅するより広範なnf-testとfunctional
CIの見直しは[Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1) #Gで
引き続き追跡しています。

### Variant QC summarizer

`bin/summarize_variant_qc.py`は1つの`bcftools stats`レポートを、下記の`variant_qc.tsv`、
`sample_qc.tsv`、`summary.txt`へ変換します。これは`modules/local/bcftools_stats.nf`に
組み込まれていた約230行のAWKスクリプトを置き換えたもので(Issue #16)、集約ロジックを
Nextflow・Groovy・shellエスケープから独立して単体テストできるようにするためです。

pin済みの`bcftools:1.24` containerにはPython interpreterが含まれておらず、GATK containerに
同梱の`bcftools`を再利用するとpin済みのbcftoolsバージョンが1.24から同梱の1.13へ黙って
ダウングレードしてしまいます。どちらも受け入れる代わりに、variant QCは今や2つのNextflow
processとして実行されます: `BCFTOOLS_STATS`(変更なし、`bcftools stats`のみ実行)に続いて
`SUMMARIZE_VARIANT_QC`(専用のdigest-pin済み`python:3.12` container内で`bin/summarize_variant_qc.py`
を実行)。7つのstage/type組み合わせと28個の公開QC成果物は変更されていませんが、このQCステップは
Nextflow traceにおいて7ではなく14のタスク実行としてカウントされるようになりました。フル
(非`-slim`)のPython imageを使用しているのは、traceが有効な実行(nf-testと`-with-trace`は
いずれもこれに該当)ではNextflowがresource使用量metricsを収集するためtaskコンテナ内で`ps`を
shell outするためです。`python:3.12-slim`にはこれが含まれていません。

summarizerのunit testは標準ライブラリのみで、追加の依存なしに実行できます。

```bash
python3 -m unittest discover -s tests/bin -v
```

このrefactorの出力は、変更前後でクリーンなsynthetic Docker smoke testを実行し、
`qc/variants/`配下の全ファイルをdiffすることで以前のAWK実装と照合済みです: 7つのstage/type
組み合わせすべてにわたる28個の成果物すべてがbyte-identicalでした。

### Variant出力

| 出力 | 説明 |
| --- | --- |
| `variants/gvcf/<sample_id>.g.vcf.gz` | HaplotypeCallerが生成するサンプル単位のGVCF |
| `variants/gvcf/<sample_id>.g.vcf.gz.tbi` | サンプルGVCFのTabix index |
| `variants/raw/cohort.raw.vcf.gz` | 参照ゲノムのcontig順でgatherされた未フィルタの複数サンプルcohort VCF |
| `variants/raw/cohort.raw.vcf.gz.tbi` | raw cohort VCFのTabix index |
| `variants/by_type/cohort.snp.vcf.gz` | SNPとして選択されたrawコホートレコード |
| `variants/by_type/cohort.indel.vcf.gz` | indelとして選択されたrawコホートレコード |
| `variants/filtered/cohort.<snp-or-indel>.filtered.vcf.gz` | ハードフィルタラベルを適用したvariant type別VCF |
| `variants/pass/cohort.<snp-or-indel>.pass.vcf.gz` | 設定済みハードフィルタすべてをFILTER値がPASSしたレコード |
| `qc/variants/cohort.<stage>.<type>.bcftools.stats.tsv` | 完全な`bcftools stats`出力 |
| `qc/variants/cohort.<stage>.<type>.variant_qc.tsv` | 機械可読なcohortおよびvariant単位のQC指標 |
| `qc/variants/cohort.<stage>.<type>.sample_qc.tsv` | 機械可読なサンプル単位の遺伝子型・欠損率指標 |
| `qc/variants/cohort.<stage>.<type>.summary.txt` | 人間可読なQCサマリー |
| `qc/variants/cohort.filtered.<type>.filter_breakdown.tsv` | filtered VCFのFILTER値会計: 全体/PASS/非PASS、組み合わせ別・タグ別レコード数、複数タグレコード数 |
| `qc/variants/cohort.filtered.<type>.annotation_qc.tsv` | `QD`、`QUAL`、`SOR`、`FS`、`MQ`、`MQRankSum`、`ReadPosRankSum`のannotation別存在・フィルタヒット会計 |
| `qc/variants/cohort.filtered.<type>.filter_qc.summary.txt` | 人間可読なFILTER・annotationカバレッジサマリー |
| `qc/variants/cohort.variant_type_accounting.tsv` | `raw/all`を`raw/snp`・`raw/indel`選択に対してreconcileしたcohort全体の会計 |
| `qc/variants/cohort.variant_type_accounting.summary.txt` | 人間可読なvariant type会計サマリー |
| `variants/gs_normalized/cohort_gs.normalized.vcf.gz` | `raw/all`をleft-normalize・multiallelic分割したもの(`bcftools norm -m-`)。すべてのvariant typeをまだ含む |
| `qc/variants/cohort_gs.classification_accounting.tsv` | GS lineageのsplit後REF/ALT形状分類と重複key会計 |
| `qc/variants/cohort_gs.classification_accounting.summary.txt` | 人間可読な分類会計サマリー |
| `variants/gs_classified/cohort_gs.classified.vcf.gz` | 正規化VCFのbiallelic SNPのみの部分集合、FILTERを`.`へリセット |
| `variants/gs_filtered/cohort_gs.snp.filtered.vcf.gz` | 同じSNPハードフィルタを再適用したclassified VCF |
| `variants/gs_pass/cohort_gs.snp.pass.vcf.gz` | GS対象のPASSレコード。GSパネル自体のsource VCF |

上記の圧縮VCFにはすべて`.tbi` indexが付随します。Variant QCは7つのstage/type組み合わせ
(`raw/all`、`raw/snp`、`raw/indel`、`filtered/snp`、`filtered/indel`、`pass/snp`、
`pass/indel`)について生成されます。上記のFILTER/annotationカバレッジファイルは
variant typeごとに`filtered` VCFからのみ生成されます(`raw`のFILTER値は常に未設定、`pass`の
FILTER値は構造上常に`PASS`であるため、どちらのstageも意味のあるFILTER内訳を持ちません)。
variant type会計ファイルはcohortごとに1回生成されます。

cohort QCテーブルはサンプル数・レコード数・variant type数・multiallelic site数、
transition/transversion/Ti-Tv比、cohortの遺伝子型欠損数・欠損率、サンプル一覧を報告します。
サンプル単位テーブルはreference-homozygous・非reference-homozygous・heterozygous・欠損
遺伝子型数、欠損率、平均depth、singleton数を報告します。

raw/filtered/PASSのVCFは中間的な科学的結果のままです。PASSラベルの存在は、そのレコードが
設定済みのルールを通過したことのみを意味し、解析可能なSNPパネルであるとか、しきい値が
生物学的コホートに適していることの証拠であると解釈してはいけません。`VariantFiltration`は
レコードを一切削除しません: `filtered` VCFはそのraw stage・type選択後の入力と全く同じレコードを
保持しており、それぞれが除外されるのではなくFILTERラベルを付与されます。下流のステップ
(`PASS`抽出、上記のFILTER/annotationカバレッジQC)はそのより小さいレコード集合ではなく
このラベルを読み取ります。

### FilterとAnnotationカバレッジQC

レコードの`FILTER`列が特定のタグを欠いている理由は、全く異なる2つの場合があり得ます:
該当annotationが評価されてしきい値を満たした場合と、そのannotationがそのレコードに
**存在せず**フィルタ式が全く評価されなかった場合です。特に`MQRankSum`と`ReadPosRankSum`は、
比較可能なreference支持リード集団とalternate支持リード集団が存在しない場合にGATKによって
省略されます(たとえばこのリポジトリの現在のsynthetic fixtureのすべてのレコードは、両
annotationとも欠損しています)。`MQRankSum`/`ReadPosRankSum`の欠損は、対応する
`SNP_MQRANKSUM_LOW`/`SNP_READPOSRANKSUM_LOW`(または`INDEL_READPOSRANKSUM_LOW`)フィルタが
満たされたことを意味せず、単に評価できなかったことを意味します。`annotation_qc.tsv`の
`evaluable_rate`(存在数/全体数)と`filter_hit_rate`(タグ付き数/存在数、レコード0件の場合は
`0`ではなく`NA`として報告)は、まさにこの2つの状況を区別可能にするために存在します。
`evaluable_rate`も確認せずに、低いまたは0の`filter_hit_rate`をしきい値が適切に較正されている
証拠として読んではいけません。

`QUAL`はVCF自身の固定QUAL列であり、INFO annotationではなく、その大きさはコホートサイズと
depthにスケールします: 固定の`snp_filter_qual_min`/`indel_filter_qual_min`しきい値は、
サイズやシーケンシングdepthが異なるコホート間で同じ判別的意味を持ちません。
`annotation_qc.tsv`は一貫性のためINFO由来のannotationと同じ存在/欠損の慣例でこれを追跡
しますが、同じ方法で計算されているからではありません。

すべてのannotationがすべてのvariant typeにハードフィルタを持つわけではありません: indelには
`SOR`、`MQ`、`MQRankSum`フィルタがありません(下記パラメータ表参照)。そのため
`annotation_qc.tsv`は`variant_type=indel`のこれら3行について`filter_tag`、
`filter_tagged_records`、`filter_hit_rate`を`NA`として報告します。これはそのvariant typeに
対してそのannotationのフィルタリングがスコープ外であることを意味し、「欠損」
(`present_records`/`missing_records`で別途報告)や「評価されたがヒット0件」とは異なる状況です。

`GATK SelectVariants`の`--select-type-to-include`は各レコード全体を1つの overall type
(SNP、INDEL、MIXED、MNP、...)へ分類し、その単一typeとの完全一致でのみ選択します
(pin済みGATK 4.6.2.0 containerに対して`tests/modules/gatk_selectvariants.nf.test`で確認済み)。
そのためMIXED型のレコード(SNP型とindel型の両方のALT alleleを持つ1つのsite)は`SNP`にも
`INDEL`にもマッチせず、`cohort.snp.vcf.gz`と`cohort.indel.vcf.gz`の**両方**から除外されます。
両方に選択されることはありません。`records_not_selected`
(`variant_type_accounting.tsv`内の`raw_all_records - raw_snp_records - raw_indel_records`)は
通常運用では`>= 0`であることが期待されます。これは`raw/all`に存在するがtype別選択の両方から
除外されたMIXED/MNP/symbolic/その他のtypeのレコードをカウントするもので、`raw/snp`と
`raw/indel`が`raw/all`をpartitionすることは保証されていません。すべてのALT alleleが同じ
elementary type(たとえば2つのSNP allele)であるmultiallelic siteは、純粋な`SNP`または
`INDEL`として分類され、`records_not_selected`には含まれない点に注意してください。MIXED型の
siteのみが両選択から除外されるため、`number_of_multiallelic_sites`と`records_not_selected`は
正確に一致することが期待されているわけではありません。`raw/snp`と`raw/indel`が実際には
disjointでない場合、`records_not_selected`は負になり得ます。`variant_type_accounting.tsv`は
計算された値を常にそのまま報告し(隠したりclampしたりしない)、これをwarningとしてフラグ
します。あわせて`snp_indel_duplicate_records`(同一の`CHROM`/`POS`/`REF`/`ALT`が両方の選択に
存在するレコード)を、この不変条件が実際に破られ得る唯一の方法である、真に重複した出力
レコードの直接的な診断として報告します。

### ゲノミックセレクション(GS)パネル

このワークフローは`raw/all`から独立してゲノミックセレクション(GS) SNPパネルを導出します。
multiallelicおよびMIXED型のレコードを正規化・再分類してから、独立したlineage(`cohort_gs.*`)で
SNPハードフィルタを再適用するため、上記の主要な`raw`/`filtered`/`pass`出力はこのプロセスに
一切触れられません。dosage matrixスキーマ(v1)は設計上diploid専用であり、`sample_ploidy`が
`2`でない場合は(そうでなければすべての遺伝子型コールが黙って欠損として符号化されてしまうため)
即座に失敗します(出力なし)。完全なデータ契約 -- 入力/正規化/分類/対象選定ルール、dosage
符号化と`int8`が却下された理由、欠損/非標準遺伝子型ポリシー、ディスク上/メモリ内matrix形状、
空panel契約 -- は[`docs/gs_panel_data_contract.md`](docs/gs_panel_data_contract.md)に
文書化されており、ここでは出力一覧のみを繰り返します。

GS lineageはデフォルトで実行されます(`enable_gs_panel = true`)が、`--enable_gs_panel false`
で完全に無効化できます(Issue #20)。これは特に非diploidのvariant callingをend-to-endで実行
できるようにするために存在します: `sample_ploidy != 2`と(デフォルトの)`enable_gs_panel = true`
の組み合わせは、variant calling lineage全体を実行してGSパネル自身のdiploid専用チェックに
到達してから失敗する代わりに、どのprocess開始前にも即座に失敗します。`enable_gs_panel = false`
では、`variants/gs_*`や`gs_panel/`出力は一切生成されません。主要な`raw`/`filtered`/`pass`/QC
lineageはどちらの場合も影響を受けません。非diploidの実行が成功したこと自体は、ハード
フィルタしきい値がそのploidyに対して生物学的に適切であることの証拠ではありません --
[Variant callingおよびfilteringパラメータ](#variant-callingおよびfilteringパラメータ)参照。

| 出力 | 説明 |
| --- | --- |
| `gs_panel/cohort.gs_panel.genotype_matrix.tsv.gz` | Dosage matrix(variant行 x sample列)。`-1`/`0`/`+1`/`nan` |
| `gs_panel/cohort.gs_panel.sample_metadata.tsv` | サンプル単位の欠損/非標準遺伝子型数と率 |
| `gs_panel/cohort.gs_panel.variant_metadata.tsv` | Variant単位の`CHROM`/`POS`/`REF`/`ALT`/`QUAL`と欠損率 |
| `gs_panel/cohort.gs_panel.genotype_encoding_accounting.tsv` | Cohort全体の遺伝子型分類件数(標準dosage対欠損対非diploid対非biallelic-index、および独立したphased-call件数 -- phasingそれ自体はdosageに影響しない) |
| `gs_panel/cohort.gs_panel.genotype_encoding_accounting.summary.txt` | 人間可読な遺伝子型符号化サマリー |
| `gs_panel/cohort.gs_panel.record_accounting.tsv` | full-lineageのレコードreconciliation(`raw_all` -> `normalized` -> `classified` -> `gs_pass` -> matrix/metadata、matrixファイル自体との相互チェック)と`panel_status`(`empty`/`populated`)。不一致はwarningではなくhard errorであり、その場合ファイルは一切書き込まれない |
| `gs_panel/cohort.gs_panel.record_accounting.summary.txt` | 人間可読なレコード会計サマリー |
| `gs_panel/cohort.gs_panel.manifest.json` | スキーマversion付きの再現性manifest: run ID、container digest、`sample_ploidy`と遺伝子型符号化スキーマ、そしてパネル成果物すべてとそれを構築したraw/all VCF・参照FASTA/FAIのchecksum |

0件の対象レコードを持つGSパネルは正常な結果であり、エラーではありません -- 実際、このパイプ
ラインのデフォルトsynthetic fixtureが現在生成しているものです。両方のsynthetic SNPが
`SNP_SOR_HIGH`で失敗するためです(主要lineageの`pass/snp`出力と同じ結果)。matrix header
(サンプル一覧)とサンプルmetadataは完全に書き込まれます。欠けているのはvariant行のみで、
`record_accounting.tsv`の`panel_status`フィールドが明示的にそう述べます。

本文書執筆時点で、[`genomic-prediction-resnet-hybrid`](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid)
にはアズキ固有・VCF固有のingestionコードはありません。このパネルの規約(dosage符号化、dtype、
ディスク上の形状、manifest構造)は、その規約に合わせて構築できるよう、そのリポジトリ自身の
テスト済みSoyNAMロード規約に合わせて選択されたものであり、対応するloaderが既に存在するから
ではありません。

### Samplesheet契約

入力はCSVファイルである必要があります。

| 列 | 必須 | 説明 |
| --- | --- | --- |
| `sample_id` | Yes | 生物学的サンプル識別子。複数read groupのために値の重複が許容される |
| `read_group_id` | Yes | samplesheet全体で一意でなければならない識別子 |
| `fastq_1` | Yes | `.fq.gz`または`.fastq.gz`で終わる実在のread 1ファイル |
| `fastq_2` | Yes | `.fq.gz`または`.fastq.gz`で終わる実在のread 2ファイル |
| `library_id` | Yes | シーケンシングlibrary識別子 |
| `platform` | Yes | シーケンシングplatform。初期契約は`ILLUMINA`を受け付ける |
| `platform_unit` | No | Flowcell・lane・sample barcode識別子 |

`platform_unit`を指定する場合、`FLOWCELL.LANE.SAMPLE_BARCODE`のような値でread groupを区別
すべきです。異なるread group間で再利用すべきではありません。

例:

```csv
sample_id,read_group_id,fastq_1,fastq_2,library_id,platform,platform_unit
sample_a,sample_a_L001,reads/a_L001_R1.fastq.gz,reads/a_L001_R2.fastq.gz,lib_a,ILLUMINA,flowcell1.L001.ATCACG
sample_a,sample_a_L002,reads/a_L002_R1.fastq.gz,reads/a_L002_R2.fastq.gz,lib_a,ILLUMINA,flowcell1.L002.ATCACG
```

想定外の列、重複する`read_group_id`値、存在しないファイル、read 1/read 2の同一パス、複数
read groupで再利用されたFASTQファイルは、解析processが開始する前に拒否されます。

### 参照ゲノムBundle契約

以下のパラメータが参照ゲノムbundleを定義します。

| パラメータ | 必須 | 説明 |
| --- | --- | --- |
| `reference_id` | Yes | 参照ゲノムbundleの安定した識別子 |
| `reference_name` | Yes | 人間可読なassembly名 |
| `reference_fasta` | Yes | 実在の非圧縮参照FASTA |
| `reference_accession` | No | 公開データベースaccession |
| `reference_species` | No | 参照ゲノムが表す種 |
| `reference_cultivar` | No | 参照ゲノムが表すcultivar |
| `reference_fai` | No | `<reference_fasta>.fai`と命名された互換prebuilt FASTA index |
| `reference_dict` | No | `<reference_basename>.dict`と命名された互換prebuilt sequence dictionary |
| `bwa_index_prefix` | No | basenameが`reference_fasta`と一致する互換prebuilt BWA-MEM2 index prefix |

synthetic test設定は`conf/test.config`で定義されています。Longxiaodou 4の例は
`conf/references/longxiaodou4.config.example`にあります。

オプションのindexパラメータを省略すると、ワークフローはSAMtools・GATK・BWA-MEM2で対応する
indexを生成します。指定した場合、indexパスと期待されるファイル名はprocess実行前に検証
されます。

Issue #11で、生成/prebuilt両方の参照ゲノムパスに対して、`.fai`と`.dict`が実際に同じcontigを
同じ順序・同じ長さで記述していることを検証するゲートが追加されました。詳細は
[Joint Genotyping scale hardening](#joint-genotyping-scale-hardeningissue-11)を参照してください。

### Variant-callingおよびFilteringパラメータ

| パラメータ | デフォルト | 説明 |
| --- | ---: | --- |
| `sample_ploidy` | `2` | 正の整数: run全体の**global**なploidyで、HaplotypeCallerとGenotypeGVCFsの両方に渡される。すべてのサンプル・すべてのcontigに一様に適用され、このパイプラインは1回のrun内でのmixed ploidyをサポートしない。[GSパネル](#ゲノミックセレクションgsパネル)のdosage matrixスキーマ(v1)はdiploid専用であり、これが`2`でなく`enable_gs_panel`が`true`(デフォルト)の場合は即座に失敗する -- 下記参照 |
| `enable_gs_panel` | `true` | GS panel lineageを実行するかどうか。`sample_ploidy`が`2`以外の場合は`false`にする必要がある(Issue #20) |
| `optical_duplicate_pixel_distance` | `100` | 非負整数: GATK MarkDuplicatesの`--OPTICAL_DUPLICATE_PIXEL_DISTANCE`へ渡す、2つの重複がoptical(PCR/library重複ではなく)重複と呼ばれるための、readクラスタ間の最大pixel距離。`100`はGATK自身の未変更デフォルトと一致し、unpatterned flowcell(例: HiSeq 2500)に適している。Patterned flowcell(例: HiSeq X、NovaSeq)はクラスタをはるかに密に詰め込むため、GATK/Picardのガイダンスに従いより大きな値(しばしば約2500)を必要とすることが多い。これはlibrary size推定に用いるoptical対PCR重複の*下位分類*にのみ影響し、座標/CIGARベースの基盤となる重複FLAG自体には影響しない(Issue #8) |
| `genomicsdb_batch_size` | `50` | 正の整数: GATK GenomicsDBImportの`--batch-size`へ渡す。各contig単位workspaceを構築する際に同時に開くサンプル数を制御する。`50`は**現時点の初期運用値**であり、このパイプラインの目標である327検体コホートや20〜30検体規模での最適値として検証された値ではない(Issue #11)。5検体規模での実データ動作は確認済み(Issue #26) -- [`docs/joint_genotyping_scaling.md`](docs/joint_genotyping_scaling.md)参照 |
| `snp_filter_qd_min` | `2.0` | `QD`がこの値未満のSNPを`SNP_QD_LOW`としてマークする |
| `snp_filter_qual_min` | `30.0` | `QUAL`がこの値未満のSNPを`SNP_QUAL_LOW`としてマークする |
| `snp_filter_sor_max` | `3.0` | `SOR`がこの値を超えるSNPを`SNP_SOR_HIGH`としてマークする |
| `snp_filter_fs_max` | `60.0` | `FS`がこの値を超えるSNPを`SNP_FS_HIGH`としてマークする |
| `snp_filter_mq_min` | `40.0` | `MQ`がこの値未満のSNPを`SNP_MQ_LOW`としてマークする |
| `snp_filter_mq_rank_sum_min` | `-12.5` | `MQRankSum`がこの値未満のSNPを`SNP_MQRANKSUM_LOW`としてマークする |
| `snp_filter_read_pos_rank_sum_min` | `-8.0` | `ReadPosRankSum`がこの値未満のSNPを`SNP_READPOSRANKSUM_LOW`としてマークする |
| `indel_filter_qd_min` | `2.0` | `QD`がこの値未満のindelを`INDEL_QD_LOW`としてマークする |
| `indel_filter_qual_min` | `30.0` | `QUAL`がこの値未満のindelを`INDEL_QUAL_LOW`としてマークする |
| `indel_filter_fs_max` | `200.0` | `FS`がこの値を超えるindelを`INDEL_FS_HIGH`としてマークする |
| `indel_filter_read_pos_rank_sum_min` | `-20.0` | `ReadPosRankSum`がこの値未満のindelを`INDEL_READPOSRANKSUM_LOW`としてマークする |

Ploidyは1回のpipeline runにおいて、すべての生物学的サンプル・すべてのcontigに一貫して適用
されます。mixed-ploidyコホートは現在サポートされておらず、サンプル単位・contig単位のploidy
上書きもありません。

`--sample-ploidy`はHaplotypeCallerとGenotypeGVCFsの両方に明示的に渡されます(Issue #20)。
これはバグ修正ではなく契約強化です: pin済みGATK 4.6.2.0 containerに対して、GenotypeGVCFsから
`--sample-ploidy`を省略した場合と明示的に渡した場合(実際のhaploid gVCFから構築した同一の
GenomicsDBワークスペースに対して`--sample-ploidy 1`で検証)は、比較したすべてのfield --
`CHROM`/`POS`/`REF`/`ALT`、サンプル順、`GT`、`AC`、`AN`、`AF`、`PL` -- でbyte-identicalな
出力を生成しました。GenotypeGVCFsは自身のGenomicsDB/gVCF入力から正しいploidyを既に自力で
導出します。それでも`--sample-ploidy`を明示的に渡すのは、GATK自身の(ここでは正しいが
私たちが強制していない)推論挙動 -- 将来のGATKバージョンで警告なく変わり得る -- に依存する
のではなく、コマンドライン自体からパイプラインの意図を監査可能にするためです。

非diploidのrunが成功したこと -- このリポジトリのhaploid end-to-end nf-testに合格することを
含む -- は、ploidyが正しく伝播しパイプラインの機構が動作することを検証するものです。上記の
SNP/indelハードフィルタしきい値が非diploidサンプルに対して生物学的に適切であることの証拠
ではありません。これらのしきい値はこのリポジトリのdiploid synthetic fixtureに対してのみ
検証されてきました。

フィルタリングのデフォルト値は設定可能な運用上の出発点です。QUAL・SOR・indel固有のフィルタで
歴史的pilotのルールを拡張していますが、参照ゲノム・検体・シーケンシングdepth・コホートサイズに
わたる妥当性は検証されていません。しきい値を変更することは科学的出力を変更することを意味し、
run parameterと共に記録すべきです。

### ワークフローの構成

```text
main.nf
nextflow.config
nextflow_schema.json
nf-test.config
assets/schema_input.json
conf/test.config
conf/references/
bin/
workflows/
subworkflows/local/
modules/local/
tests/data/
tests/scripts/
tests/pipeline/
tests/modules/
tests/bin/
tests/nextflow.config
```

---

## Joint Genotyping Scale Hardening(Issue #11)

Issue #8 Phase 5がSeedcore-01実機に残していた、実際のLongxiaodou 4参照ゲノムとSRR29909135
検体の単一accession runは、`GATK_GENOMICSDBIMPORT`で失敗していました。根本原因は、Nextflowの
`path` input qualifierが、そのchannel slotに対して解決されるファイルがちょうど1つの場合に
single-element Listを裸のscalar `Path`/`File`へ暗黙的にunwrapするという、文書化された
Nextflowの挙動でした(このパイプライン固有のバグではありません)。表示されたエラー
(`"1176577035 gVCFs, 642455 indexes"`)は、実際にはgVCFファイルとその`.tbi` indexファイルの
byte単位のサイズそのものでした(`.size()`をscalar `File`に対して呼び出すとbyte長が返り、
"1"は返らないため)。このリポジトリのすべてのsynthetic fixtureは2サンプル以上を使用するため、
実際の単一accession runがこれに遭遇するまで、この問題は表面化していませんでした。

この発見をもとに、以下を実施しました。

- **List-safety修正**: `GATK_GENOMICSDBIMPORT`と、同じ実データ証跡から発見された同一パターンを
  持つ`GATK_GATHERVCFS`の両方で、`gvcfs`/`gvcf_indexes`(および`vcfs`/`vcf_indexes`)を
  `.size()`や`.collect{}`を呼び出す前にListへ正規化するようにした。
- **メモリ契約**: `GenomicsDBImport`のJVM heapを、`task.memory`の安全な80%(GiB丸めではなく
  MiB単位で計算)へ上限設定した。以前の固定1GiB予約は、既存の16GiB割り当てにおいて
  GenomicsDBImportのnative TileDBストレージ層にわずか6.25%の余裕しか残していなかった。
  必要な20%以上のnative層余裕はこの80%上限から自動的に導かれるものであり、別途assertされた
  数値ではない。専用の`process_genomicsdb`resourceラベルを追加し、値は`process_high`から
  変更せず引き継いだ(`cpus=8`、`memory = 16.GB * task.attempt`、`time='24h'`)。
  `task.attempt`ベースのOOM retryスケーリングは維持されており、80%比率はすべてのretry
  attemptで`task.memory`から新たに計算し直される。
- **`genomicsdb_batch_size`パラメータ**: 新規追加(デフォルト`50`、最小値`1`)。詳細は
  上記の[Variant-callingおよびFilteringパラメータ](#variant-callingおよびfilteringパラメータ)
  を参照。
- **参照ゲノム契約**: `bin/validate_reference_contigs.py`が、`.fai`と`.dict`のcontig名・長さ・
  **順序**を位置ベース(setベースではなく)で比較する。同一contig集合だが順序が異なる組は
  naive setベースの比較では見逃されるが、これは明示的にrequiredかつ合格しているnegative test
  ケースであり、module levelとpipeline levelの両方でカバーされている。`modules/local/validate_reference_contigs.nf`
  がpass-throughゲートとして`workflows/adzuki_snp_pipeline.nf`に組み込まれており、生成・
  prebuilt両方の参照ゲノムパスが収束した直後に位置する。不一致はどのGATK processが開始
  する前にも、不一致の最初の位置を示すactionableなメッセージで実行全体を失敗させる。

`--consolidate`は意図的に有効化しておらず、新たな`genomicsdb_consolidate`パラメータも
追加していません。再評価すべき具体的な条件を含む完全な文書は
[`docs/joint_genotyping_scaling.md`](docs/joint_genotyping_scaling.md)にあります。

実データによるターゲットsmoke test(Issue #8 Phase 5から得た同じ実gVCFを再利用し、
FASTQ→mapping→HaplotypeCallerの再実行は行わなかった)では、修正後のコードが実際の
染色体規模contig(`NC_068970.1`、65.4 Mb)と小さなunplaced scaffold(`NW_026294847.1`、
2,353 bp)の両方で成功することを確認しました。これは1検体分の証跡であり、複数サンプルや
コホート規模での検証ではありません -- Issue #26がこれを5検体規模まで拡張しています
(下記参照)。

---

## 実データコホートE2E検証(Issue #26)

<!-- 実行結果は進行中: Seedcore-01での実行完了後にこのセクションを更新する -->

Issue #8(mapping hardening)とIssue #11(Joint Genotyping scale hardening)が完了した後、
Issue #26はSeedcore-01実機上で、公開WGSデータを用いた複数サンプルのSNPコーリングを段階的に
検証し、受託解析で必要となるVCF・QC・遺伝子型行列・処理時間・メモリ使用量・再現性証跡を
確認するものです。327検体のフル実行はこのIssueのスコープ外であり、3〜5検体規模の検証結果を
もとに20〜30検体への拡張条件を評価します。

このセクションの実行結果は現在準備中です。最新の内容は
[`docs/real_cohort_e2e.md`](docs/real_cohort_e2e.md)を参照してください。

---

## Pilot Workflow

![Manual pilot workflow](docs/workflow.png)

以下のコマンドは、歴史的な単一サンプルpilotを文書化したものです。実行可能なNextflow実装の
技術的証跡・歴史的背景として保持しています。

依存関係のバージョン、checksum、resource要件、すべての中間検証ステップが固定されていないため、
これらは**クリーンな環境での再現手順ではまだありません**。

可読性とshellの正しさのため、記録されたコマンドには出力ディレクトリ作成の追加、参照パスの
quote、`--native-pair-hmm-threads`の長いoption名の綴り修正といった軽微な編集を加えています。
これらの編集は、この文書更新のために完全な手順を再実行したことを意味しません。

---

## 歴史的な環境構築

```bash
conda create -n bioinfo -c conda-forge -c bioconda \
  fastqc sra-tools bwa samtools gatk4 bcftools ncbi-datasets-cli -y
conda activate bioinfo
```

---

## 歴史的な単一サンプル手順

### 1. 参照ゲノムのダウンロード

```bash
datasets download genome accession GCF_016808095.1 \
  --include genome \
  --filename adzuki_reference.zip

unzip adzuki_reference.zip -d adzuki_reference

REF=adzuki_reference/ncbi_dataset/data/GCF_016808095.1/GCF_016808095.1_ASM1680809v1_genomic.fna
```

### 2. 参照ゲノムのIndex作成

```bash
bwa index "$REF"
gatk CreateSequenceDictionary -R "$REF"
samtools faidx "$REF"
```

### 3. デモンストレーション用WGSデータのダウンロード

```bash
prefetch SRR29909135 --output-directory ./raw_data

fasterq-dump ./raw_data/SRR29909135/SRR29909135.sra \
  --outdir ./raw_data \
  --threads 4 \
  --progress
```

### 4. 初期品質管理

```bash
mkdir -p fastqc_results

fastqc \
  ./raw_data/SRR29909135_1.fastq \
  ./raw_data/SRR29909135_2.fastq \
  --outdir ./fastqc_results \
  --threads 4
```

この歴史的手順はFastQCによる検査を行いますが、readのtrimmingは行いません。実行可能な
Nextflowワークフローは現在、trimming前後のFastQCとペアエンドfastp trimmingを実行します。

### 5. Mapping

```bash
bwa mem -t 4 \
  -R "@RG\tID:SRR29909135\tSM:SRR29909135\tPL:ILLUMINA" \
  "$REF" \
  ./raw_data/SRR29909135_1.fastq \
  ./raw_data/SRR29909135_2.fastq \
  | samtools sort -@ 4 -o SRR29909135.bam

samtools index SRR29909135.bam
```

### 6. 歴史的pilotでの重複除去

歴史的なコマンドは`samtools markdup -r -d 2500`を使用していました。`-r`オプションは重複
リードを単にフラグ設定するのではなく削除していました。この不可逆な挙動は、pilot手順の記録
としてのみここに保持しています。実行可能なNextflowワークフローは代わりに
`REMOVE_DUPLICATES=false`のGATK MarkDuplicatesを使用し、重複レコードをDUPフラグ付きで
保持することで、下流のツールが元の証跡を破壊せずにそれらを除外できるようにしています。

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

### 7. GVCFモードでのVariant Calling

```bash
gatk HaplotypeCaller \
  -R "$REF" \
  -I SRR29909135.markdup.bam \
  -O SRR29909135.g.vcf.gz \
  -ERC GVCF \
  --tmp-dir /tmp \
  --native-pair-hmm-threads 4
```

### 8. デモンストレーションサンプルのGenotyping

```bash
gatk GenotypeGVCFs \
  -R "$REF" \
  -V SRR29909135.g.vcf.gz \
  -O SRR29909135.vcf.gz
```

このコマンドは単一サンプルpilotでのみ使用されました。実行可能な複数サンプルpipelineは
HaplotypeCaller gVCFに続いてGenomicsDBImportとGenotypeGVCFsをJoint Genotypingに使用します。
Joint Genotypingは30検体しきい値に条件付けられるのではなく、複数サンプルコホートの
デフォルトです。

### 9. SNPの選択とフィルタリング

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

上記のハードフィルタしきい値は歴史的pilotを文書化したものです。参照ゲノム・検体・コホート
サイズにわたる妥当性はまだ検証されていません。

---

## 設計判断

### Base Quality Score Recalibration (BQSR)

BQSRは文書化された歴史的手順と実行可能なワークフローの両方から意図的に除外されています。
GATK BaseRecalibratorはknown-sites VCFを使用して、既知の多型とシーケンシングエラーモデル化に
使うミスマッチを区別します。このリポジトリは、Longxiaodou 4参照ゲノムbundleに対して適切な
known-sitesリソースを特定・検証していません。

Bootstrapped BQSRも、known-sites構築とそれが結果に与える影響がここでは検証されていないため、
現在のスコープ外です。これは実装漏れではなく意図的な設計判断です。
[GATK BaseRecalibratorの文書](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
を参照してください。

### 参照ゲノムBundleポリシー

参照アセンブリは明示的なpipeline入力です。異なるcultivarやassemblyに対して生成された結果を、
別途の検証なしに直接互換なものとして扱ってはいけません。Longxiaodou 4は文書化された実データの
例として保持し、test profileは再配布に安全な小さなsynthetic参照ゲノムを使用します。

---

## 単一サンプルPilotの結果

![Mapping statistics from the single-sample pilot](docs/mapping_stats.png)

| ステップ | 件数 | 備考 |
| --- | ---: | --- |
| Total reads | 57,597,756 | Paired-end, 150 bp |
| Mapping rate | 99.35% | GCF_016808095.1に対するBWA-MEM |
| Properly paired | 93.59% | 単一サンプルpilot |
| Duplicate rate | 約10% | 歴史的run |
| Total variants | 783,836 | SNPとindel |
| SNPs extracted | 678,212 | SelectVariants後 |
| PASS SNPs | 610,790 | 文書化されたハードフィルタ後 |

![Variant counts from the single-sample pilot](docs/variant_counts.png)

![SNP summary from the single-sample pilot](docs/snp_summary.png)

これらの値はSRR29909135に対する1回の歴史的実行を記述したものです。327検体のフルコホートに
対する推定値ではなく、他のワークフローに対するSNPコーリング精度や優位性を実証するもの
でもありません。

---

## 実装ロードマップ

実装ロードマップは[Issue #1](https://github.com/hoso-jpn/adzuki-snp-pipeline/issues/1)で
追跡しています。

以下の基盤・解析・QC機能を実装済みです。

- strict構文パーサ対応のNextflow DSL2ワークフロー
- パラメータとsamplesheetの検証
- index生成/再利用に対応した設定可能な参照ゲノムbundle
- raw/trimmed FastQCとペアエンドfastp trimming
- read-groupを意識したBWA-MEM2マッピングと座標ソート
- サンプル単位のread-group merge
- 重複レコードを削除しないlibrary単位のGATK重複マーク
- BAM indexとSAMtools flagstat/stats/idxstatsレポート
- サンプル単位のGATK HaplotypeCaller GVCF生成
- contig単位のGenomicsDBImportとGenotypeGVCFs
- indexed raw cohort VCFへのreference順gathering
- indexedなvariant type別出力を伴うSNP/indel分離
- 設定可能なGATKハードフィルタリングとindexedなPASSのみの出力
- raw/filtered/PASSの`bcftools stats` QC
- 機械可読なcohort・サンプル単位QCテーブル
- 人間可読なvariant-QCサマリー
- 再配布可能なdeterministic fixtureと機能的なDocker smoke-test profile
- 両サンプルが両方のdeterministic SNP locusをカバーし、欠損コールではなく確信度の高い
  `1/1`/`0/0`遺伝子型へ解決されるsynthetic fixture
- mainへのすべてのpush/pull requestでGitHub Actions上で実行される、その遺伝子型/annotation
  契約を検証するnf-test pipeline-level test
- filtered VCFのFILTER値・annotation別evaluable rate会計。annotation欠損としきい値未達を
  区別し、cohort全体の`raw/all`対`raw/snp`/`raw/indel`のreconciliationも行う
- 文書化されたGSパネルデータ契約(Issue #1 #F): `raw/all`の正規化とsplit後のmultiallelic・
  MIXED型レコード再分類、独立したGS固有ハードフィルタ再適用、sample/variant metadataを伴う
  dosage遺伝子型matrix、full-lineageレコード会計、software version・parameter・checksum
  情報を伴うスキーマversion付き再現性manifest。実際にpinされたbcftools containerに対して
  MIXEDレコード分割を検証するmodule-level nf-testも含む
- 実データ投入に先立つmapping工程のhardening(Issue #8): 中間`.sam`を作らない単一のpiped
  `BWA_MEM2_MEM_SORT`タスク、文書化された非oversubscribe CPU/メモリ計算式、`groupKey`
  経由のread-groupを意識した早期サンプル単位merge、`task.attempt`でスケールするOOM retry、
  設定可能な`optical_duplicate_pixel_distance`、GATKのread-name解析を実際に検証する
  CASAVA形式のsynthetic read名、そして上記すべての恒久的なnf-testカバレッジ。**Phase 5**
  (Issue #8)はさらに`BWA_MEM2_INDEX`/`BWA_MEM2_MEM_SORT`を実際のLongxiaodou 4参照ゲノムと
  実際の約19.3xカバレッジのWGS検体に対して実機(Seedcore-01)で測定し
  ([`docs/mapping_real_reference_profile.md`](docs/mapping_real_reference_profile.md)参照)、
  上記の80/20 CPU分割とメモリ分割計算式がその1つの参照ゲノム/検体/マシンに対して、OOMや
  retryなく安全であることを確認した。これは1台のマシン・1検体の測定であり、パフォーマンス
  保証やコホート規模の検証ではない(Issue #11/#26)
- 複数サンプルコホート投入に先立つJoint Genotyping工程のscale hardening(Issue #11):
  実際の単一accession run(Issue #8 Phase 5と同じSeedcore-01/Longxiaodou 4/SRR29909135の
  gVCF)が、正確に1つのサンプルが解決された場合にのみ表面化するNextflowの`path`-input
  List/scalar collapseバグを`GATK_GENOMICSDBIMPORT`(および同じ証跡から発見された同一
  パターンを持つ`GATK_GATHERVCFS`)で露呈させ、これを修正した。`GenomicsDBImport`のJVM
  heapは、既存の16GiB割り当てにおいてわずか6.25%の余裕しかnative TileDBストレージ層に
  残していなかった以前の固定1GiB予約の代わりに、`task.memory`の安全な80%(GiB丸めでは
  なくMiB単位で計算)へ上限設定されるようになった。新しい`genomicsdb_batch_size`パラメータ
  (現時点の初期運用値として文書化、検証済みの値ではない)、そしてどちらのGATK process開始
  前にも生成・prebuilt両方の参照ゲノムbundleパスをゲートするFAI/sequence dictionary
  contig名・長さ・**順序**検証(`bin/validate_reference_contigs.py`)を追加した。完全な
  契約、実際の単一accession観察結果(同じgVCFを再利用したターゲット実データ
  GenomicsDBImport smoke testを含む)、20〜30/327検体規模で未検証のまま残っている事項は
  [`docs/joint_genotyping_scaling.md`](docs/joint_genotyping_scaling.md)を参照
- Seedcore-01実機上での5検体実データE2Eコホート検証(Issue #26): 同一BioProjectに属する
  5件の公開WGS検体を、QC→mapping→重複マーク→gVCF→Joint Genotyping→hard filtering→
  variant QC→GSパネルまで通した。実行時間・peak RSS・storage使用量・sample/variant
  accounting・再現性証跡を実測で記録した。20〜30検体・327検体規模への拡張条件を、この
  実測に基づいて評価した(実行を伴わない)。詳細は
  [実データコホートE2E検証](#実データコホートe2e検証issue-26)と
  [`docs/real_cohort_e2e.md`](docs/real_cohort_e2e.md)を参照

以下の解析・再現性機能は計画中です。

- MultiQCレポート集約
- 全モジュールを網羅するnf-testカバレッジと、Pythonでテスト可能なvariant-QCモジュール
  (Issue #1 #G)
- `genomic-prediction-resnet-hybrid`におけるアズキ/VCFパネルのingestion経路(そのリポジトリ側の
  follow-up issueとして追跡。ここではスコープ外)
- 20〜30検体・327検体規模でのコホート検証(Issue #26)

ある能力は、対応するコードと検証がこのリポジトリに実際に存在して初めて実装済みとみなされます。

---

## データ取扱いポリシー

この公開リポジトリに含まれるのは以下のみです。

- ワークフローとpipelineのソースコード
- 再現性のために必要なコマンドと設定
- 公開データセットに基づく解析
- syntheticまたは再配布可能なテストデータ
- 機密情報を含まない技術検証記録

以下は含まれません。

- 前職またはその他の非公開研究データ
- 顧客データ
- 専有のSNPパネルやマーカー選定
- 機密のビジネスロジック
- 認証情報・アクセストークン・非公開インフラの詳細
- 将来の知的財産として保護すべき資料

生のWGSデータと完全な参照ゲノムbundleはこのリポジトリにcommitされません。それぞれの
権威ある公開ソースから取得する必要があります。

---

## 関連リポジトリ

以下のリポジトリは、より長期的な研究スタックにおける異なる段階を表しています。まだ
自動化されたend-to-endワークフローでは接続されていません。

- [adzuki-snp-pipeline](https://github.com/hoso-jpn/adzuki-snp-pipeline) — 公開WGSデータから
  コホートvariantとSNP行列まで。現在再構築中
- [adzuki-gwas-analysis](https://github.com/hoso-jpn/adzuki-gwas-analysis) — 公開されている
  GWAS summary statisticsの解析
- [genomic-prediction-resnet-hybrid](https://github.com/hoso-jpn/genomic-prediction-resnet-hybrid) — 互換性のある
  個体レベルデータを用いたGBLUPとニューラルゲノム予測モデルの監査可能な比較

`adzuki-gwas-analysis`が現在使用している公開Dryadデータセットには個体レベルの遺伝子型・
表現型は含まれていません。したがって、これら3つのリポジトリを既に稼働しているSNP→GWAS→
ゲノム予測のend-to-endサービスとして読むべきではありません。

本文書執筆時点で、`genomic-prediction-resnet-hybrid`にはアズキ固有・VCF固有のingestionコード
はなく、検証済みの唯一のデータ経路は無関係なSoyNAMダイズデータを読み込むものです。この
pipelineの[GSパネル](#ゲノミックセレクションgsパネル)は、そのloaderが既に存在するからでは
なく、将来のloaderがそれに合わせて構築できるよう明確な契約を持てるように、そのリポジトリ自身の
テスト済み規約(dosage符号化、dtype、ディスク上の形状、manifest構造)に合わせて設計されました。

---

## 著者

**Yusuke Hosokawa**<br>
独立研究者・AIエンジニア<br>
長期的な農業AI構想[Florigen AI](https://florigen.ai)を構築中

Plant Genetics × Edge AI × Physical AI

- [GitHub](https://github.com/hoso-jpn)
- [researchmap](https://researchmap.jp/hosokawa-yusuke)
- [Breeding Science (2025) — イネ染色体segment substitution系統におけるQTLマッピング](https://doi.org/10.1270/jsbbs.24058)

このリポジトリのゲノミクス研究は、長期的な農業AI・ロボティクス研究の土台となる植物ドメイン・
バイオインフォマティクスの専門性を文書化するものです。

---

## License

このプロジェクトはMIT Licenseの下でライセンスされています。[LICENSE](LICENSE)を参照して
ください。
