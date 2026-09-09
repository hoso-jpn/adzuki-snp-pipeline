# Docker実行のhost user所有契約

Issue #51で導入した実行時契約です。scientific semantics、threshold、schema、
artifact naming、pin済みtool/container version、resource policyは変更していません。

## 何が問題だったか

Dockerのcontainerは既定でcontainer内rootとしてprocessを実行するため、taskが
work directoryへ書いたfileはすべて`root:root`所有で生成されていました。published
artifact自体は`publishDir`が`mode: 'copy'`でありNextflow（host process）がcopyする
ため元々host user所有ですが、work directoryと`.nf-test`作業データは起動ユーザー自身が
削除できず、失敗runやtest runのたびに回収不能なデータが残り続けます。

synthetic profileの1 runで、実際に208 entriesがroot所有で残っていました
（`-profile test,docker`、2 sample / 3 read group fixture）。

## 変更内容

`nextflow.config`の`docker.runOptions`へ`-u $(id -u):$(id -g)`を追加し、container
taskを起動したhost userとして実行します。UID/GIDはNextflowではなくtask launcherの
shellが評価するため、実行hostの実ユーザーになります。

### profile間の上書き関係

`docker.runOptions`は追記ではなく単一keyへの代入です。`docker_amd64`は
`--platform linux/amd64`のために同じkeyへ代入するため、`docker` profile側にだけ
`-u`を書くと`docker_amd64`選択時に無効化されます。

そこでbase assignmentを`profiles`の外へ置き、`docker_amd64`では同じoptionを明示的に
併記しています。Nextflow 26.04には`docker.platform`設定が存在せず、記述してもparseは
通るものの無視され`--platform`が一切出力されないため、`--platform`を`runOptions`から
外すことはできません。config構文は共有変数宣言（`def`）も
`Variable declarations cannot be mixed with config statements`として拒否します。
この二重記述は`tests/bin/test_docker_host_user_contracts.py`が固定しています。

実効値は次のとおりです（`nextflow config -profile ...`で確認）。

| profile | 実効`docker.runOptions` |
| --- | --- |
| `docker` | `-u $(id -u):$(id -g)` |
| `docker_amd64` | `-u $(id -u):$(id -g) --platform linux/amd64` |
| `docker,docker_amd64` | `-u $(id -u):$(id -g) --platform linux/amd64` |
| `docker_amd64,docker` | `-u $(id -u):$(id -g) --platform linux/amd64` |
| `test,docker` | `-u $(id -u):$(id -g)` |
| `test,docker_amd64` | `-u $(id -u):$(id -g) --platform linux/amd64` |

### `docker.fixOwnership`を採用しなかった理由

`docker.fixOwnership`はtask commandの後にcontainer内で`chown`を実行する仕組みで、
task正常終了時にしか適用されません。Issue #51が回収対象としている滞留データの主因は
むしろ失敗run・中断runのwork directoryであり、そこには適用されません。`-u`は
そもそもroot所有のfileを作らないため、失敗の有無に依存しません。

### MULTIQCのstaging root

MULTIQCはcontainer絶対path `/multiqc` 配下へinput symlinkを張り、そこだけをscanします。
これはIssue #38が、MultiQCがmultiqc_data.json / multiqc_sources.jsonへ記録するsource
pathをcontainer pathに閉じ、host work directoryの絶対pathを公開artifactへ漏らさない
ために選んだ設計です（[`multiqc.md`](multiqc.md)）。

host user実行では、非特権ユーザーがcontainerのfilesystem rootへ`/multiqc`を作成できず
`mkdir: cannot create directory '/multiqc': Permission denied`で失敗します。published
provenanceのsource path契約を変えないため、staging root pathは変更せず
`containerOptions '--tmpfs /multiqc:rw'`で書き込み可能なmountを与えています。

task work directoryを`/multiqc`へbindする案は、staging rootとtask directoryが同一に
なりscript内の`cp multiqc_config.yaml /multiqc/multiqc_config.yaml`がsame-file copy
errorになるため採用していません。tmpfsに置かれるのはMultiQCのstagingとreport outputで、
scriptがtask directoryへcopyし直します。この分はdiskではなくtaskのmemory budgetへ計上
されます。sizeはsynthetic fixture規模でしか観測していません。

## 検証

Issue #51の検証順序（MultiQC → FastQC/GATK → GenomicsDBImport）に沿って、synthetic
profileの`-profile test,docker`をmodified/unmodifiedで実行して比較しました。

- MULTIQC: 上記のとおり`-u`のみでは失敗し、tmpfs mount追加後に成功。
  `report_data_sources`は23 sourceすべてが`/multiqc/input/`配下、`report_saved_raw_data`
  のmodule keyもunmodified runと一致。multiqc_data.json / multiqc_sources.json /
  multiqc.logへhost pathの漏洩なし。`multiqc_config.yaml`と`multiqc_version.txt`は
  byte-identical。
- FastQC: 12 report ZIPの`fastqc_data.txt`はunmodified runとbyte-identical。
- GATK / GenomicsDBImport: 73 processすべて成功。GenomicsDB workspace directoryを含め
  work directoryにroot所有entryは0（unmodified runでは208）。published artifactは
  203 entryすべてhost user所有。

`-profile test,docker_amd64`でも同じsyntheticを実行し、73 processすべて成功、work
directoryのroot所有entry 0、published artifactすべてhost user所有を確認しました。生成された
73件の`.command.run`すべてに`-u $(id -u):$(id -g)`と`--platform linux/amd64`の双方が
含まれています。ただしこれはamd64 host上の実行であり、Apple Silicon実機でのemulation経路の
検証ではありません。

### scientific output

publishされたVCF / TSV / TXT 68件のうち53件はunmodified runとbyte-identicalでした。
差分のある15件（VCF 13件とMarkDuplicates metrics 2件）は、unmodified runを2回実行した
場合にも**同一の15件**が差分になります。内訳は次のとおりで、この変更に起因するものでは
ありません。

- GATK / bcftoolsがheaderへ書く実行時刻（`Date=`、`# Started on:`）
- `GenomicsDBImport`のcommand line headerに記録される`--variant`引数順（`.collect()`は
  task完了順に依存するため既存のrun間非決定性）

全13 VCFについて、`##` meta headerを除いたdata body（`#CHROM`行のsample列順を含む）は
unmodified run 2回分と完全に一致しました。GS panel matrix / variant metadata /
sample metadata / accounting TSVはbyte-identicalです。

## Known limitations

- containerの`/etc/passwd`に存在しないUIDで実行するため、`$HOME`は`/`（書き込み不可）に
  解決されます。実害は次の2点のみで、いずれもwork directory内に閉じ、published
  artifactとFastQC計測値には影響しません。
  - FastQC taskのstderrに`Fontconfig error: No writable cache directories`が出る。
  - JavaがuserハンドルとしてFastQC task work directoryへ`?`という名のdirectoryを作る
    （`user.home`は`$HOME`ではなく`getpwuid`由来のため、`-e HOME=...`を渡しても解消
    しません）。これらはhost user所有であり通常権限で削除できます。
- 既に存在するroot所有のwork directory / `.nf-test`作業データの回収はこの変更の対象では
  ありません（Issue #51のnon-goal）。
- Docker以外のcontainer engine（Singularity/Apptainer等）は対象外です。
  `--tmpfs`はDocker固有optionです。
- 検証はsynthetic fixture（2 sample / 3 read group）でのみ実施しています。real cohortや
  327検体規模での再実行は行っていません。MULTIQC tmpfsのmemory影響も同fixture規模でしか
  観測していません。
- `-profile docker_amd64`はamd64 host上でのrun（両option併存とownership）までを確認して
  おり、Apple Silicon実機でのemulation経路の検証は行っていません。
