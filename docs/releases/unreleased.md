# Unreleased

次のversionを切る際のrelease notes素材です。`v0.2.0`以降にmergeされた、release notesへ
載せるべき変更を記録します。versionを切る時点で`docs/releases/<version>.md`へ取り込みます。

## Runtime behaviour

- Docker実行のcontainer taskが、container内rootではなく起動したhost user
  (`-u $(id -u):$(id -g)`) として実行されるようになりました。work directoryへ書かれる
  fileが実行ユーザー所有になり、失敗runや`.nf-test`作業データを通常権限で回収できます。
  `docker`と`docker_amd64`の双方に適用され、`docker_amd64`では`--platform linux/amd64`と
  併存します。MULTIQCのcontainer絶対staging root `/multiqc` は、published provenanceの
  source path契約を維持したままtmpfs mountで書き込み可能にしています。scientific
  semantics、threshold、schema、artifact naming、pin済みtool/container version、resource
  policyは変更していません。実行契約・検証内容・既知の制約は
  [`../docker_host_user_execution.md`](../docker_host_user_execution.md)を参照してください。
  既に存在するroot所有artifactの回収はこの変更の対象外です。(Issue #51)
