#!/bin/bash
set -uo pipefail
for m in 1 2; do
  url=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/067/SRR11787767/SRR11787767_${m}.fastq.gz
  echo "start $m $(date -Is)"
  for attempt in 1 2 3 4 5 6 7 8 9 10; do
    have=$(stat -c %s partial_${m}.fastq.gz 2>/dev/null || echo 0)
    [ "$have" -ge 1600000000 ] && break
    curl -sS --connect-timeout 30 --retry 3 -r ${have}-1599999999 "$url" >> partial_${m}.fastq.gz
    echo "attempt $attempt exit=$? bytes=$(stat -c %s partial_${m}.fastq.gz) $(date -Is)"
  done
done
sha256sum partial_1.fastq.gz partial_2.fastq.gz > partial.sha256
echo "done $(date -Is)"
