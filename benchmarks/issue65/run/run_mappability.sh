#!/usr/bin/env bash
# Issue #65 reference self-alignment mappability: 150 bp error-free tiles every 50 bp.
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
REF=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
IDX=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521/results/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
PY='python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'
BWA='community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113@sha256:5ebd1290d9680195817ce75915b79ae2e608834c017824b7e2bc7b141509b242'
U="$(id -u):$(id -g)"
RO="-v /data/adzuki-issue45-20260910/runs:/data/adzuki-issue45-20260910/runs:ro"
cd "$W"
date -Is > mappability.start
docker run --rm -i -u "$U" $RO -v "$W:$W" -w "$W/tools" "$PY" python3 mappability.py reads --reference-fasta "$REF" --read-length 150 --step 50 \
 | docker run --rm -i -u "$U" $RO -v "$W:$W" "$BWA" bash -c "set -o pipefail; bwa-mem2 mem -t 30 $IDX /dev/stdin 2>/data/adzuki-issue45-20260910/issue65-callable-region-20260915/mappability.bwa-mem2.log | samtools view -F 0x900" \
 | docker run --rm -i -u "$U" $RO -v "$W:$W" -w "$W/tools" "$PY" python3 mappability.py summarize --reference-fai "$REF.fai" --read-length 150 --step 50 --min-mapq 30 --output-dir "$W/mappability"
date -Is > mappability.end
