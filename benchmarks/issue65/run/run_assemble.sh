#!/usr/bin/env bash
# Issue #65: stratify the cross-platform calls, then assemble evaluations and the delivery matrix.
# usage: run_assemble.sh <git sha of the benchmark code in tools/>
set -euo pipefail
GIT_SHA=$1
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
REF=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
REF_SHA=e9838db1b048b54b21534285aaf95eae64cb3019d85af270b44e20b3d545f383
BCFTOOLS='quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'
PY='python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'
U="$(id -u):$(id -g)"
V=(-u "$U" -v /data/adzuki-issue45-20260910:/data/adzuki-issue45-20260910:ro -v "$W:$W")
cd "$W"
if [ -f crossplatform/SRR11787767.vcf.gz ]; then
  docker run --rm "${V[@]}" "$BCFTOOLS" bcftools query -r NC_068975.1:1-20000000 \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n' "$W/crossplatform/SRR11787767.vcf.gz" \
    | docker run --rm -i "${V[@]}" -w "$W/tools" "$PY" python3 stratify_callset.py \
        --config "$W/configs/crossplatform_window.json" --reference-fai "$REF.fai" --reference-sha256 "$REF_SHA" \
        --calls - --output "$W/crossplatform/stratify.json"
fi
docker run --rm "${V[@]}" -w "$W/tools" "$PY" python3 assemble_evidence.py assemble --workdir "$W" \
  --reference-fasta "$REF" --reference-sha256 "$REF_SHA" --git-sha "$GIT_SHA"
