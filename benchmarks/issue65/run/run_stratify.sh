#!/usr/bin/env bash
# Issue #65 descriptive stratification of the Issue #64 GS SNP PASS panel (site filter only)
# and its genotype-masked derivative (site filter + illustrative DP/GQ mask), genome-wide and
# over the callable window.
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
REF=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
PASS=/data/adzuki-issue45-20260910/issue64-genotype-quality-20260914/runs/gs-lineage-e1/results/variants/gs_pass/cohort_gs.snp.pass.vcf.gz
MASKED=/data/adzuki-issue45-20260910/issue64-genotype-quality-20260914/runs/real-checks-2f5ca34/enabled/cohort.gs_panel.quality_masked.vcf.gz
REF_SHA=e9838db1b048b54b21534285aaf95eae64cb3019d85af270b44e20b3d545f383
BCFTOOLS='quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'
PY='python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'
U="$(id -u):$(id -g)"
V=(-u "$U" -v /data/adzuki-issue45-20260910:/data/adzuki-issue45-20260910:ro -v "$W:$W")
S=$W/stratify
cd "$S"
sha256sum "$(readlink -f "$PASS")" "$(readlink -f "$MASKED")" > inputs.sha256
FORMAT='%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n'
for scope in ${SCOPES:-genome_wide window}; do
  region_arg=""; [ "$scope" = window ] && region_arg="-r NC_068975.1:1-20000000"
  rm -f "$S/$scope.masked.fifo"; mkfifo "$S/$scope.masked.fifo"
  docker run --rm "${V[@]}" "$BCFTOOLS" bcftools query $region_arg -f "$FORMAT" -o "$S/$scope.masked.fifo" "$MASKED" &
  masked_pid=$!
  start=$(date +%s)
  docker run --rm "${V[@]}" "$BCFTOOLS" bcftools query $region_arg -f "$FORMAT" "$PASS" \
    | docker run --rm -i "${V[@]}" -w "$W/tools" "$PY" python3 stratify_callset.py \
        --config "$W/configs/gs_snp_pass_$scope.json" --reference-fai "$REF.fai" --reference-sha256 "$REF_SHA" \
        --calls - --masked-calls "$S/$scope.masked.fifo" --output "$S/$scope.json"
  wait "$masked_pid"
  echo -e "$scope\t$(( $(date +%s) - start ))" >> wall_seconds.tsv
  rm -f "$S/$scope.masked.fifo"
done
sha256sum genome_wide.json window.json > outputs.sha256 2>/dev/null || true
