#!/usr/bin/env bash
# Issue #65 callable/depth over the Issue #45 scaling window from the existing 51 markdup BAMs.
# bams.txt (sample order = samples.csv) was written beforehand from the Issue #45 run's work directory:
#   for s in $(tail -n +2 $RUN/samples.csv | cut -d, -f1); do
#     b=$(find work -name "$s.markdup.bam.bai" | head -1); echo "$(dirname $b)/$s.markdup.bam"; done > bams.txt
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
RUN=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521
REGION=NC_068975.1:1-20000000
SAMTOOLS='quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'
PY='python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'
U="$(id -u):$(id -g)"
RO="-v /data/adzuki-issue45-20260910/runs:/data/adzuki-issue45-20260910/runs:ro"
cd "$W"
sed "s#^#$RUN/#" bams.txt > bams.abs.txt
SAMPLES=$(tail -n +2 "$RUN/samples.csv" | cut -d, -f1 | paste -sd,)
date -Is > callable.start
/usr/bin/time -v -o callable.depth.time bash -c "docker run --rm -u '$U' $RO -v '$W:$W' '$SAMTOOLS' samtools depth -a -Q 20 -q 10 -r $REGION -f $W/bams.abs.txt | gzip -1 > depth.window.tsv.gz"
/usr/bin/time -v -o callable.label.time docker run --rm -u "$U" -v "$W:$W" -w "$W/tools" "$PY" python3 callable_regions.py \
  --depth "$W/depth.window.tsv.gz" --samples "$SAMPLES" --region "$REGION" \
  --min-depth 5 --max-depth-factor 2.5 --min-sample-fraction 0.8 --mapq 20 --baseq 10 --output-dir "$W/callable"
date -Is > callable.end
