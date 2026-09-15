#!/usr/bin/env bash
# Issue #65 sequence-derived region assets on GCF_016808095.1 (as run on seedcore-01, host python3).
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
REF=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
FTP=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/808/095/GCF_016808095.1_ASM1680809v1
cd "$W"
mkdir -p public
( cd public && curl -sS -O "$FTP/GCF_016808095.1_ASM1680809v1_genomic_gaps.txt.gz" && curl -sS -O "$FTP/md5checksums.txt" )
( /usr/bin/time -v python3 tools/build_region_assets.py --reference-fasta "$REF" --output-dir regions \
    --gc-window 1000 --homopolymer-min 10 \
    --ncbi-gaps public/GCF_016808095.1_ASM1680809v1_genomic_gaps.txt.gz ) 2> regions.time
