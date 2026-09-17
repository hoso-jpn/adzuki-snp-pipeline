#!/usr/bin/env bash
# Issue #65: downsampling stability and caller concordance on one existing BAM.
# Sample: SRR29908806, the most deeply sequenced of the 51 (samtools stats bases mapped).
# Window: NC_068975.1:1-20000000. Nothing here reprocesses FASTQ.
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
RUN=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521
REF=$RUN/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
SAMPLE=SRR29908806
REGION=NC_068975.1:1-20000000
SAMTOOLS='quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'
GATK='broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'
BCFTOOLS='quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'
U="$(id -u):$(id -g)"
D=$W/downsampling
M=$W/measured.sh
V=(-u "$U" -v /data/adzuki-issue45-20260910/runs:/data/adzuki-issue45-20260910/runs:ro -v "$W:$W" -w "$D")
cd "$D"
BAM=$(sed "s#^#$RUN/#" "$W/bams.txt" | grep "/$SAMPLE.markdup.bam$")
sha256sum "$(readlink -f "$BAM")" "$(readlink -f "$BAM.bai")" "$REF" > inputs.sha256
date -Is > start

# 1. window BAMs: full, and read-name-hash subsamples (mates stay together) at fixed seeds
declare -A FRACTION=([full]=1 [f050_s65]=0.5 [f050_s66]=0.5 [f025_s65]=0.25 [f025_s66]=0.25)
declare -A SEED=([full]=0 [f050_s65]=65 [f050_s66]=66 [f025_s65]=65 [f025_s66]=66)
for id in full f050_s65 f050_s66 f025_s65 f025_s66; do
  if [ "$id" = full ]; then sub=""; else sub="--subsample ${FRACTION[$id]} --subsample-seed ${SEED[$id]}"; fi
  $M "$id.view" "${V[@]}" "$SAMTOOLS" sh -c "samtools view -b $sub -o $id.bam $BAM $REGION && samtools index $id.bam && samtools view -c $id.bam > $id.reads.txt && samtools depth -a -Q 20 -q 10 -r $REGION $id.bam | awk '{s+=\$3} END {printf \"%.4f\\n\", s/NR}' > $id.mean_depth.txt"
done

# 2. production HaplotypeCaller arguments, restricted to the window, 4 threads each, in parallel
pids=()
for id in full f050_s65 f050_s66 f025_s65 f025_s66; do
  $M "$id.hc" "${V[@]}" --cpus 4 "$GATK" gatk --java-options -Xmx15g HaplotypeCaller \
    --reference "$REF" --input "$id.bam" --output "$id.g.vcf.gz" --emit-ref-confidence GVCF \
    --sample-ploidy 2 --native-pair-hmm-threads 4 --create-output-variant-index true --intervals "$REGION" &
  pids+=($!)
done
bcf_call() {
  $M full.bcftools "${V[@]}" --cpus 2 "$BCFTOOLS" sh -c "bcftools mpileup --no-version -f $REF -r $REGION -q 20 -Q 10 -a FORMAT/AD,FORMAT/DP -Ou full.bam | bcftools call --no-version -m -v -Oz -o full.bcftools.vcf.gz && bcftools index -t full.bcftools.vcf.gz"
}
bcf_call &
pids+=($!)
for p in "${pids[@]}"; do wait "$p"; done

# 3. single-sample GenotypeGVCFs on each GVCF, same window
pids=()
for id in full f050_s65 f050_s66 f025_s65 f025_s66; do
  $M "$id.gg" "${V[@]}" --cpus 2 "$GATK" gatk --java-options -Xmx7g GenotypeGVCFs \
    --reference "$REF" --variant "$id.g.vcf.gz" --intervals "$REGION" --sample-ploidy 2 \
    --output "$id.vcf.gz" --create-output-variant-index true &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
sha256sum *.vcf.gz *.g.vcf.gz > outputs.sha256
du -cb *.bam *.bai *.g.vcf.gz* *.vcf.gz* | tail -1 > storage_bytes.txt
date -Is > end
