#!/usr/bin/env bash
# Issue #65 cross-platform self-consistency: Illumina WGS reads (SRR11787767, a leading byte range)
# of BioSample SAMN14776547, the same BioSample as the PacBio-derived reference GCF_016808095.1,
# through the production fastp / BWA-MEM2 / MarkDuplicates / HaplotypeCaller arguments, called on
# NC_068975.1:1-20000000. Calls describe where these reads disagree with the assembly consensus.
set -euo pipefail
W=/data/adzuki-issue45-20260910/issue65-callable-region-20260915
RUN=/data/adzuki-issue45-20260910/runs/issue45-51gvcf-20260910-200521
REF=$RUN/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
IDX=$RUN/results/reference/GCF_016808095.1_ASM1680809v1_genomic.fna
REGION=NC_068975.1:1-20000000
ID=SRR11787767
FASTP='quay.io/biocontainers/fastp:1.3.6--h43da1c4_0@sha256:cbbe2402b6b6704df470d7d77dcb498eefd5bcd01f4c38be0ec69899e79ac134'
BWA='community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113@sha256:5ebd1290d9680195817ce75915b79ae2e608834c017824b7e2bc7b141509b242'
SAMTOOLS='quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'
GATK='broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'
U="$(id -u):$(id -g)"
M=$W/measured.sh
C=$W/crossplatform
IN=$W/input/$ID
V=(-u "$U" -v /data/adzuki-issue45-20260910/runs:/data/adzuki-issue45-20260910/runs:ro -v "$W:$W" -w "$C")
cd "$C"
until grep -q '^done' "$IN/download.log"; do sleep 30; done
date -Is > start
sha256sum "$IN/partial_1.fastq.gz" "$IN/partial_2.fastq.gz" > partial.sha256

# 1. keep whole records present in both truncated mates, and prove the mates pair up
n1=$( (gzip -dc "$IN/partial_1.fastq.gz" 2>/dev/null || true) | wc -l)
n2=$( (gzip -dc "$IN/partial_2.fastq.gz" 2>/dev/null || true) | wc -l)
pairs=$(( (n1 < n2 ? n1 : n2) / 4 - 1 ))
for m in 1 2; do (gzip -dc "$IN/partial_$m.fastq.gz" 2>/dev/null || true) | head -n $((pairs * 4)) > "$ID.R$m.fastq"; done
paste <(awk 'NR%4==1{print $1}' "$ID.R1.fastq") <(awk 'NR%4==1{print $1}' "$ID.R2.fastq") | awk '$1!=$2{bad++} END{if (bad) {print "unpaired names:", bad; exit 1}}'
awk 'NR%4==2{b+=length($0)} END{print NR/4, b}' "$ID.R1.fastq" > read_pairs_and_r1_bases.txt
sha256sum "$ID.R1.fastq" "$ID.R2.fastq" > fastq_used.sha256

# 2. production fastp arguments
$M fastp "${V[@]}" --cpus 8 "$FASTP" fastp --in1 "$ID.R1.fastq" --in2 "$ID.R2.fastq" \
  --out1 "$ID.trimmed_R1.fastq.gz" --out2 "$ID.trimmed_R2.fastq.gz" --json "$ID.fastp.json" --html "$ID.fastp.html" \
  --detect_adapter_for_pe --thread 8
rm -f "$ID.R1.fastq" "$ID.R2.fastq"

# 3. production BWA-MEM2 | samtools sort pipe
$M bwa "${V[@]}" --cpus 28 "$BWA" bash -c "set -o pipefail; bwa-mem2 mem -t 24 -R '@RG\tID:$ID\tSM:$ID\tLB:1-ye-1\tPL:ILLUMINA' $IDX $ID.trimmed_R1.fastq.gz $ID.trimmed_R2.fastq.gz 2> $ID.bwa-mem2.log | samtools sort -@ 4 -m 2000M -o $ID.sorted.bam -"

# 4. production MarkDuplicates arguments, then index
$M markdup "${V[@]}" --cpus 4 "$GATK" gatk --java-options -Xmx15g MarkDuplicates --INPUT "$ID.sorted.bam" \
  --OUTPUT "$ID.markdup.bam" --METRICS_FILE "$ID.markduplicates.metrics.txt" --REMOVE_DUPLICATES false \
  --CREATE_INDEX false --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --VALIDATION_STRINGENCY STRICT
rm -f "$ID.sorted.bam"
$M index "${V[@]}" "$SAMTOOLS" sh -c "samtools index $ID.markdup.bam && samtools depth -a -Q 20 -q 10 -r $REGION $ID.markdup.bam | awk '{s+=\$3; if (\$3>=5) c++} END {printf \"%.4f\t%d\n\", s/NR, c}' > $ID.window_mean_depth_and_ge5.txt"

# 5. production HaplotypeCaller arguments on the window, then GenotypeGVCFs
$M hc "${V[@]}" --cpus 4 "$GATK" gatk --java-options -Xmx15g HaplotypeCaller --reference "$REF" --input "$ID.markdup.bam" \
  --output "$ID.g.vcf.gz" --emit-ref-confidence GVCF --sample-ploidy 2 --native-pair-hmm-threads 4 \
  --create-output-variant-index true --intervals "$REGION"
$M gg "${V[@]}" --cpus 2 "$GATK" gatk --java-options -Xmx7g GenotypeGVCFs --reference "$REF" --variant "$ID.g.vcf.gz" \
  --intervals "$REGION" --sample-ploidy 2 --output "$ID.vcf.gz" --create-output-variant-index true
sha256sum "$ID.markdup.bam" "$ID.g.vcf.gz" "$ID.vcf.gz" > outputs.sha256
du -cb "$ID".* | tail -1 > storage_bytes.txt
date -Is > end
