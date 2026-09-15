#!/bin/bash
# Issue #64: real-cohort checks on the 51-sample GS-eligible PASS VCF, run in
# order and one at a time so no measurement contends with another.
#
#   1. disabled compatibility: the branch's builder with the mask off must
#      reproduce main's five panel outputs byte for byte;
#   2. manifest compatibility: main's and the branch's GS manifest for the same
#      unmasked panel must agree on everything except the schema v3 change;
#   3. the default reject policy against real calls;
#   4. one illustrative enabled build (not a recommendation), indexed and
#      verified, with time and peak RSS;
#   5. the threshold-grid sensitivity analysis, cross-checked against step 4.
#
# Usage: run_real_cohort_checks.sh <issue64 base dir> <main checkout> <feature checkout> <out dir>
set -euo pipefail
BASE=$1; MAIN=$2; FEATURE=$3; OUT=$4
LINEAGE=$BASE/runs/gs-lineage-e1
PASS_VCF=$LINEAGE/results/variants/gs_pass/cohort_gs.snp.pass.vcf.gz
BASELINE=$LINEAGE/results/gs_panel
PY=$(grep -o "python  : '[^']*'" "$FEATURE/conf/containers.config" | cut -d"'" -f2)
BC=$(grep -o "bcftools: '[^']*'" "$FEATURE/conf/containers.config" | cut -d"'" -f2)
mkdir -p "$OUT"
cd "$OUT"

run_py() {  # run_py <label> <checkout> <script and args...>; measured inside the pinned container
  local label=$1 checkout=$2; shift 2
  docker run --rm -u "$(id -u):$(id -g)" --network none \
    -v "$BASE:$BASE" -v "$OUT:$OUT" -w "$OUT" \
    -v "$FEATURE/benchmarks/issue45/measure_process.py:/measure_process.py:ro" \
    "$PY" python3 /measure_process.py "$OUT/$label.metrics.json" python3 "$checkout/bin/$@" \
    > "$OUT/$label.stdout" 2> "$OUT/$label.stderr"
}

panel_args() {  # panel_args <dir>
  echo --gs-pass-vcf "$PASS_VCF" --cohort-id cohort --sample-ploidy 2 \
    --matrix-output "$1/cohort.gs_panel.genotype_matrix.tsv.gz" \
    --sample-metadata-output "$1/cohort.gs_panel.sample_metadata.tsv" \
    --variant-metadata-output "$1/cohort.gs_panel.variant_metadata.tsv" \
    --genotype-accounting-output "$1/cohort.gs_panel.genotype_encoding_accounting.tsv" \
    --genotype-accounting-summary-output "$1/cohort.gs_panel.genotype_encoding_accounting.summary.txt"
}

echo "== 1 disabled compatibility $(date -Is)"
mkdir -p disabled
run_py disabled_build "$FEATURE" build_gs_panel.py $(panel_args "$OUT/disabled")
for name in genotype_matrix.tsv.gz sample_metadata.tsv variant_metadata.tsv \
            genotype_encoding_accounting.tsv genotype_encoding_accounting.summary.txt; do
  if cmp -s "$BASELINE/cohort.gs_panel.$name" "disabled/cohort.gs_panel.$name"; then
    echo "identical cohort.gs_panel.$name"
  else
    echo "DIFFERENT cohort.gs_panel.$name"; exit 1
  fi
done

echo "== 2 manifest compatibility $(date -Is)"
MANIFEST_COMMON=(--cohort-id cohort --pipeline-version 0.2.0 --git-commit ''
  --container-gs-normalize-variants a --container-classify-normalized-variants b
  --container-gs-index-classified-variants c --container-gatk-variantfiltration-gs d
  --container-gatk-selectpassvariants-gs e --container-build-gs-panel f
  --container-reconcile-gs-panel-accounting g --container-build-gs-panel-manifest h
  --sample-ploidy 2 --snp-filter-qd-min 2.0 --snp-filter-qual-min 30.0 --snp-filter-sor-max 3.0
  --snp-filter-fs-max 60.0 --snp-filter-mq-min 40.0 --snp-filter-mq-rank-sum-min -12.5
  --snp-filter-read-pos-rank-sum-min -8.0
  --record-accounting "$BASELINE/cohort.gs_panel.record_accounting.tsv")
for name in genotype_matrix.tsv.gz sample_metadata.tsv variant_metadata.tsv \
            genotype_encoding_accounting.tsv genotype_encoding_accounting.summary.txt \
            record_accounting.tsv record_accounting.summary.txt; do
  MANIFEST_COMMON+=(--checksum-file "$BASELINE/cohort.gs_panel.$name")
done
MANIFEST_COMMON+=(--checksum-file "$PASS_VCF")
run_py manifest_main "$MAIN" build_gs_panel_manifest.py "${MANIFEST_COMMON[@]}" --output "$OUT/manifest_main.json"
run_py manifest_feature "$FEATURE" build_gs_panel_manifest.py "${MANIFEST_COMMON[@]}" --no-genotype-quality-mask --output "$OUT/manifest_feature.json"
python3 - "$OUT/manifest_main.json" "$OUT/manifest_feature.json" > manifest_compatibility.json <<'EOF'
import json, sys
main, feature = (json.load(open(p)) for p in sys.argv[1:])
volatile = {"run_id", "generated_at", "manifest_hash"}
changed = {"schema_version", "genotype_quality_mask"}
m = {k: v for k, v in main.items() if k not in volatile}
f = {k: v for k, v in feature.items() if k not in volatile | changed}
result = {
    "main_schema_version": main["schema_version"],
    "feature_schema_version": feature["schema_version"],
    "keys_only_in_main": sorted(set(m) - set(f) - changed),
    "keys_only_in_feature": sorted(set(feature) - set(main) - volatile),
    "every_other_field_equal": {k: v for k, v in m.items() if k not in changed} == f,
    "checksums_equal": main["checksums"] == feature["checksums"],
    "feature_genotype_quality_mask": feature["genotype_quality_mask"],
}
print(json.dumps(result, indent=2))
sys.exit(0 if result["every_other_field_equal"] and result["keys_only_in_main"] == [] and result["keys_only_in_feature"] == ["genotype_quality_mask"] else 1)
EOF
cat manifest_compatibility.json

echo "== 3 default reject policy $(date -Is)"
mkdir -p reject
set +e
run_py reject_build "$FEATURE" build_gs_panel.py $(panel_args "$OUT/reject") \
  --genotype-quality-mask --genotype-min-dp 10 --genotype-min-gq 20 \
  --quality-masked-vcf-output "$OUT/reject/cohort.gs_panel.quality_masked.vcf.gz" \
  --genotype-quality-policy-output "$OUT/reject/cohort.gs_panel.genotype_quality_policy.json"
echo "reject exit=$?"
set -e
tail -2 reject_build.stderr || true
ls reject

echo "== 4 illustrative enabled build $(date -Is)"
mkdir -p enabled
run_py enabled_build "$FEATURE" build_gs_panel.py $(panel_args "$OUT/enabled") \
  --genotype-quality-mask --genotype-min-dp 10 --genotype-min-gq 20 \
  --genotype-missing-format-field-policy unevaluated --genotype-missing-value-policy unevaluated \
  --genotype-malformed-value-policy unevaluated \
  --quality-masked-vcf-output "$OUT/enabled/cohort.gs_panel.quality_masked.vcf.gz" \
  --genotype-quality-policy-output "$OUT/enabled/cohort.gs_panel.genotype_quality_policy.json"
/usr/bin/time -v -o index.time docker run --rm -u "$(id -u):$(id -g)" --network none -v "$OUT:$OUT" -w "$OUT/enabled" \
  "$BC" bcftools index --tbi cohort.gs_panel.quality_masked.vcf.gz
run_py enabled_verify "$FEATURE" verify_gs_genotype_quality_mask.py --cohort-id cohort \
  --genotype-quality-policy "$OUT/enabled/cohort.gs_panel.genotype_quality_policy.json" \
  --gs-pass-vcf "$PASS_VCF" --quality-masked-vcf "$OUT/enabled/cohort.gs_panel.quality_masked.vcf.gz" \
  --matrix "$OUT/enabled/cohort.gs_panel.genotype_matrix.tsv.gz" \
  --variant-metadata "$OUT/enabled/cohort.gs_panel.variant_metadata.tsv" \
  --sample-metadata "$OUT/enabled/cohort.gs_panel.sample_metadata.tsv" \
  --genotype-accounting "$OUT/enabled/cohort.gs_panel.genotype_encoding_accounting.tsv" \
  --output "$OUT/enabled/cohort.gs_panel.genotype_quality_mask_verification.tsv" \
  --summary-output "$OUT/enabled/cohort.gs_panel.genotype_quality_mask_verification.summary.txt"
cat enabled/cohort.gs_panel.genotype_quality_mask_verification.tsv
# The GS-eligible PASS VCF itself must be untouched by any of this.
sha256sum "$PASS_VCF" > pass_vcf_after.sha256

echo "== 5 sensitivity grid $(date -Is)"
docker run --rm -u "$(id -u):$(id -g)" --network none -v "$BASE:$BASE" -v "$OUT:$OUT" -w "$OUT" \
  -v "$FEATURE/benchmarks/issue45/measure_process.py:/measure_process.py:ro" \
  "$PY" python3 /measure_process.py "$OUT/sensitivity.metrics.json" \
  python3 "$FEATURE/benchmarks/issue64/genotype_quality_sensitivity.py" \
  --gs-pass-vcf "$PASS_VCF" --dp-thresholds 3,5,8,10,15,20 --gq-thresholds 10,20,30 \
  --output "$OUT/sensitivity.json" \
  --cross-check-accounting "$OUT/enabled/cohort.gs_panel.genotype_encoding_accounting.tsv" \
  --cross-check-min-dp 10 --cross-check-min-gq 20 > sensitivity.stdout 2> sensitivity.stderr
echo "== done $(date -Is)"
