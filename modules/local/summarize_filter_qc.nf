// Issue #30: dedicated resource label (see nextflow.config for the
// real-data benchmark this was sized from). Sized for the SNP
// invocation's real footprint (indel is far lighter -- 1.4 GiB true
// peak vs. SNP's 8.4 GiB on Issue #26's real 5-sample cohort -- but
// both invocations share this one process/label; splitting further by
// meta.variant_type was considered and rejected as unnecessary
// complexity for a single label difference not yet shown to matter).
process SUMMARIZE_FILTER_QC {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_variant_qc_summary'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(filtered_vcf), path(filtered_vcf_index)

    output:
    tuple(
        val(meta),
        path("${meta.id}.filtered.${meta.variant_type}.filter_breakdown.tsv"),
        path("${meta.id}.filtered.${meta.variant_type}.annotation_qc.tsv"),
        path("${meta.id}.filtered.${meta.variant_type}.filter_qc.summary.txt"),
        emit: qc
    )

    script:
    prefix = "${meta.id}.filtered.${meta.variant_type}"

    """
    summarize_filter_qc.py \
        --filtered-vcf ${filtered_vcf} \
        --cohort-id '${meta.id}' \
        --stage 'filtered' \
        --variant-type '${meta.variant_type}' \
        --filter-breakdown-output ${prefix}.filter_breakdown.tsv \
        --annotation-qc-output ${prefix}.annotation_qc.tsv \
        --summary-output ${prefix}.filter_qc.summary.txt
    """
}
