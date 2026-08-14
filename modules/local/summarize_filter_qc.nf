process SUMMARIZE_FILTER_QC {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_low'

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
