process RECONCILE_GS_PANEL_ACCOUNTING {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(raw_all_vcf), path(raw_all_vcf_index)
    tuple val(normalized_meta), path(normalized_vcf), path(normalized_vcf_index)
    path(normalization_accounting)
    tuple val(gs_pass_meta), path(gs_pass_vcf), path(gs_pass_vcf_index)
    tuple val(matrix_meta), path(matrix)
    path(variant_metadata)
    path(sample_metadata)

    output:
    tuple(
        val(meta),
        path("${meta.id}.gs_panel.record_accounting.tsv"),
        path("${meta.id}.gs_panel.record_accounting.summary.txt"),
        emit: accounting
    )

    script:
    """
    reconcile_gs_panel_accounting.py \
        --cohort-id '${meta.id}' \
        --raw-all-vcf ${raw_all_vcf} \
        --normalized-vcf ${normalized_vcf} \
        --normalization-accounting ${normalization_accounting} \
        --gs-pass-vcf ${gs_pass_vcf} \
        --matrix ${matrix} \
        --variant-metadata ${variant_metadata} \
        --sample-metadata ${sample_metadata} \
        --output ${meta.id}.gs_panel.record_accounting.tsv \
        --summary-output ${meta.id}.gs_panel.record_accounting.summary.txt
    """
}
