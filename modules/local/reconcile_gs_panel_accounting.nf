process RECONCILE_GS_PANEL_ACCOUNTING {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

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
    // Issue #52: see modules/local/gs_normalize_variants.nf for why this
    // records task.container (the resolved, post-override effective
    // container) rather than trusting the `container` directive above.
    val(task.container), emit: container_id

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
