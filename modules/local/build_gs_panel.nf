process BUILD_GS_PANEL {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(gs_pass_vcf), path(gs_pass_vcf_index)

    output:
    tuple(val(meta), path("${meta.id}.gs_panel.genotype_matrix.tsv.gz"), emit: matrix)
    path("${meta.id}.gs_panel.sample_metadata.tsv"), emit: sample_metadata
    path("${meta.id}.gs_panel.variant_metadata.tsv"), emit: variant_metadata
    path("${meta.id}.gs_panel.genotype_encoding_accounting.tsv"), emit: genotype_accounting
    path("${meta.id}.gs_panel.genotype_encoding_accounting.summary.txt"), emit: genotype_accounting_summary

    script:
    prefix = "${meta.id}.gs_panel"

    """
    build_gs_panel.py \
        --gs-pass-vcf ${gs_pass_vcf} \
        --cohort-id '${meta.id}' \
        --sample-ploidy ${params.sample_ploidy} \
        --matrix-output ${prefix}.genotype_matrix.tsv.gz \
        --sample-metadata-output ${prefix}.sample_metadata.tsv \
        --variant-metadata-output ${prefix}.variant_metadata.tsv \
        --genotype-accounting-output ${prefix}.genotype_encoding_accounting.tsv \
        --genotype-accounting-summary-output ${prefix}.genotype_encoding_accounting.summary.txt
    """
}
