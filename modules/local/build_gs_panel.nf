// Issue #30: dedicated resource label (see nextflow.config for the
// real-data benchmark this was sized from). On Issue #26's real
// 5-sample cohort this process's true peak RSS (5.34 GiB, measured at
// a generous memory ceiling) exceeded process_low's previous 4 GiB
// first-attempt allocation -- the original run's own report of
// "peak_rss: 4 GB" was the cgroup ceiling itself, not genuine headroom.
process BUILD_GS_PANEL {
    tag "${meta.id}"
    label 'process_gs_panel'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

    input:
    tuple val(meta), path(gs_pass_vcf), path(gs_pass_vcf_index)

    output:
    tuple(val(meta), path("${meta.id}.gs_panel.genotype_matrix.tsv.gz"), emit: matrix)
    path("${meta.id}.gs_panel.sample_metadata.tsv"), emit: sample_metadata
    path("${meta.id}.gs_panel.variant_metadata.tsv"), emit: variant_metadata
    path("${meta.id}.gs_panel.genotype_encoding_accounting.tsv"), emit: genotype_accounting
    path("${meta.id}.gs_panel.genotype_encoding_accounting.summary.txt"), emit: genotype_accounting_summary
    // Issue #52: see modules/local/gs_normalize_variants.nf for why this
    // records task.container (the resolved, post-override effective
    // container) rather than trusting the `container` directive above.
    val(task.container), emit: container_id

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
