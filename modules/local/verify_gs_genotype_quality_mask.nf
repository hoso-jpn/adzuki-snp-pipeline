// Issue #64: check that the quality-masked VCF, the matrix, both metadata
// files and the genotype accounting all describe one masked state, and that
// each call was masked exactly when the recorded policy masks it. A mismatch
// fails the run; BUILD_GS_PANEL_MANIFEST waits for this, so no manifest is
// ever written for a masked panel that did not verify.
process VERIFY_GS_GENOTYPE_QUALITY_MASK {
    tag "${meta.id}"
    // Bounded memory (one row of each input plus per-sample counters), so the
    // generic low tier; its real 51-sample wall time is recorded in
    // docs/gs_panel_genotype_quality_mask.md.
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

    input:
    tuple val(meta), path(gs_pass_vcf), path(gs_pass_vcf_index)
    tuple val(masked_meta), path(quality_masked_vcf), path(quality_masked_vcf_index)
    tuple val(matrix_meta), path(matrix)
    path(sample_metadata)
    path(variant_metadata)
    path(genotype_accounting)
    path(genotype_quality_policy)

    output:
    tuple(
        val(meta),
        path("${meta.id}.gs_panel.genotype_quality_mask_verification.tsv"),
        path("${meta.id}.gs_panel.genotype_quality_mask_verification.summary.txt"),
        emit: verification
    )
    // Issue #52: see gs_normalize_variants.nf for why this records
    // task.container (the resolved, post-override effective container).
    val(task.container), emit: container_id

    script:
    """
    verify_gs_genotype_quality_mask.py \
        --cohort-id '${meta.id}' \
        --genotype-quality-policy ${genotype_quality_policy} \
        --gs-pass-vcf ${gs_pass_vcf} \
        --quality-masked-vcf ${quality_masked_vcf} \
        --matrix ${matrix} \
        --variant-metadata ${variant_metadata} \
        --sample-metadata ${sample_metadata} \
        --genotype-accounting ${genotype_accounting} \
        --output ${meta.id}.gs_panel.genotype_quality_mask_verification.tsv \
        --summary-output ${meta.id}.gs_panel.genotype_quality_mask_verification.summary.txt
    """
}
