// Issue #64: index the quality-masked VCF that BUILD_GS_PANEL wrote.
//
// Deliberately index-only. BUILD_GS_PANEL already writes the VCF as BGZF in
// the same pass as the matrix, so re-serializing it here (bcftools view -Oz)
// would put a second tool between the masking decision and the published
// bytes. `bcftools index` also rejects anything that is not valid BGZF, so
// this step doubles as a format check inside the pinned bcftools container.
process GS_INDEX_QUALITY_MASKED_VCF {
    tag "${meta.id}"
    label 'process_low'

    container params.containers.bcftools

    input:
    tuple val(meta), path(quality_masked_vcf)

    output:
    tuple(
        val(meta),
        path(quality_masked_vcf, includeInputs: true),
        path("${quality_masked_vcf}.tbi"),
        emit: vcf
    )
    // Issue #52: see gs_normalize_variants.nf for why this records
    // task.container (the resolved, post-override effective container).
    val(task.container), emit: container_id

    script:
    """
    bcftools index --tbi ${quality_masked_vcf}
    """
}
