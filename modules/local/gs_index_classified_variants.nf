process GS_INDEX_CLASSIFIED_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    container params.containers.bcftools

    input:
    tuple val(meta), path(classified_vcf)

    output:
    tuple(
        val(meta),
        path("${meta.id}.classified.vcf.gz"),
        path("${meta.id}.classified.vcf.gz.tbi"),
        emit: vcf
    )
    // Issue #52: see gs_normalize_variants.nf for why this records
    // task.container (the resolved, post-override effective container)
    // rather than trusting the `container` directive above.
    val(task.container), emit: container_id

    script:
    """
    bcftools view --output-type z --output ${meta.id}.classified.vcf.gz ${classified_vcf}
    bcftools index --tbi ${meta.id}.classified.vcf.gz
    """
}
