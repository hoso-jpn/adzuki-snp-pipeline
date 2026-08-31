process GATK_SELECTPASSVARIANTS {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container params.containers.gatk

    input:
    tuple(
        val(meta),
        path(vcf),
        path(vcf_index)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.variant_type}.pass.vcf.gz"),
        path("${meta.id}.${meta.variant_type}.pass.vcf.gz.tbi"),
        emit: vcf
    )
    // Issue #52: see gatk_variantfiltration.nf for why this records
    // task.container (the resolved, per-alias effective container) rather
    // than the `container` directive's default -- this module is aliased
    // as both GATK_SELECTPASSVARIANTS and GATK_SELECTPASSVARIANTS_GS.
    val(task.container), emit: container_id

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    """
    gatk --java-options "-Xmx${memory_gb}g" SelectVariants \
        --variant ${vcf} \
        --exclude-filtered true \
        --output ${meta.id}.${meta.variant_type}.pass.vcf.gz \
        --create-output-variant-index true
    """
}
