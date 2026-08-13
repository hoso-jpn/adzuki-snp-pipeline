process GATK_SELECTPASSVARIANTS {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

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
