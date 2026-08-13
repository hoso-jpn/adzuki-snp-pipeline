process GATK_SELECTVARIANTS {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple(
        val(meta),
        path(vcf),
        path(vcf_index),
        val(variant_type)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.variant_type}.vcf.gz"),
        path("${meta.id}.${meta.variant_type}.vcf.gz.tbi"),
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
        --select-type-to-include ${variant_type} \
        --output ${meta.id}.${meta.variant_type}.vcf.gz \
        --create-output-variant-index true
    """
}
