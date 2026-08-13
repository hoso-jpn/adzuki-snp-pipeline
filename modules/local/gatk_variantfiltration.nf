process GATK_VARIANTFILTRATION {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple(
        val(meta),
        path(vcf),
        path(vcf_index),
        val(filters)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.variant_type}.filtered.vcf.gz"),
        path("${meta.id}.${meta.variant_type}.filtered.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    if (filters.isEmpty()) {
        error(
            "no hard filters were configured for ${meta.variant_type}"
        )
    }

    filter_arguments = filters
        .collect { filter ->
            def name = filter['name']
            def expression = filter['expression']

            if (!name || !expression) {
                error(
                    "invalid hard filter for ${meta.variant_type}: ${filter}"
                )
            }

            "--filter-name '${name}' " +
                "--filter-expression '${expression}'"
        }
        .join(" \\\n        ")

    """
    gatk --java-options "-Xmx${memory_gb}g" VariantFiltration \
        --variant ${vcf} \
        ${filter_arguments} \
        --output ${meta.id}.${meta.variant_type}.filtered.vcf.gz \
        --create-output-variant-index true
    """
}
