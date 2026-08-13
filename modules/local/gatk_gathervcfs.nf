process GATK_GATHERVCFS {
    tag "${meta.id}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    publishDir(
        "${params.outdir}/variants/raw",
        mode: 'copy',
        pattern: "${meta.id}.raw.vcf.gz*"
    )

    input:
    tuple(
        val(meta),
        path(vcfs),
        path(vcf_indexes)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.raw.vcf.gz"),
        path("${meta.id}.raw.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )
    input_arguments = vcfs
        .collect { vcf -> "--INPUT ${vcf}" }
        .join(" \\\n        ")

    if (vcfs.size() != vcf_indexes.size()) {
        error(
            'the number of VCFs and indexes must match: ' +
            "${vcfs.size()} VCFs, " +
            "${vcf_indexes.size()} indexes"
        )
    }

    """
    gatk --java-options "-Xmx${memory_gb}g" GatherVcfs \
        ${input_arguments} \
        --OUTPUT ${meta.id}.raw.vcf.gz \
        --CREATE_INDEX true
    """
}
