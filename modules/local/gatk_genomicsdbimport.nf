process GATK_GENOMICSDBIMPORT {
    tag "${interval_meta.id}"
    label 'process_high'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple val(interval_meta), val(interval)
    path(gvcfs)
    path(gvcf_indexes)

    output:
    tuple(
        val(interval_meta),
        val(interval),
        path("${interval_meta.id}.genomicsdb"),
        emit: genomicsdb
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )
    variant_arguments = gvcfs
        .collect { gvcf -> "--variant ${gvcf}" }
        .join(" \\\n        ")

    if (gvcfs.size() != gvcf_indexes.size()) {
        error(
            'the number of gVCFs and indexes must match: ' +
            "${gvcfs.size()} gVCFs, " +
            "${gvcf_indexes.size()} indexes"
        )
    }

    """
    gatk --java-options "-Xmx${memory_gb}g" GenomicsDBImport \
        ${variant_arguments} \
        --genomicsdb-workspace-path ${interval_meta.id}.genomicsdb \
        --intervals '${interval}' \
        --reader-threads ${task.cpus} \
        --tmp-dir .
    """
}
