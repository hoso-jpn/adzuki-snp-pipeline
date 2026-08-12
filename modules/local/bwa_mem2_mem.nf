process BWA_MEM2_MEM {
    tag "${meta.id}:${meta.read_group_id}"
    label 'process_high'

    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0@sha256:374f4b910b4b04d32772fccd4cab1cdcd0758356856a960cfa8b1edebfd38c9f'

    publishDir(
        "${params.outdir}/logs/mapping",
        mode: 'copy',
        pattern: '*.bwa-mem2.log'
    )

    input:
    tuple val(meta), path(read1), path(read2)
    tuple val(reference_meta), path(fasta)
    tuple val(index_meta), path(indexes)

    output:
    tuple val(meta), path("${meta.read_group_id}.sam"), emit: sam
    tuple(
        val(meta),
        path("${meta.read_group_id}.bwa-mem2.log"),
        emit: log
    )

    script:
    platform_unit = meta.platform_unit
        ? "\\tPU:${meta.platform_unit}"
        : ''
    read_group = (
        "@RG\\tID:${meta.read_group_id}" +
        "\\tSM:${meta.id}" +
        "\\tLB:${meta.library_id}" +
        "\\tPL:${meta.platform}" +
        platform_unit
    )

    """
    bwa-mem2 mem \
        -t ${task.cpus} \
        -R '${read_group}' \
        ${fasta} \
        ${read1} \
        ${read2} \
        > ${meta.read_group_id}.sam \
        2> ${meta.read_group_id}.bwa-mem2.log
    """
}
