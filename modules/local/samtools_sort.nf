process SAMTOOLS_SORT {
    tag "${meta.id}:${meta.read_group_id}"
    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'

    input:
    tuple val(meta), path(sam)

    output:
    tuple(
        val(meta),
        path("${meta.read_group_id}.sorted.bam"),
        emit: bam
    )

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${meta.read_group_id}.sorted.bam \
        ${sam}
    """
}
