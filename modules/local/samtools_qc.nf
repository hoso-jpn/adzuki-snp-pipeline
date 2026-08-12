process SAMTOOLS_QC {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'

    publishDir(
        "${params.outdir}/qc/samtools",
        mode: 'copy',
        pattern: '*.txt'
    )

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    tuple val(meta), path("${meta.id}.stats.txt"), emit: stats
    tuple val(meta), path("${meta.id}.idxstats.txt"), emit: idxstats

    script:
    """
    samtools quickcheck --verbose ${bam}

    samtools flagstat \
        --threads ${task.cpus} \
        ${bam} \
        > ${meta.id}.flagstat.txt

    samtools stats \
        --threads ${task.cpus} \
        ${bam} \
        > ${meta.id}.stats.txt

    samtools idxstats \
        ${bam} \
        > ${meta.id}.idxstats.txt
    """
}
