process SAMTOOLS_FAIDX {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'

    publishDir(
        "${params.outdir}/reference",
        mode: 'copy',
        pattern: '*.fai'
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.name}.fai"), emit: fai

    script:
    """
    samtools faidx ${fasta}
    """
}
