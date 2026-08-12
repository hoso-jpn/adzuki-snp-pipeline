process BWA_MEM2_INDEX {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0@sha256:374f4b910b4b04d32772fccd4cab1cdcd0758356856a960cfa8b1edebfd38c9f'

    publishDir(
        "${params.outdir}/reference",
        mode: 'copy',
        pattern: '*.{0123,amb,ann,bwt.2bit.64,pac}'
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.name}.*"), emit: indexes

    script:
    """
    bwa-mem2 index ${fasta}
    """
}
