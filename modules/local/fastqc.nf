process FASTQC {
    tag "${stage}:${meta.id}:${meta.read_group_id}"
    label 'process_low'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0@sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900'

    publishDir(
        "${params.outdir}/qc/fastqc/${stage}",
        mode: 'copy',
        pattern: '*_fastqc.{html,zip}'
    )

    input:
    tuple val(meta), path(reads)
    val stage

    output:
    tuple val(meta), path('*_fastqc.html'), emit: html
    tuple val(meta), path('*_fastqc.zip'), emit: zip

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads}
    """
}
