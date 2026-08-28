process FASTQC {
    tag "${stage}:${meta.id}:${meta.read_group_id}"
    label 'process_low'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0@sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900'

    input:
    tuple val(meta), path(reads)
    val stage

    output:
    tuple val(meta), path('*_fastqc.html'), emit: html
    tuple val(meta), path('*_fastqc.zip'), emit: zip

    script:
    def read1_name = "${meta.read_group_id}.${stage}.R1.fastq.gz"
    def read2_name = "${meta.read_group_id}.${stage}.R2.fastq.gz"

    """
    ln -s ${reads[0]} ${read1_name}
    ln -s ${reads[1]} ${read2_name}

    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${read1_name} \
        ${read2_name}
    """
}
