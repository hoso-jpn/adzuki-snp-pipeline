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

    // Issue #42: this process's *effective* container -- Nextflow's own
    // task.container, resolved after any withName/alias/fully-qualified-
    // selector/profile override on top of the `container` directive above.
    // This module is included twice, as FASTQC_RAW and FASTQC_TRIMMED, and
    // a withName selector matches the *included* name, so the two aliases
    // can be overridden independently; the run manifest records them as
    // two separate process identities rather than one "fastqc" entry. See
    // workflows/adzuki_snp_pipeline.nf for the canonical process keys and
    // docs/run_manifest_data_contract.md for the schema v2 contract.
    val(task.container), emit: container_id

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
