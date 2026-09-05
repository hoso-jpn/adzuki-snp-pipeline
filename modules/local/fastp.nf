process FASTP {
    tag "${meta.id}:${meta.read_group_id}"
    label 'process_low'

    container 'quay.io/biocontainers/fastp:1.3.6--h43da1c4_0@sha256:cbbe2402b6b6704df470d7d77dcb498eefd5bcd01f4c38be0ec69899e79ac134'

    publishDir(
        "${params.outdir}/reads/trimmed",
        mode: 'copy',
        pattern: '*.fastq.gz'
    )
    publishDir(
        "${params.outdir}/qc/fastp",
        mode: 'copy',
        pattern: '*.{json,html}'
    )

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple(
        val(meta),
        path("${meta.read_group_id}.trimmed_R1.fastq.gz"),
        path("${meta.read_group_id}.trimmed_R2.fastq.gz"),
        emit: reads
    )
    tuple(
        val(meta),
        path("${meta.read_group_id}.fastp.json"),
        path("${meta.read_group_id}.fastp.html"),
        emit: reports
    )

    // Issue #42: this process's *effective* container -- Nextflow's own
    // task.container, resolved after any withName/alias/fully-qualified-
    // selector/profile override on top of the `container` directive above
    // -- so the run-level provenance manifest records what this task
    // actually ran in rather than a default the pipeline assumed. See
    // workflows/adzuki_snp_pipeline.nf for the canonical process key this
    // invocation is recorded under, and docs/run_manifest_data_contract.md
    // for the schema v2 contract.
    val(task.container), emit: container_id

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${meta.read_group_id}.trimmed_R1.fastq.gz \
        --out2 ${meta.read_group_id}.trimmed_R2.fastq.gz \
        --json ${meta.read_group_id}.fastp.json \
        --html ${meta.read_group_id}.fastp.html \
        --detect_adapter_for_pe \
        --thread ${task.cpus}
    """
}
