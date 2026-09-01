process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1@sha256:a130447589651ed09252aa95a5e4f4132942cdb54d835d81a04a9a930d656561'

    publishDir(
        "${params.outdir}/alignment",
        mode: 'copy',
        pattern: '*.markdup.bam*'
    )

    input:
    tuple val(meta), path(bam)

    output:
    tuple(
        val(meta),
        path(bam),
        path("${bam}.bai"),
        emit: bam
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
    samtools index \
        -@ ${task.cpus} \
        ${bam}
    """
}
