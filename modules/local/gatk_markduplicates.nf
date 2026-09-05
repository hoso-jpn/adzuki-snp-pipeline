process GATK_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    publishDir(
        "${params.outdir}/qc/markduplicates",
        mode: 'copy',
        pattern: '*.metrics.txt'
    )

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"), emit: bam
    tuple(
        val(meta),
        path("${meta.id}.markduplicates.metrics.txt"),
        emit: metrics
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
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    """
    gatk --java-options "-Xmx${memory_gb}g" MarkDuplicates \
        --INPUT ${bam} \
        --OUTPUT ${meta.id}.markdup.bam \
        --METRICS_FILE ${meta.id}.markduplicates.metrics.txt \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX false \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.optical_duplicate_pixel_distance} \
        --VALIDATION_STRINGENCY STRICT
    """
}
