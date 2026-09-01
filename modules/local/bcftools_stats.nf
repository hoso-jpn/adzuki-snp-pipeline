process BCFTOOLS_STATS {
    tag "${meta.id}:${meta.qc_stage}:${meta.variant_type}"
    label 'process_low'

    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.bcftools.stats.tsv"),
        emit: stats
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
    prefix = "${meta.id}.${meta.qc_stage}.${meta.variant_type}"

    """
    bcftools stats \
        --threads ${task.cpus} \
        --samples - \
        ${vcf} \
        > ${prefix}.bcftools.stats.tsv
    """
}
