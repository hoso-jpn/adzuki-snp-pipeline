process SUMMARIZE_VARIANT_QC {
    tag "${meta.id}:${meta.qc_stage}:${meta.variant_type}"
    label 'process_low'

    // The "-slim" Python image does not include `ps`, which Nextflow
    // needs inside the task container to collect resource-usage
    // metrics whenever a run is traced (`-with-trace`, as nf-test
    // always does); the full image below includes it.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(bcftools_stats_tsv)

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.variant_qc.tsv"),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.sample_qc.tsv"),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.summary.txt"),
        emit: qc
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
    summarize_variant_qc.py \
        --bcftools-stats ${bcftools_stats_tsv} \
        --cohort-id '${meta.id}' \
        --stage '${meta.qc_stage}' \
        --variant-type '${meta.variant_type}' \
        --variant-qc-output ${prefix}.variant_qc.tsv \
        --sample-qc-output ${prefix}.sample_qc.tsv \
        --summary-output ${prefix}.summary.txt
    """
}
