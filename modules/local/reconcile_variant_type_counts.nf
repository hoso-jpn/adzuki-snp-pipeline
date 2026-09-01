process RECONCILE_VARIANT_TYPE_COUNTS {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(raw_all_vcf), path(raw_all_vcf_index)
    tuple val(raw_snp_meta), path(raw_snp_vcf), path(raw_snp_vcf_index)
    tuple val(raw_indel_meta), path(raw_indel_vcf), path(raw_indel_vcf_index)
    path(raw_all_variant_qc)

    output:
    tuple(
        val(meta),
        path("${meta.id}.variant_type_accounting.tsv"),
        path("${meta.id}.variant_type_accounting.summary.txt"),
        emit: accounting
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
    reconcile_variant_type_counts.py \
        --cohort-id '${meta.id}' \
        --raw-all-vcf ${raw_all_vcf} \
        --raw-snp-vcf ${raw_snp_vcf} \
        --raw-indel-vcf ${raw_indel_vcf} \
        --raw-all-variant-qc ${raw_all_variant_qc} \
        --output ${meta.id}.variant_type_accounting.tsv \
        --summary-output ${meta.id}.variant_type_accounting.summary.txt
    """
}
