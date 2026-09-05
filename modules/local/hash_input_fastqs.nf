process HASH_INPUT_FASTQS {
    tag "${meta.read_group_id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required. Sourced from
    // conf/containers.config like every other process this repository has
    // moved there (Issue #52), rather than re-declaring the digest here.
    container params.containers.python

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple(
        val(meta),
        path("${meta.read_group_id}.input_provenance.tsv"),
        emit: provenance
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
    // Issue #42: `meta.rank` is the read group's zero-based position in
    // the samplesheet, attached in main.nf. Nextflow gives no ordering
    // guarantee across these parallel tasks, so the rank travels with the
    // row and bin/build_run_manifest.py sorts on it -- then drops it --
    // to present the manifest's read groups in samplesheet order. It is
    // never published.
    """
    hash_input_fastqs.py \
        --rank ${meta.rank} \
        --sample-id '${meta.sample_id}' \
        --read-group-id '${meta.read_group_id}' \
        --library-id '${meta.library_id}' \
        --platform '${meta.platform}' \
        --platform-unit '${meta.platform_unit}' \
        --fastq-1 ${read1} \
        --fastq-2 ${read2} \
        --output ${meta.read_group_id}.input_provenance.tsv
    """
}
