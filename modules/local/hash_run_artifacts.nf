process HASH_RUN_ARTIFACTS {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

    input:
    tuple val(meta), path(artifacts)

    output:
    tuple(
        val(meta),
        path("${meta.id}.artifact_checksums.tsv"),
        emit: checksums
    )
    // Issue #42: see modules/local/hash_input_fastqs.nf for why every
    // containerized process emits its own effective task.container.
    val(task.container), emit: container_id

    script:
    // Issue #42: invoked once per artifact *group* (per sample for
    // gVCFs, once for the cohort-level VCFs) rather than once for the
    // whole run, so hashing fans out with the rest of the pipeline
    // instead of funnelling every deliverable through the single
    // terminal task that writes the manifest. At this pipeline's target
    // scale that is the difference between hashing a few files per task
    // and re-staging hundreds of gigabytes of gVCFs into one task.
    def artifact_args = artifacts
        .collect { artifact -> "--artifact ${artifact}" }
        .join(' \\\n        ')

    """
    hash_run_artifacts.py \
        ${artifact_args} \
        --output ${meta.id}.artifact_checksums.tsv
    """
}
