process BWA_MEM2_INDEX {
    tag "${meta.id}"
    label 'process_bwa_index'

    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0@sha256:374f4b910b4b04d32772fccd4cab1cdcd0758356856a960cfa8b1edebfd38c9f'

    publishDir(
        "${params.outdir}/reference",
        mode: 'copy',
        pattern: '*.{0123,amb,ann,bwt.2bit.64,pac}'
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.name}.*"), emit: indexes

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
    bwa-mem2 index ${fasta}
    """
}
