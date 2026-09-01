process GATK_CREATE_SEQUENCE_DICTIONARY {
    tag "${meta.id}"
    label 'process_low'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    publishDir(
        "${params.outdir}/reference",
        mode: 'copy',
        pattern: '*.dict'
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}.dict"), emit: dict

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
    gatk CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}
