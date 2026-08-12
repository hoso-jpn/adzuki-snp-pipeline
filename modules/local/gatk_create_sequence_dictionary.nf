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

    script:
    """
    gatk CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}
