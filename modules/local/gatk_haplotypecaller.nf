// Issue #30: dedicated resource label, independent of process_high
// (which GATK_GENOTYPEGVCFS also uses) -- see nextflow.config for the
// real 8-vs-4-cpu benchmark this cpus value is based on.
process GATK_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_haplotypecaller'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(reference_meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dict)

    output:
    tuple(
        val(meta),
        path("${meta.id}.g.vcf.gz"),
        path("${meta.id}.g.vcf.gz.tbi"),
        emit: gvcf
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    """
    gatk --java-options "-Xmx${memory_gb}g" HaplotypeCaller \
        --reference ${fasta} \
        --input ${bam} \
        --output ${meta.id}.g.vcf.gz \
        --emit-ref-confidence GVCF \
        --sample-ploidy ${params.sample_ploidy} \
        --native-pair-hmm-threads ${task.cpus} \
        --create-output-variant-index true
    """
}
