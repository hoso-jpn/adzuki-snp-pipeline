process GATK_GENOTYPEGVCFS {
    tag "${interval_meta.id}"
    label 'process_high'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple(
        val(interval_meta),
        val(interval),
        path(genomicsdb)
    )
    tuple val(reference_meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dict)

    output:
    tuple(
        val(interval_meta),
        path("${interval_meta.id}.raw.vcf.gz"),
        path("${interval_meta.id}.raw.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    """
    gatk --java-options "-Xmx${memory_gb}g" GenotypeGVCFs \
        --reference ${fasta} \
        --variant gendb://${genomicsdb} \
        --intervals '${interval}' \
        --output ${interval_meta.id}.raw.vcf.gz \
        --create-output-variant-index true
    """
}
