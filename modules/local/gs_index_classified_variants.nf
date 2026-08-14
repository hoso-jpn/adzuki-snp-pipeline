process GS_INDEX_CLASSIFIED_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'

    input:
    tuple val(meta), path(classified_vcf)

    output:
    tuple(
        val(meta),
        path("${meta.id}.classified.vcf.gz"),
        path("${meta.id}.classified.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    """
    bcftools view --output-type z --output ${meta.id}.classified.vcf.gz ${classified_vcf}
    bcftools index --tbi ${meta.id}.classified.vcf.gz
    """
}
