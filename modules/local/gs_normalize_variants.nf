process GS_NORMALIZE_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(reference_meta), path(fasta)
    tuple val(fai_meta), path(fai)

    output:
    tuple(
        val(meta),
        path("${meta.id}.normalized.vcf.gz"),
        path("${meta.id}.normalized.vcf.gz.tbi"),
        emit: vcf
    )
    path("${meta.id}.normalize.report.txt"), emit: report

    script:
    """
    bcftools norm \
        --fasta-ref ${fasta} \
        --multiallelics -both \
        --check-ref e \
        --output-type z \
        --output ${meta.id}.normalized.vcf.gz \
        ${vcf} \
        2> ${meta.id}.normalize.report.txt

    bcftools index --tbi ${meta.id}.normalized.vcf.gz
    """
}
