process GS_NORMALIZE_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    container params.containers.bcftools

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
    // Issue #52: this process's *effective* container -- Nextflow's own
    // task.container, resolved after any withName/alias/fully-qualified-
    // selector/profile override on top of the `container` directive above
    // -- so a consumer (BUILD_GS_PANEL_MANIFEST) can record what actually
    // ran rather than trusting the directive's default value.
    val(task.container), emit: container_id

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
