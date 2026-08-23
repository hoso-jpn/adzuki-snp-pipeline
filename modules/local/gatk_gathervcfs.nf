// Issue #11: List-safety fix only (see GATK_GENOMICSDBIMPORT for the full
// rationale -- the same real-data run exposed the identical bug pattern
// here: Nextflow silently unwraps a single-element List into a bare
// scalar Path/File for `path(vcfs)`/`path(vcf_indexes)`, and `.size()` on
// that scalar returns a byte count, not "1"). This process's Xmx formula
// and 'process_medium' label are unchanged -- out of this Issue's scope,
// which is GenomicsDBImport-specific memory hardening.
process GATK_GATHERVCFS {
    tag "${meta.id}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple(
        val(meta),
        path(vcfs),
        path(vcf_indexes)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.raw.vcf.gz"),
        path("${meta.id}.raw.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    vcf_list = vcfs instanceof List ? vcfs : [vcfs]
    vcf_index_list = vcf_indexes instanceof List ? vcf_indexes : [vcf_indexes]

    if (vcf_list.size() != vcf_index_list.size()) {
        error(
            'the number of VCFs and indexes must match: ' +
            "${vcf_list.size()} VCFs, " +
            "${vcf_index_list.size()} indexes"
        )
    }

    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )
    input_arguments = vcf_list
        .collect { vcf -> "--INPUT ${vcf}" }
        .join(" \\\n        ")

    """
    gatk --java-options "-Xmx${memory_gb}g" GatherVcfs \
        ${input_arguments} \
        --OUTPUT ${meta.id}.raw.vcf.gz

    gatk --java-options "-Xmx${memory_gb}g" IndexFeatureFile \
        --input ${meta.id}.raw.vcf.gz
    """
}
