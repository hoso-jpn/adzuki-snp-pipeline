// Issue #20: --sample-ploidy is passed explicitly here even though GATK
// 4.6.2.0 was empirically confirmed (against this exact pinned container,
// comparing omitted vs. explicit --sample-ploidy against an identical
// GenomicsDB workspace) to already derive the correct per-sample ploidy
// from the gVCF/GenomicsDB input on its own -- the two invocations
// produced byte-identical CHROM/POS/REF/ALT, sample order, GT, AC, AN,
// AF, and PL. This is therefore a contract-strengthening change, not a
// bug fix: it makes params.sample_ploidy's effect on this process
// explicit in its own command line, rather than relying on GATK's
// undocumented-to-this-pipeline inference behavior, which could change
// in a future GATK version without warning.
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
        --sample-ploidy ${params.sample_ploidy} \
        --output ${interval_meta.id}.raw.vcf.gz \
        --create-output-variant-index true
    """
}
