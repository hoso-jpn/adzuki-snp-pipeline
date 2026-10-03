include { GATK_CREATE_SEQUENCE_DICTIONARY } from '../../../modules/local/gatk_create_sequence_dictionary'
include { BUILD_SAMPLE_NAME_MAP } from '../../../modules/local/build_sample_name_map'
include { GATK_GENOMICSDBIMPORT } from '../../../modules/local/gatk_genomicsdbimport'
include { GATK_GENOTYPEGVCFS } from '../../../modules/local/gatk_genotypegvcfs'

process MAKE_BOUNDARY_GVCF {
    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'
    label 'process_low'
    input:
    tuple val(meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dictionary)
    path(source)
    output:
    path('boundary.g.vcf.gz'), emit: gvcf
    path('boundary.g.vcf.gz.tbi'), emit: index
    script:
    """
    gatk --java-options '-Xmx1g' SelectVariants -R ${fasta} -V ${source} \
        -O boundary.g.vcf.gz --create-output-variant-index true
    """
}

workflow JOINT_WINDOW_BOUNDARY {
    take:
    reference
    fai
    source
    main:
    GATK_CREATE_SEQUENCE_DICTIONARY(reference)
    MAKE_BOUNDARY_GVCF(reference, fai, GATK_CREATE_SEQUENCE_DICTIONARY.out.dict, source)
    BUILD_SAMPLE_NAME_MAP(
        [[sample_id: 'boundary', gvcf: 'boundary.g.vcf.gz', index: 'boundary.g.vcf.gz.tbi']],
        MAKE_BOUNDARY_GVCF.out.gvcf, MAKE_BOUNDARY_GVCF.out.index)
    intervals = channel.of(
        tuple([id: 'interval_000001', rank: 0, contig: 'chrSynthetic1'], 'chrSynthetic1:1-2000'),
        tuple([id: 'interval_000002', rank: 1, contig: 'chrSynthetic1'], 'chrSynthetic1:2001-5000'),
    )
    GATK_GENOMICSDBIMPORT(intervals, BUILD_SAMPLE_NAME_MAP.out.gvcfs.first(),
        BUILD_SAMPLE_NAME_MAP.out.indexes.first(), BUILD_SAMPLE_NAME_MAP.out.sample_map.first())
    GATK_GENOTYPEGVCFS(GATK_GENOMICSDBIMPORT.out.genomicsdb, reference, fai,
        GATK_CREATE_SEQUENCE_DICTIONARY.out.dict)
    emit:
    vcf = GATK_GENOTYPEGVCFS.out.vcf
}
