include { BWA_MEM2_INDEX } from '../../../modules/local/bwa_mem2_index'
include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools_faidx'
include { GATK_CREATE_SEQUENCE_DICTIONARY } from '../../../modules/local/gatk_create_sequence_dictionary'
include { VALIDATE_REFERENCE_BUNDLE } from '../../../modules/local/validate_reference_bundle'

workflow BUILD_REFERENCE_FIXTURE {
    take:
    reference

    main:
    BWA_MEM2_INDEX(reference)
    SAMTOOLS_FAIDX(reference)
    GATK_CREATE_SEQUENCE_DICTIONARY(reference)
    origin = BWA_MEM2_INDEX.out.build_receipt.map { _meta, receipt, version, container ->
        tuple([mode: 'generated', container: container], receipt, version)
    }
    VALIDATE_REFERENCE_BUNDLE(reference, SAMTOOLS_FAIDX.out.fai,
        GATK_CREATE_SEQUENCE_DICTIONARY.out.dict, BWA_MEM2_INDEX.out.indexes, origin)

    emit:
    fasta = reference
    fai = VALIDATE_REFERENCE_BUNDLE.out.fai
    dict = VALIDATE_REFERENCE_BUNDLE.out.dict
    indexes = VALIDATE_REFERENCE_BUNDLE.out.indexes
    manifest = VALIDATE_REFERENCE_BUNDLE.out.manifest
}
