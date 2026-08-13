include {
    FASTQC as FASTQC_RAW
} from '../modules/local/fastqc'

include {
    FASTQC as FASTQC_TRIMMED
} from '../modules/local/fastqc'

include {
    FASTP
} from '../modules/local/fastp'

include {
    SAMTOOLS_FAIDX
} from '../modules/local/samtools_faidx'

include {
    GATK_CREATE_SEQUENCE_DICTIONARY
} from '../modules/local/gatk_create_sequence_dictionary'

include {
    BWA_MEM2_INDEX
} from '../modules/local/bwa_mem2_index'

include {
    BWA_MEM2_MEM
} from '../modules/local/bwa_mem2_mem'

include {
    SAMTOOLS_SORT
} from '../modules/local/samtools_sort'

include {
    SAMTOOLS_MERGE
} from '../modules/local/samtools_merge'

include {
    GATK_MARKDUPLICATES
} from '../modules/local/gatk_markduplicates'

include {
    SAMTOOLS_INDEX
} from '../modules/local/samtools_index'

include {
    SAMTOOLS_QC
} from '../modules/local/samtools_qc'

include {
    GATK_HAPLOTYPECALLER
} from '../modules/local/gatk_haplotypecaller'

include {
    GATK_GENOMICSDBIMPORT
} from '../modules/local/gatk_genomicsdbimport'

include {
    GATK_GENOTYPEGVCFS
} from '../modules/local/gatk_genotypegvcfs'

include {
    GATK_GATHERVCFS
} from '../modules/local/gatk_gathervcfs'

workflow ADZUKI_SNP_PIPELINE {
    take:
    samples_ch
    reference_ch

    main:
    raw_reads_ch = samples_ch.map {
        meta,
        read1,
        read2 ->
        tuple(meta, [read1, read2])
    }

    FASTQC_RAW(
        raw_reads_ch,
        'raw'
    )

    FASTP(samples_ch)

    trimmed_reads_for_qc_ch = FASTP.out.reads.map {
        meta,
        read1,
        read2 ->
        tuple(meta, [read1, read2])
    }

    FASTQC_TRIMMED(
        trimmed_reads_for_qc_ch,
        'trimmed'
    )

    if (params.reference_fai) {
        reference_fai_ch = reference_ch.map {
            meta,
            _fasta ->
            tuple(
                meta,
                file(
                    params.reference_fai,
                    checkIfExists: true
                )
            )
        }
    } else {
        SAMTOOLS_FAIDX(reference_ch)
        reference_fai_ch = SAMTOOLS_FAIDX.out.fai
    }

    if (params.reference_dict) {
        reference_dict_ch = reference_ch.map {
            meta,
            _fasta ->
            tuple(
                meta,
                file(
                    params.reference_dict,
                    checkIfExists: true
                )
            )
        }
    } else {
        GATK_CREATE_SEQUENCE_DICTIONARY(reference_ch)
        reference_dict_ch =
            GATK_CREATE_SEQUENCE_DICTIONARY.out.dict
    }

    if (params.bwa_index_prefix) {
        bwa_indexes_ch = reference_ch.map {
            meta,
            _fasta ->
            def indexes = [
                file(
                    "${params.bwa_index_prefix}.0123",
                    checkIfExists: true
                ),
                file(
                    "${params.bwa_index_prefix}.amb",
                    checkIfExists: true
                ),
                file(
                    "${params.bwa_index_prefix}.ann",
                    checkIfExists: true
                ),
                file(
                    "${params.bwa_index_prefix}.bwt.2bit.64",
                    checkIfExists: true
                ),
                file(
                    "${params.bwa_index_prefix}.pac",
                    checkIfExists: true
                )
            ]

            tuple(meta, indexes)
        }
    } else {
        BWA_MEM2_INDEX(reference_ch)
        bwa_indexes_ch = BWA_MEM2_INDEX.out.indexes
    }

    BWA_MEM2_MEM(
        FASTP.out.reads,
        reference_ch,
        bwa_indexes_ch
    )

    SAMTOOLS_SORT(BWA_MEM2_MEM.out.sam)

    sample_bams_ch = SAMTOOLS_SORT.out.bam
        .map {
            meta,
            bam ->
            tuple(meta.id, bam)
        }
        .groupTuple()
        .map {
            sample_id,
            bams ->
            def sorted_bams = bams.sort {
                bam ->
                bam.name
            }

            tuple(
                [id: sample_id],
                sorted_bams
            )
        }

    SAMTOOLS_MERGE(sample_bams_ch)
    GATK_MARKDUPLICATES(SAMTOOLS_MERGE.out.bam)
    SAMTOOLS_INDEX(GATK_MARKDUPLICATES.out.bam)
    SAMTOOLS_QC(SAMTOOLS_INDEX.out.bam)

    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out.bam,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch,
    )

    intervals_ch = reference_fai_ch.flatMap {
        _reference_meta,
        fai ->
        fai.readLines()
            .findAll { line -> !line.isBlank() }
            .withIndex()
            .collect { line, index ->
                def fields = line.split('\\t')

                if (fields.size() < 2) {
                    error("Invalid FASTA index entry: ${line}")
                }

                def contig = fields[0]
                def safe_contig = contig.replaceAll(
                    '[^A-Za-z0-9._-]',
                    '_',
                )

                tuple(
                    [
                        id: String.format(
                            'interval_%06d_%s',
                            index + 1,
                            safe_contig,
                        ),
                        rank: index,
                        contig: contig,
                    ],
                    contig,
                )
            }
    }

    gvcfs_ch = GATK_HAPLOTYPECALLER.out.gvcf
        .map { _meta, gvcf, _gvcf_index -> gvcf }
        .collect()

    gvcf_indexes_ch = GATK_HAPLOTYPECALLER.out.gvcf
        .map { _meta, _gvcf, gvcf_index -> gvcf_index }
        .collect()

    GATK_GENOMICSDBIMPORT(
        intervals_ch,
        gvcfs_ch,
        gvcf_indexes_ch,
    )

    GATK_GENOTYPEGVCFS(
        GATK_GENOMICSDBIMPORT.out.genomicsdb,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch,
    )

    cohort_vcfs_ch = GATK_GENOTYPEGVCFS.out.vcf
        .collect(flat: false)
        .map { entries ->
            def sorted_entries = entries.sort {
                entry -> entry[1].getFileName().toString()
            }

            tuple(
                [id: 'cohort'],
                sorted_entries.collect { entry -> entry[1] },
                sorted_entries.collect { entry -> entry[2] },
            )
        }

    GATK_GATHERVCFS(cohort_vcfs_ch)

    emit:
    raw_fastqc_html = FASTQC_RAW.out.html
    raw_fastqc_zip = FASTQC_RAW.out.zip
    trimmed_reads = FASTP.out.reads
    fastp_reports = FASTP.out.reports
    trimmed_fastqc_html = FASTQC_TRIMMED.out.html
    trimmed_fastqc_zip = FASTQC_TRIMMED.out.zip
    mapping_logs = BWA_MEM2_MEM.out.log
    duplicate_metrics = GATK_MARKDUPLICATES.out.metrics
    marked_bams = SAMTOOLS_INDEX.out.bam
    sample_gvcfs = GATK_HAPLOTYPECALLER.out.gvcf
    genomicsdb = GATK_GENOMICSDBIMPORT.out.genomicsdb
    interval_vcfs = GATK_GENOTYPEGVCFS.out.vcf
    raw_vcf = GATK_GATHERVCFS.out.vcf
    flagstat = SAMTOOLS_QC.out.flagstat
    stats = SAMTOOLS_QC.out.stats
    idxstats = SAMTOOLS_QC.out.idxstats
    reference_fai = reference_fai_ch
    reference_dict = reference_dict_ch
    bwa_indexes = bwa_indexes_ch
}
