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

include {
    GATK_SELECTVARIANTS
} from '../modules/local/gatk_selectvariants'

include {
    GATK_VARIANTFILTRATION
} from '../modules/local/gatk_variantfiltration'

include {
    GATK_SELECTPASSVARIANTS
} from '../modules/local/gatk_selectpassvariants'

include {
    BCFTOOLS_STATS
} from '../modules/local/bcftools_stats'

include {
    SUMMARIZE_VARIANT_QC
} from '../modules/local/summarize_variant_qc'

include {
    SUMMARIZE_FILTER_QC
} from '../modules/local/summarize_filter_qc'

include {
    RECONCILE_VARIANT_TYPE_COUNTS
} from '../modules/local/reconcile_variant_type_counts'

include {
    GS_NORMALIZE_VARIANTS
} from '../modules/local/gs_normalize_variants'

include {
    CLASSIFY_NORMALIZED_VARIANTS
} from '../modules/local/classify_normalized_variants'

include {
    GS_INDEX_CLASSIFIED_VARIANTS
} from '../modules/local/gs_index_classified_variants'

include {
    GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_GS
} from '../modules/local/gatk_variantfiltration'

include {
    GATK_SELECTPASSVARIANTS as GATK_SELECTPASSVARIANTS_GS
} from '../modules/local/gatk_selectpassvariants'

// Hard-filter definitions shared by the primary SNP/indel filtering
// lineage and the GS-panel-specific re-filtering lineage (see
// GATK_VARIANTFILTRATION_GS below), so the configured thresholds are
// never duplicated between the two call sites.
def snpHardFilters() {
    return [
        [
            name: 'SNP_QD_LOW',
            expression: "QD < ${params.snp_filter_qd_min}",
        ],
        [
            name: 'SNP_QUAL_LOW',
            expression: "QUAL < ${params.snp_filter_qual_min}",
        ],
        [
            name: 'SNP_SOR_HIGH',
            expression: "SOR > ${params.snp_filter_sor_max}",
        ],
        [
            name: 'SNP_FS_HIGH',
            expression: "FS > ${params.snp_filter_fs_max}",
        ],
        [
            name: 'SNP_MQ_LOW',
            expression: "MQ < ${params.snp_filter_mq_min}",
        ],
        [
            name: 'SNP_MQRANKSUM_LOW',
            expression: "MQRankSum < ${params.snp_filter_mq_rank_sum_min}",
        ],
        [
            name: 'SNP_READPOSRANKSUM_LOW',
            expression: "ReadPosRankSum < ${params.snp_filter_read_pos_rank_sum_min}",
        ],
    ]
}

def indelHardFilters() {
    return [
        [
            name: 'INDEL_QD_LOW',
            expression: "QD < ${params.indel_filter_qd_min}",
        ],
        [
            name: 'INDEL_QUAL_LOW',
            expression: "QUAL < ${params.indel_filter_qual_min}",
        ],
        [
            name: 'INDEL_FS_HIGH',
            expression: "FS > ${params.indel_filter_fs_max}",
        ],
        [
            name: 'INDEL_READPOSRANKSUM_LOW',
            expression: "ReadPosRankSum < ${params.indel_filter_read_pos_rank_sum_min}",
        ],
    ]
}

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

    variant_types_ch = GATK_GATHERVCFS.out.vcf
        .flatMap { meta, vcf, vcf_index ->
            [
                tuple(
                    meta + [variant_type: 'snp'],
                    vcf,
                    vcf_index,
                    'SNP',
                ),
                tuple(
                    meta + [variant_type: 'indel'],
                    vcf,
                    vcf_index,
                    'INDEL',
                ),
            ]
        }

    GATK_SELECTVARIANTS(variant_types_ch)

    hard_filter_inputs_ch = GATK_SELECTVARIANTS.out.vcf
        .map { meta, vcf, vcf_index ->
            def filters

            if (meta['variant_type'] == 'snp') {
                filters = snpHardFilters()
            } else if (meta['variant_type'] == 'indel') {
                filters = indelHardFilters()
            } else {
                error(
                    "unsupported variant type: ${meta['variant_type']}"
                )
            }

            tuple(meta, vcf, vcf_index, filters)
        }

    GATK_VARIANTFILTRATION(hard_filter_inputs_ch)
    GATK_SELECTPASSVARIANTS(GATK_VARIANTFILTRATION.out.vcf)

    raw_all_qc_inputs_ch = GATK_GATHERVCFS.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(
                meta + [
                    qc_stage: 'raw',
                    variant_type: 'all',
                ],
                vcf,
                vcf_index,
            )
        }

    raw_by_type_qc_inputs_ch = GATK_SELECTVARIANTS.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(
                meta + [qc_stage: 'raw'],
                vcf,
                vcf_index,
            )
        }

    filtered_qc_inputs_ch = GATK_VARIANTFILTRATION.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(
                meta + [qc_stage: 'filtered'],
                vcf,
                vcf_index,
            )
        }

    pass_qc_inputs_ch = GATK_SELECTPASSVARIANTS.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(
                meta + [qc_stage: 'pass'],
                vcf,
                vcf_index,
            )
        }

    variant_qc_inputs_ch = raw_all_qc_inputs_ch
        .mix(raw_by_type_qc_inputs_ch)
        .mix(filtered_qc_inputs_ch)
        .mix(pass_qc_inputs_ch)

    BCFTOOLS_STATS(variant_qc_inputs_ch)
    SUMMARIZE_VARIANT_QC(BCFTOOLS_STATS.out.stats)

    SUMMARIZE_FILTER_QC(GATK_VARIANTFILTRATION.out.vcf)

    raw_snp_only_ch = GATK_SELECTVARIANTS.out.vcf
        .filter { meta, _vcf, _vcf_index -> meta['variant_type'] == 'snp' }

    raw_indel_only_ch = GATK_SELECTVARIANTS.out.vcf
        .filter { meta, _vcf, _vcf_index -> meta['variant_type'] == 'indel' }

    raw_all_variant_qc_ch = SUMMARIZE_VARIANT_QC.out.qc
        .filter { meta, _variant_qc_tsv, _sample_qc_tsv, _summary_txt ->
            meta['qc_stage'] == 'raw' && meta['variant_type'] == 'all'
        }
        .map { _meta, variant_qc_tsv, _sample_qc_tsv, _summary_txt -> variant_qc_tsv }

    RECONCILE_VARIANT_TYPE_COUNTS(
        GATK_GATHERVCFS.out.vcf,
        raw_snp_only_ch,
        raw_indel_only_ch,
        raw_all_variant_qc_ch,
    )

    // GS panel: normalize raw/all (the only stage that can still contain
    // a GATK-MIXED record, since GATK_SELECTVARIANTS above already
    // excludes MIXED from both cohort.snp.vcf.gz and
    // cohort.indel.vcf.gz), reclassify by post-split REF/ALT shape, and
    // re-run the existing SNP hard filters on a distinct 'cohort_gs'
    // lineage so the primary raw/filtered/pass outputs above are
    // untouched.
    gs_normalize_input_ch = GATK_GATHERVCFS.out.vcf
        .map { meta, vcf, vcf_index -> tuple(meta + [id: 'cohort_gs'], vcf, vcf_index) }

    GS_NORMALIZE_VARIANTS(
        gs_normalize_input_ch,
        reference_ch,
        reference_fai_ch,
    )

    CLASSIFY_NORMALIZED_VARIANTS(GS_NORMALIZE_VARIANTS.out.vcf)
    GS_INDEX_CLASSIFIED_VARIANTS(CLASSIFY_NORMALIZED_VARIANTS.out.vcf)

    gs_hard_filter_inputs_ch = GS_INDEX_CLASSIFIED_VARIANTS.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(meta + [variant_type: 'snp'], vcf, vcf_index, snpHardFilters())
        }

    GATK_VARIANTFILTRATION_GS(gs_hard_filter_inputs_ch)
    GATK_SELECTPASSVARIANTS_GS(GATK_VARIANTFILTRATION_GS.out.vcf)

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
    variant_type_vcfs = GATK_SELECTVARIANTS.out.vcf
    filtered_vcfs = GATK_VARIANTFILTRATION.out.vcf
    pass_vcfs = GATK_SELECTPASSVARIANTS.out.vcf
    bcftools_stats = BCFTOOLS_STATS.out.stats
    variant_qc = SUMMARIZE_VARIANT_QC.out.qc
    filter_qc = SUMMARIZE_FILTER_QC.out.qc
    variant_type_accounting = RECONCILE_VARIANT_TYPE_COUNTS.out.accounting
    gs_normalized_vcf = GS_NORMALIZE_VARIANTS.out.vcf
    gs_classified_vcf = GS_INDEX_CLASSIFIED_VARIANTS.out.vcf
    gs_filtered_vcf = GATK_VARIANTFILTRATION_GS.out.vcf
    gs_pass_vcf = GATK_SELECTPASSVARIANTS_GS.out.vcf
    flagstat = SAMTOOLS_QC.out.flagstat
    stats = SAMTOOLS_QC.out.stats
    idxstats = SAMTOOLS_QC.out.idxstats
    reference_fai = reference_fai_ch
    reference_dict = reference_dict_ch
    bwa_indexes = bwa_indexes_ch
}
