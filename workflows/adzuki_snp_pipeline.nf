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
    VALIDATE_REFERENCE_CONTIGS
} from '../modules/local/validate_reference_contigs'

include {
    BWA_MEM2_INDEX
} from '../modules/local/bwa_mem2_index'

include {
    BWA_MEM2_MEM_SORT
} from '../modules/local/bwa_mem2_mem_sort'

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
    MULTIQC
} from '../modules/local/multiqc'

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

include {
    BUILD_GS_PANEL
} from '../modules/local/build_gs_panel'

include {
    RECONCILE_GS_PANEL_ACCOUNTING
} from '../modules/local/reconcile_gs_panel_accounting'

include {
    BUILD_GS_PANEL_MANIFEST
} from '../modules/local/build_gs_panel_manifest'

include {
    HASH_INPUT_FASTQS
} from '../modules/local/hash_input_fastqs'

include {
    HASH_REFERENCE_BUNDLE
} from '../modules/local/hash_reference_bundle'

include {
    HASH_RUN_ARTIFACTS
} from '../modules/local/hash_run_artifacts'

include {
    BUILD_RUN_MANIFEST
} from '../modules/local/build_run_manifest'

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

// Issue #17: the single place that maps a canonical variant type to its
// hard-filter list, used by both the primary SNP/indel filtering
// lineage and the GS lineage below. This is a total lookup over the two
// representable variant types rather than an if/else-if with a third
// `else error("unsupported variant type: ...")` branch, because
// GATK_SELECTVARIANTS already fail-fasts at its own process boundary on
// any meta.variant_type outside {snp, indel} (see
// modules/local/gatk_selectvariants.nf) -- and every tuple reaching the
// primary call site is that process's own output, so the branch was
// unreachable and duplicated a validation that now has one owner. A
// lookup miss would still be caught, since GATK_VARIANTFILTRATION
// rejects a null/empty `filters` at its own boundary.
def hardFiltersForVariantType(variant_type) {
    def filters_by_variant_type = [
        snp: snpHardFilters(),
        indel: indelHardFilters(),
    ]

    return filters_by_variant_type[variant_type?.toString()]
}

// Issue #42: label one process invocation's effective container
// identities with the canonical process key the run manifest records
// them under.
//
// The key is assigned here, explicitly, rather than derived at run time
// from something like `task.process`. Two of this pipeline's modules are
// included twice under different aliases (FASTQC as FASTQC_RAW and
// FASTQC_TRIMMED; GATK_VARIANTFILTRATION and GATK_SELECTPASSVARIANTS
// each also as their `_GS` alias), and each alias can be overridden
// independently -- so a key derived from the module's own process name
// would collapse pairs that are allowed to differ, and one derived from
// a runtime-formatted qualified name would tie a published schema to how
// a given Nextflow version happens to spell it. The keys below are part
// of the schema v2 contract (docs/run_manifest_data_contract.md) and
// change only when this workflow deliberately changes them.
def containerProvenance(container_id_ch, String process_key) {
    return container_id_ch.map { container -> "${process_key}\t${container}" }
}

workflow ADZUKI_SNP_PIPELINE {
    take:
    samples_ch
    reference_ch
    read_group_counts_by_sample
    // Issue #42: the same samplesheet rows as samples_ch, each carrying
    // its zero-based samplesheet position, built in main.nf so this
    // provenance branch cannot alter samples_ch's own cardinality or
    // ordering contract.
    input_provenance_rows_ch
    // Issue #42: the full 40-character commit this run executed from,
    // resolved in main.nf (see resolvePipelineCommit there for why
    // workflow.commitId alone is not enough).
    run_git_commit

    main:
    // Issue #42: container identities for processes that only run under
    // some configurations -- reference indexing when no prebuilt bundle
    // was supplied, and the whole GS lineage when it is enabled. They are
    // accumulated in the branches that actually invoke them, so the run
    // manifest lists the processes this run *executed* rather than every
    // process it could have. Recording a default for a process that never
    // ran would be a claim about software that never touched the data.
    optional_container_provenance_ch = channel.empty()

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
        optional_container_provenance_ch = optional_container_provenance_ch.mix(
            containerProvenance(SAMTOOLS_FAIDX.out.container_id, 'samtools_faidx')
        )
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
        optional_container_provenance_ch = optional_container_provenance_ch.mix(
            containerProvenance(
                GATK_CREATE_SEQUENCE_DICTIONARY.out.container_id,
                'gatk_create_sequence_dictionary',
            )
        )
    }

    // Issue #11 / #41: validate once here for both the generated and
    // prebuilt reference-bundle paths. The pass-through FAI/dict outputs
    // gate variant calling directly. In addition, reference_validation_gate_ch
    // below gates the FASTP reads before sample-dependent mapping, so source
    // order is not being mistaken for an execution barrier. Reference-only
    // BWA index construction may still proceed in parallel because it consumes
    // the FASTA, not the FAI/dict pair.
    VALIDATE_REFERENCE_CONTIGS(
        reference_fai_ch,
        reference_dict_ch,
    )
    reference_fai_ch = VALIDATE_REFERENCE_CONTIGS.out.fai
    reference_dict_ch = VALIDATE_REFERENCE_CONTIGS.out.dict

    // `.first()` turns the validator's single successful emission into a
    // reusable value channel. Every FASTP tuple is combined with that value,
    // preserving read-group cardinality while making validation success an
    // explicit dataflow prerequisite for every BWA_MEM2_MEM_SORT task. A
    // validation failure emits no value, so mapping cannot be scheduled.
    reference_validation_gate_ch = reference_fai_ch
        .map { meta, _fai -> meta.id }
        .first()
    validated_mapping_reads_ch = FASTP.out.reads
        .combine(reference_validation_gate_ch)
        .map {
            meta,
            read1,
            read2,
            validated_reference_id ->
            if (meta == null || validated_reference_id == null) {
                error('reference validation gate received an invalid metadata value')
            }
            tuple(meta, read1, read2)
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
        optional_container_provenance_ch = optional_container_provenance_ch.mix(
            containerProvenance(BWA_MEM2_INDEX.out.container_id, 'bwa_mem2_index')
        )
    }

    BWA_MEM2_MEM_SORT(
        validated_mapping_reads_ch,
        reference_ch,
        bwa_indexes_ch
    )

    // Issue #8: groupTuple() with no explicit size buffers every
    // tuple until BWA_MEM2_MEM_SORT's *entire* upstream channel
    // closes -- i.e. until every read group of every sample has
    // finished mapping -- before emitting even a single completed
    // sample group, even if that sample's own read groups all
    // finished long ago. read_group_counts_by_sample (computed once
    // in main.nf, directly from the samplesheet, before any channel
    // work begins) tells groupKey() exactly how many tuples to expect
    // per sample_id, so groupTuple() can emit each sample's group the
    // moment its own read groups are all in, independent of any other
    // sample's progress. `groupKey.target` recovers the original
    // sample_id string afterwards (see Nextflow's grouping docs for
    // this pattern) -- the sorted-bam-name ordering contract below is
    // otherwise unchanged.
    sample_bams_ch = BWA_MEM2_MEM_SORT.out.bam
        .map {
            meta,
            bam ->
            def expected_read_groups = read_group_counts_by_sample[meta.id]
            tuple(groupKey(meta.id, expected_read_groups), bam)
        }
        .groupTuple()
        .map {
            sample_key,
            bams ->
            def sorted_bams = bams.sort {
                bam ->
                bam.name
            }

            tuple(
                [id: sample_key.target],
                sorted_bams
            )
        }

    SAMTOOLS_MERGE(sample_bams_ch)
    GATK_MARKDUPLICATES(SAMTOOLS_MERGE.out.bam)
    SAMTOOLS_INDEX(GATK_MARKDUPLICATES.out.bam)
    SAMTOOLS_QC(SAMTOOLS_INDEX.out.bam)

    raw_fastqc_for_multiqc_ch = FASTQC_RAW.out.zip
        .map { _meta, report -> report }
        .collect()
    trimmed_fastqc_for_multiqc_ch = FASTQC_TRIMMED.out.zip
        .map { _meta, report -> report }
        .collect()
    fastp_for_multiqc_ch = FASTP.out.reports
        .map { _meta, json, _html -> json }
        .collect()
    markduplicates_for_multiqc_ch = GATK_MARKDUPLICATES.out.metrics
        .map { _meta, report -> report }
        .collect()
    flagstat_for_multiqc_ch = SAMTOOLS_QC.out.flagstat
        .map { _meta, report -> report }
        .collect()
    stats_for_multiqc_ch = SAMTOOLS_QC.out.stats
        .map { _meta, report -> report }
        .collect()
    idxstats_for_multiqc_ch = SAMTOOLS_QC.out.idxstats
        .map { _meta, report -> report }
        .collect()
    multiqc_config_ch = channel.value(
        file("${projectDir}/conf/multiqc_config.yaml", checkIfExists: true)
    )

    MULTIQC(
        raw_fastqc_for_multiqc_ch,
        trimmed_fastqc_for_multiqc_ch,
        fastp_for_multiqc_ch,
        markduplicates_for_multiqc_ch,
        flagstat_for_multiqc_ch,
        stats_for_multiqc_ch,
        idxstats_for_multiqc_ch,
        multiqc_config_ch,
    )

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

    // Issue #17: `meta.variant_type` alone -- the uppercase 'SNP'/'INDEL'
    // GATK argument that used to be passed alongside it as a fourth
    // tuple element is now derived from this same canonical tag inside
    // modules/local/gatk_selectvariants.nf, so the two can no longer
    // disagree (e.g. [variant_type: 'snp'] paired with 'INDEL').
    variant_types_ch = GATK_GATHERVCFS.out.vcf
        .flatMap { meta, vcf, vcf_index ->
            [
                tuple(
                    meta + [variant_type: 'snp'],
                    vcf,
                    vcf_index,
                ),
                tuple(
                    meta + [variant_type: 'indel'],
                    vcf,
                    vcf_index,
                ),
            ]
        }

    GATK_SELECTVARIANTS(variant_types_ch)

    hard_filter_inputs_ch = GATK_SELECTVARIANTS.out.vcf
        .map { meta, vcf, vcf_index ->
            tuple(
                meta,
                vcf,
                vcf_index,
                hardFiltersForVariantType(meta['variant_type']),
            )
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

    // See main.nf for why this is `.toString().toBoolean()` rather than
    // a bare truthiness check: a CLI-provided `--enable_gs_panel false`
    // resolves to the String "false", which is truthy in Groovy.
    gs_panel_enabled = params.enable_gs_panel.toString().toBoolean()

    if (gs_panel_enabled) {
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
                tuple(
                    meta + [variant_type: 'snp'],
                    vcf,
                    vcf_index,
                    hardFiltersForVariantType('snp'),
                )
            }

        GATK_VARIANTFILTRATION_GS(gs_hard_filter_inputs_ch)
        GATK_SELECTPASSVARIANTS_GS(GATK_VARIANTFILTRATION_GS.out.vcf)

        // Revert meta.id from 'cohort_gs' back to 'cohort' now that the GS
        // lineage's own provenance trail (normalized/classified/filtered)
        // is complete, so the final panel deliverables use clean filenames
        // (cohort.gs_panel.*) rather than a redundant cohort_gs.gs_panel.* .
        gs_pass_for_panel_ch = GATK_SELECTPASSVARIANTS_GS.out.vcf
            .map { meta, vcf, vcf_index -> tuple(meta + [id: 'cohort'], vcf, vcf_index) }

        BUILD_GS_PANEL(gs_pass_for_panel_ch)

        RECONCILE_GS_PANEL_ACCOUNTING(
            GATK_GATHERVCFS.out.vcf,
            GS_NORMALIZE_VARIANTS.out.vcf,
            CLASSIFY_NORMALIZED_VARIANTS.out.accounting,
            gs_pass_for_panel_ch,
            BUILD_GS_PANEL.out.matrix,
            BUILD_GS_PANEL.out.variant_metadata,
            BUILD_GS_PANEL.out.sample_metadata,
        )

        // Issue #52: no literal container digests here any more. Each GS
        // process below emits its own `container_id` -- task.container,
        // captured from inside that exact task after Nextflow has already
        // resolved any withName/alias/fully-qualified-selector/profile
        // override on top of its `container` directive (itself sourced
        // from conf/containers.config, the single default-value source of
        // truth; see that file's header for why the default alone is not
        // sufficient). BUILD_GS_PANEL_MANIFEST therefore records what each
        // task actually ran in, not what the pipeline assumed it would.
        //
        // Seven `container_id` channels are wired below; the manifest's
        // eighth `containers` entry is BUILD_GS_PANEL_MANIFEST's own
        // container, which that process reads from its own task.container
        // inside its own script rather than from a channel routed through
        // here (see modules/local/build_gs_panel_manifest.nf).
        BUILD_GS_PANEL_MANIFEST(
            RECONCILE_GS_PANEL_ACCOUNTING.out.accounting,
            gs_pass_for_panel_ch,
            BUILD_GS_PANEL.out.matrix.map { _meta, matrix -> matrix },
            BUILD_GS_PANEL.out.sample_metadata,
            BUILD_GS_PANEL.out.variant_metadata,
            BUILD_GS_PANEL.out.genotype_accounting,
            BUILD_GS_PANEL.out.genotype_accounting_summary,
            GATK_GATHERVCFS.out.vcf,
            reference_ch,
            reference_fai_ch,
            workflow.manifest.version,
            workflow.commitId ?: '',
            GS_NORMALIZE_VARIANTS.out.container_id,
            CLASSIFY_NORMALIZED_VARIANTS.out.container_id,
            GS_INDEX_CLASSIFIED_VARIANTS.out.container_id,
            GATK_VARIANTFILTRATION_GS.out.container_id,
            GATK_SELECTPASSVARIANTS_GS.out.container_id,
            BUILD_GS_PANEL.out.container_id,
            RECONCILE_GS_PANEL_ACCOUNTING.out.container_id,
        )

        gs_normalized_vcf_ch = GS_NORMALIZE_VARIANTS.out.vcf
        gs_classified_vcf_ch = GS_INDEX_CLASSIFIED_VARIANTS.out.vcf
        gs_filtered_vcf_ch = GATK_VARIANTFILTRATION_GS.out.vcf
        gs_pass_vcf_ch = GATK_SELECTPASSVARIANTS_GS.out.vcf
        gs_panel_matrix_ch = BUILD_GS_PANEL.out.matrix
        gs_panel_sample_metadata_ch = BUILD_GS_PANEL.out.sample_metadata
        gs_panel_variant_metadata_ch = BUILD_GS_PANEL.out.variant_metadata
        gs_panel_genotype_accounting_ch = BUILD_GS_PANEL.out.genotype_accounting
        gs_panel_record_accounting_ch = RECONCILE_GS_PANEL_ACCOUNTING.out.accounting
        gs_panel_manifest_ch = BUILD_GS_PANEL_MANIFEST.out.manifest

        // Issue #42: the GS lineage's own eight processes. These keys
        // exist in the run manifest only when the GS panel ran; with
        // enable_gs_panel=false none of these processes is invoked and
        // none of these keys appears. Note that the two `_GS` aliases are
        // recorded separately from the primary lineage's invocations of
        // the same modules, which can be overridden independently of
        // them.
        optional_container_provenance_ch = optional_container_provenance_ch
            .mix(containerProvenance(GS_NORMALIZE_VARIANTS.out.container_id, 'gs_normalize_variants'))
            .mix(containerProvenance(CLASSIFY_NORMALIZED_VARIANTS.out.container_id, 'classify_normalized_variants'))
            .mix(containerProvenance(GS_INDEX_CLASSIFIED_VARIANTS.out.container_id, 'gs_index_classified_variants'))
            .mix(containerProvenance(GATK_VARIANTFILTRATION_GS.out.container_id, 'gatk_variantfiltration_gs'))
            .mix(containerProvenance(GATK_SELECTPASSVARIANTS_GS.out.container_id, 'gatk_selectpassvariants_gs'))
            .mix(containerProvenance(BUILD_GS_PANEL.out.container_id, 'build_gs_panel'))
            .mix(containerProvenance(RECONCILE_GS_PANEL_ACCOUNTING.out.container_id, 'reconcile_gs_panel_accounting'))
            .mix(containerProvenance(BUILD_GS_PANEL_MANIFEST.out.container_id, 'build_gs_panel_manifest'))
    } else {
        // enable_gs_panel=false: skip the entire GS lineage (normalization
        // through the reproducibility manifest) without starting a single
        // one of its processes, so a non-diploid sample_ploidy can be
        // exercised end to end without ever reaching the GS panel's own
        // diploid-only fail-fast in bin/build_gs_panel.py. No
        // variants/gs_*/gs_panel/ output directories are created, since
        // Nextflow only creates a publishDir subdirectory when a process
        // actually publishes into it.
        gs_normalized_vcf_ch = channel.empty()
        gs_classified_vcf_ch = channel.empty()
        gs_filtered_vcf_ch = channel.empty()
        gs_pass_vcf_ch = channel.empty()
        gs_panel_matrix_ch = channel.empty()
        gs_panel_sample_metadata_ch = channel.empty()
        gs_panel_variant_metadata_ch = channel.empty()
        gs_panel_genotype_accounting_ch = channel.empty()
        gs_panel_record_accounting_ch = channel.empty()
        gs_panel_manifest_ch = channel.empty()
    }

    // ---------------------------------------------------------------
    // Issue #42: run-level provenance.
    //
    // Everything below turns what this run actually did into one
    // published artifact, as an ordinary required process rather than an
    // onComplete/afterScript best-effort hook. A workflow.onComplete
    // handler cannot fail the run it is reporting on, which is exactly
    // the failure mode this Issue exists to remove: "the analysis
    // succeeded but there is no provenance record" must not be a
    // successful outcome.
    // ---------------------------------------------------------------

    // Input FASTQ provenance is computed per read group, in parallel,
    // rather than by handing every raw FASTQ in the cohort to the
    // manifest process: at this pipeline's 327-sample target that design
    // would re-stage the entire input dataset into a single task purely
    // to write a JSON file.
    HASH_INPUT_FASTQS(input_provenance_rows_ch)

    // The reference bundle as it was actually used for mapping -- FASTA,
    // FAI, dictionary and the five BWA-MEM2 index files, whether this run
    // generated the index or was given a prebuilt one.
    HASH_REFERENCE_BUNDLE(
        reference_ch,
        reference_fai_ch,
        reference_dict_ch,
        bwa_indexes_ch,
    )

    // The run's own scientific deliverables, checksummed in per-group
    // tasks (one per sample gVCF, one per cohort-level VCF) for the same
    // scaling reason. This is the artifact set the historical schema v1
    // manifests recorded, unchanged.
    run_artifact_groups_ch = GATK_HAPLOTYPECALLER.out.gvcf
        .map { meta, gvcf, _gvcf_index -> tuple([id: "${meta.id}.gvcf"], [gvcf]) }
        .mix(
            GATK_GATHERVCFS.out.vcf.map { meta, vcf, _vcf_index ->
                tuple([id: "${meta.id}.raw"], [vcf])
            }
        )
        .mix(
            GATK_SELECTPASSVARIANTS.out.vcf.map { meta, vcf, _vcf_index ->
                tuple([id: "${meta.id}.${meta.variant_type}.pass"], [vcf])
            }
        )

    HASH_RUN_ARTIFACTS(run_artifact_groups_ch)

    // One `process_key<TAB>effective_container` row per executed task.
    // Every containerized process contributes its own task.container --
    // Nextflow's post-override effective value -- so the manifest records
    // what ran, not what the `container` directives defaulted to. The
    // processes listed here run in every configuration; conditional ones
    // were accumulated into optional_container_provenance_ch by the
    // branches that invoked them. BUILD_RUN_MANIFEST adds its own
    // identity from inside its own task.
    //
    // Collecting this is also what makes the manifest a genuine terminal
    // barrier: BUILD_RUN_MANIFEST consumes the aggregate, so it cannot
    // start until every one of these processes -- MULTIQC and the rest of
    // the QC side branch included -- has finished.
    runtime_container_provenance_ch = channel.empty()
        .mix(containerProvenance(FASTQC_RAW.out.container_id, 'fastqc_raw'))
        .mix(containerProvenance(FASTQC_TRIMMED.out.container_id, 'fastqc_trimmed'))
        .mix(containerProvenance(FASTP.out.container_id, 'fastp'))
        .mix(containerProvenance(VALIDATE_REFERENCE_CONTIGS.out.container_id, 'validate_reference_contigs'))
        .mix(containerProvenance(BWA_MEM2_MEM_SORT.out.container_id, 'bwa_mem2_mem_sort'))
        .mix(containerProvenance(SAMTOOLS_MERGE.out.container_id, 'samtools_merge'))
        .mix(containerProvenance(GATK_MARKDUPLICATES.out.container_id, 'gatk_markduplicates'))
        .mix(containerProvenance(SAMTOOLS_INDEX.out.container_id, 'samtools_index'))
        .mix(containerProvenance(SAMTOOLS_QC.out.container_id, 'samtools_qc'))
        .mix(containerProvenance(MULTIQC.out.container_id, 'multiqc'))
        .mix(containerProvenance(GATK_HAPLOTYPECALLER.out.container_id, 'gatk_haplotypecaller'))
        .mix(containerProvenance(GATK_GENOMICSDBIMPORT.out.container_id, 'gatk_genomicsdbimport'))
        .mix(containerProvenance(GATK_GENOTYPEGVCFS.out.container_id, 'gatk_genotypegvcfs'))
        .mix(containerProvenance(GATK_GATHERVCFS.out.container_id, 'gatk_gathervcfs'))
        .mix(containerProvenance(GATK_SELECTVARIANTS.out.container_id, 'gatk_selectvariants'))
        .mix(containerProvenance(GATK_VARIANTFILTRATION.out.container_id, 'gatk_variantfiltration'))
        .mix(containerProvenance(GATK_SELECTPASSVARIANTS.out.container_id, 'gatk_selectpassvariants'))
        .mix(containerProvenance(BCFTOOLS_STATS.out.container_id, 'bcftools_stats'))
        .mix(containerProvenance(SUMMARIZE_VARIANT_QC.out.container_id, 'summarize_variant_qc'))
        .mix(containerProvenance(SUMMARIZE_FILTER_QC.out.container_id, 'summarize_filter_qc'))
        .mix(containerProvenance(RECONCILE_VARIANT_TYPE_COUNTS.out.container_id, 'reconcile_variant_type_counts'))
        .mix(containerProvenance(HASH_INPUT_FASTQS.out.container_id, 'hash_input_fastqs'))
        .mix(containerProvenance(HASH_REFERENCE_BUNDLE.out.container_id, 'hash_reference_bundle'))
        .mix(containerProvenance(HASH_RUN_ARTIFACTS.out.container_id, 'hash_run_artifacts'))
        .mix(optional_container_provenance_ch)
        .collectFile(
            name: 'runtime_container_provenance.tsv',
            newLine: true,
            sort: true,
        )

    // The cohort and per-sample accounting the pipeline already
    // published for the raw/all cohort VCF. Both come from the same
    // SUMMARIZE_VARIANT_QC invocation, so the manifest builder can
    // cross-check them against each other instead of recomputing either.
    run_manifest_accounting_ch = SUMMARIZE_VARIANT_QC.out.qc
        .filter { meta, _variant_qc_tsv, _sample_qc_tsv, _summary_txt ->
            meta['qc_stage'] == 'raw' && meta['variant_type'] == 'all'
        }
        .map { meta, variant_qc_tsv, sample_qc_tsv, _summary_txt ->
            tuple([id: meta.id], variant_qc_tsv, sample_qc_tsv)
        }

    // Exactly one value reaches BUILD_RUN_MANIFEST in both
    // configurations. With the GS panel disabled no GS manifest exists,
    // so a placeholder is staged and `--no-gs-panel` is passed instead:
    // wiring channel.empty() into a required input would simply deadlock
    // the process and hang the run.
    if (gs_panel_enabled) {
        gs_panel_manifest_for_run_ch = BUILD_GS_PANEL_MANIFEST.out.manifest
            .map { _meta, manifest_json -> manifest_json }
    } else {
        gs_panel_manifest_for_run_ch = channel.value(
            file("${projectDir}/assets/NO_GS_PANEL_MANIFEST", checkIfExists: true)
        )
    }

    BUILD_RUN_MANIFEST(
        run_manifest_accounting_ch,
        RECONCILE_VARIANT_TYPE_COUNTS.out.accounting.map {
            meta,
            accounting_tsv,
            _summary_txt ->
            tuple(meta, accounting_tsv)
        },
        runtime_container_provenance_ch,
        HASH_INPUT_FASTQS.out.provenance.map { _meta, tsv -> tsv }.collect(),
        HASH_REFERENCE_BUNDLE.out.provenance.map { _meta, tsv -> tsv },
        HASH_RUN_ARTIFACTS.out.checksums.map { _meta, tsv -> tsv }.collect(),
        workflow.manifest.version,
        run_git_commit,
        nextflow.version.toString(),
        gs_panel_enabled,
        gs_panel_manifest_for_run_ch,
    )

    emit:
    raw_fastqc_html = FASTQC_RAW.out.html
    raw_fastqc_zip = FASTQC_RAW.out.zip
    trimmed_reads = FASTP.out.reads
    fastp_reports = FASTP.out.reports
    trimmed_fastqc_html = FASTQC_TRIMMED.out.html
    trimmed_fastqc_zip = FASTQC_TRIMMED.out.zip
    mapping_logs = BWA_MEM2_MEM_SORT.out.log
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
    gs_normalized_vcf = gs_normalized_vcf_ch
    gs_classified_vcf = gs_classified_vcf_ch
    gs_filtered_vcf = gs_filtered_vcf_ch
    gs_pass_vcf = gs_pass_vcf_ch
    gs_panel_matrix = gs_panel_matrix_ch
    gs_panel_sample_metadata = gs_panel_sample_metadata_ch
    gs_panel_variant_metadata = gs_panel_variant_metadata_ch
    gs_panel_genotype_accounting = gs_panel_genotype_accounting_ch
    gs_panel_record_accounting = gs_panel_record_accounting_ch
    gs_panel_manifest = gs_panel_manifest_ch
    flagstat = SAMTOOLS_QC.out.flagstat
    stats = SAMTOOLS_QC.out.stats
    idxstats = SAMTOOLS_QC.out.idxstats
    multiqc_report = MULTIQC.out.report
    multiqc_data = MULTIQC.out.data
    multiqc_config = MULTIQC.out.config
    multiqc_version = MULTIQC.out.version
    reference_fai = reference_fai_ch
    reference_dict = reference_dict_ch
    bwa_indexes = bwa_indexes_ch
    input_provenance = HASH_INPUT_FASTQS.out.provenance
    reference_provenance = HASH_REFERENCE_BUNDLE.out.provenance
    run_artifact_checksums = HASH_RUN_ARTIFACTS.out.checksums
    run_manifest = BUILD_RUN_MANIFEST.out.manifest
}
