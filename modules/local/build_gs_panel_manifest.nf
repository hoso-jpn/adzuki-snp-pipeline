process BUILD_GS_PANEL_MANIFEST {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(record_accounting), path(record_accounting_summary)
    tuple val(gs_pass_meta), path(gs_pass_vcf), path(gs_pass_vcf_index)
    path(matrix)
    path(sample_metadata)
    path(variant_metadata)
    path(genotype_accounting)
    path(genotype_accounting_summary)
    tuple val(raw_all_meta), path(raw_all_vcf), path(raw_all_vcf_index)
    tuple val(reference_meta), path(reference_fasta)
    tuple val(reference_fai_meta), path(reference_fai)
    val(pipeline_version)
    val(git_commit)
    // Issue #52: one effective container identity per GS-lineage process,
    // each read from that process's own `container_id` output (Nextflow's
    // task.container, resolved after any withName/alias/fully-qualified-
    // selector/profile override) rather than a shared per-tool literal --
    // see workflows/adzuki_snp_pipeline.nf for how these are wired and
    // docs/gs_panel_data_contract.md for the schema v2 rationale.
    val(gs_normalize_variants_container)
    val(classify_normalized_variants_container)
    val(gs_index_classified_variants_container)
    val(gatk_variantfiltration_gs_container)
    val(gatk_selectpassvariants_gs_container)
    val(build_gs_panel_container)
    val(reconcile_gs_panel_accounting_container)

    output:
    tuple(val(meta), path("${meta.id}.gs_panel.manifest.json"), emit: manifest)

    script:
    """
    build_gs_panel_manifest.py \
        --cohort-id '${meta.id}' \
        --pipeline-version '${pipeline_version}' \
        --git-commit '${git_commit}' \
        --container-gs-normalize-variants '${gs_normalize_variants_container}' \
        --container-classify-normalized-variants '${classify_normalized_variants_container}' \
        --container-gs-index-classified-variants '${gs_index_classified_variants_container}' \
        --container-gatk-variantfiltration-gs '${gatk_variantfiltration_gs_container}' \
        --container-gatk-selectpassvariants-gs '${gatk_selectpassvariants_gs_container}' \
        --container-build-gs-panel '${build_gs_panel_container}' \
        --container-reconcile-gs-panel-accounting '${reconcile_gs_panel_accounting_container}' \
        --sample-ploidy ${params.sample_ploidy} \
        --snp-filter-qd-min ${params.snp_filter_qd_min} \
        --snp-filter-qual-min ${params.snp_filter_qual_min} \
        --snp-filter-sor-max ${params.snp_filter_sor_max} \
        --snp-filter-fs-max ${params.snp_filter_fs_max} \
        --snp-filter-mq-min ${params.snp_filter_mq_min} \
        --snp-filter-mq-rank-sum-min ${params.snp_filter_mq_rank_sum_min} \
        --snp-filter-read-pos-rank-sum-min ${params.snp_filter_read_pos_rank_sum_min} \
        --record-accounting ${record_accounting} \
        --checksum-file ${gs_pass_vcf} \
        --checksum-file ${matrix} \
        --checksum-file ${sample_metadata} \
        --checksum-file ${variant_metadata} \
        --checksum-file ${genotype_accounting} \
        --checksum-file ${genotype_accounting_summary} \
        --checksum-file ${record_accounting} \
        --checksum-file ${record_accounting_summary} \
        --checksum-file ${raw_all_vcf} \
        --checksum-file ${reference_fasta} \
        --checksum-file ${reference_fai} \
        --output ${meta.id}.gs_panel.manifest.json
    """
}
