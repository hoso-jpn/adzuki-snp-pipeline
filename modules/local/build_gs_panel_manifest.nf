process BUILD_GS_PANEL_MANIFEST {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required. Issue #52: the digest
    // itself lives only in conf/containers.config, like every other
    // GS-lineage process's -- this process writes the manifest that
    // records those identities, so a second copy of the same default
    // here is exactly the drift this Issue exists to remove.
    container params.containers.python

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
    // This process's own effective container is deliberately NOT an
    // input: it is read below as `task.container` in this process's own
    // script. Nextflow resolves a task's container before rendering its
    // script, so a process can read its own post-override effective
    // value directly -- confirmed on Nextflow 26.04.6 with both a
    // default and a `withName` override -- and routing it through a
    // channel instead would only add a self-referential wiring step
    // without making the recorded value any more accurate.
    val(gs_normalize_variants_container)
    val(classify_normalized_variants_container)
    val(gs_index_classified_variants_container)
    val(gatk_variantfiltration_gs_container)
    val(gatk_selectpassvariants_gs_container)
    val(build_gs_panel_container)
    val(reconcile_gs_panel_accounting_container)

    output:
    tuple(val(meta), path("${meta.id}.gs_panel.manifest.json"), emit: manifest)
    // Issue #42: the same effective container this process already
    // records for itself inside the GS manifest (see `--container-build-
    // gs-panel-manifest` below), emitted as a channel as well so the
    // run-level manifest can record it the same way it records every
    // other process -- without BUILD_RUN_MANIFEST having to re-derive it
    // by reading the GS manifest back out of its own JSON.
    val(task.container), emit: container_id

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
        --container-build-gs-panel-manifest '${task.container}' \
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
