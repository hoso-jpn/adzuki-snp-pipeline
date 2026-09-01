process BUILD_RUN_MANIFEST {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required. Sourced from
    // conf/containers.config (Issue #52) rather than re-declaring the
    // digest here.
    container params.containers.python

    input:
    tuple val(meta), path(variant_qc_tsv), path(sample_qc_tsv)
    tuple val(accounting_meta), path(variant_type_accounting_tsv)
    // Issue #42: one `process_key<TAB>container` row per executed task,
    // aggregated by the workflow from every containerized process's own
    // `container_id` output. Because this single input depends on all of
    // them, BUILD_RUN_MANIFEST cannot start until every process that
    // contributes to it has finished -- including the MULTIQC side
    // branch, which would otherwise be free to still be running while
    // the run's "final" provenance artifact was already written.
    path(runtime_provenance)
    path(input_provenance)
    path(reference_provenance)
    path(artifact_checksums)
    val(pipeline_version)
    val(git_commit)
    val(nextflow_version)
    val(gs_panel_enabled)
    // The GS panel manifest when the GS lineage ran, and
    // assets/NO_GS_PANEL_MANIFEST (never read) when it did not -- see
    // that file for why a placeholder rather than an optional input.
    path(gs_panel_manifest)

    output:
    tuple(val(meta), path("${meta.id}.run_manifest.json"), emit: manifest)

    script:
    // Issue #42: this process's own effective container, read from its
    // own task.container. Nextflow resolves a task's container before
    // rendering that task's script, so a process can read its own
    // post-override value directly (verified on Nextflow 26.04.6 with
    // both a default and a withName override, the same way
    // BUILD_GS_PANEL_MANIFEST does it). Routing it back through a
    // channel would be a self-referential wiring step that makes the
    // recorded value no more accurate.
    def runtime_provenance_args = ([runtime_provenance].flatten())
        .collect { provenance -> "--runtime-provenance ${provenance}" }
        .join(' \\\n        ')
    def input_provenance_args = ([input_provenance].flatten())
        .collect { provenance -> "--input-provenance ${provenance}" }
        .join(' \\\n        ')
    def artifact_checksum_args = ([artifact_checksums].flatten())
        .collect { checksums -> "--artifact-checksums ${checksums}" }
        .join(' \\\n        ')
    def gs_panel_arg = gs_panel_enabled
        ? "--gs-panel-manifest ${gs_panel_manifest}"
        : '--no-gs-panel'
    // See main.nf for why enable_gs_panel is round-tripped through
    // toString().toBoolean(): a CLI-provided `--enable_gs_panel false`
    // arrives as the String "false", which is truthy in Groovy.
    def enable_gs_panel_flag = params.enable_gs_panel.toString().toBoolean()
        ? '--enable-gs-panel'
        : '--no-enable-gs-panel'

    """
    printf 'build_run_manifest\\t%s\\n' '${task.container}' \\
        > build_run_manifest.container.tsv

    build_run_manifest.py \
        --mode dag-v2 \
        --cohort-id '${meta.id}' \
        --pipeline-version '${pipeline_version}' \
        --git-commit '${git_commit}' \
        --nextflow-version '${nextflow_version}' \
        --reference-id '${params.reference_id}' \
        --reference-name '${params.reference_name}' \
        --reference-species '${params.reference_species}' \
        --reference-cultivar '${params.reference_cultivar}' \
        --reference-accession '${params.reference_accession}' \
        --sample-ploidy ${params.sample_ploidy} \
        --genomicsdb-batch-size ${params.genomicsdb_batch_size} \
        --optical-duplicate-pixel-distance ${params.optical_duplicate_pixel_distance} \
        ${enable_gs_panel_flag} \
        --snp-filter-qd-min ${params.snp_filter_qd_min} \
        --snp-filter-qual-min ${params.snp_filter_qual_min} \
        --snp-filter-sor-max ${params.snp_filter_sor_max} \
        --snp-filter-fs-max ${params.snp_filter_fs_max} \
        --snp-filter-mq-min ${params.snp_filter_mq_min} \
        --snp-filter-mq-rank-sum-min ${params.snp_filter_mq_rank_sum_min} \
        --snp-filter-read-pos-rank-sum-min ${params.snp_filter_read_pos_rank_sum_min} \
        --indel-filter-qd-min ${params.indel_filter_qd_min} \
        --indel-filter-qual-min ${params.indel_filter_qual_min} \
        --indel-filter-fs-max ${params.indel_filter_fs_max} \
        --indel-filter-read-pos-rank-sum-min ${params.indel_filter_read_pos_rank_sum_min} \
        ${runtime_provenance_args} \
        --runtime-provenance build_run_manifest.container.tsv \
        ${input_provenance_args} \
        --reference-provenance ${reference_provenance} \
        ${artifact_checksum_args} \
        --variant-qc-tsv ${variant_qc_tsv} \
        --sample-qc-tsv ${sample_qc_tsv} \
        --variant-type-accounting-tsv ${variant_type_accounting_tsv} \
        ${gs_panel_arg} \
        --output ${meta.id}.run_manifest.json
    """
}
