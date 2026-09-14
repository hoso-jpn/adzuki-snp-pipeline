// Issue #64: render the optional genotype quality policy from params.
//
// Values are checked here as well as by nextflow_schema.json, because this
// module is also run directly (tests/modules/build_gs_panel.nf.test) where no
// schema validation happens, and each value lands on a shell command line.
// With the mask off nothing is added, so the command -- and every output --
// is exactly what it was before this Issue.
def genotypeQualityArguments(prefix) {
    if (!params.gs_genotype_quality_mask.toString().toBoolean()) {
        return ''
    }
    def policies = ['reject', 'unevaluated', 'mask']
    def arguments = ['--genotype-quality-mask']
    [dp: params.gs_genotype_min_dp, gq: params.gs_genotype_min_gq].each { field, value ->
        if (value != null) {
            def threshold = value.toString()
            if (!(threshold ==~ /[0-9]+/)) {
                error("gs_genotype_min_${field} must be a non-negative integer, got '${threshold}'")
            }
            arguments << "--genotype-min-${field} ${threshold}"
        }
    }
    [
        'missing-format-field': params.gs_genotype_missing_format_field,
        'missing-value': params.gs_genotype_missing_value,
        'malformed-value': params.gs_genotype_malformed_value,
    ].each { option, value ->
        if (!(value.toString() in policies)) {
            error("gs_genotype_${option.replace('-', '_')} must be one of ${policies}, got '${value}'")
        }
        arguments << "--genotype-${option}-policy ${value}"
    }
    arguments << "--quality-masked-vcf-output ${prefix}.quality_masked.vcf.gz"
    arguments << "--genotype-quality-policy-output ${prefix}.genotype_quality_policy.json"
    // Leading space included here, so with the mask off the rendered
    // command has no trailing whitespace and is exactly the historical one.
    return ' ' + arguments.join(' ')
}

// Issue #30: dedicated resource label (see nextflow.config for the
// real-data benchmark this was sized from). On Issue #26's real
// 5-sample cohort this process's true peak RSS (5.34 GiB, measured at
// a generous memory ceiling) exceeded process_low's previous 4 GiB
// first-attempt allocation -- the original run's own report of
// "peak_rss: 4 GB" was the cgroup ceiling itself, not genuine headroom.
process BUILD_GS_PANEL {
    tag "${meta.id}"
    label 'process_gs_panel'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

    input:
    tuple val(meta), path(gs_pass_vcf), path(gs_pass_vcf_index)

    output:
    tuple(val(meta), path("${meta.id}.gs_panel.genotype_matrix.tsv.gz"), emit: matrix)
    path("${meta.id}.gs_panel.sample_metadata.tsv"), emit: sample_metadata
    path("${meta.id}.gs_panel.variant_metadata.tsv"), emit: variant_metadata
    path("${meta.id}.gs_panel.genotype_encoding_accounting.tsv"), emit: genotype_accounting
    path("${meta.id}.gs_panel.genotype_encoding_accounting.summary.txt"), emit: genotype_accounting_summary
    // Issue #64: produced only when params.gs_genotype_quality_mask is true.
    // The masked VCF is BGZF written by the same pass as the matrix; it is
    // indexed and published by GS_INDEX_QUALITY_MASKED_VCF.
    tuple(val(meta), path("${meta.id}.gs_panel.quality_masked.vcf.gz"), emit: quality_masked_vcf, optional: true)
    path("${meta.id}.gs_panel.genotype_quality_policy.json"), emit: quality_policy, optional: true
    // Issue #52: see modules/local/gs_normalize_variants.nf for why this
    // records task.container (the resolved, post-override effective
    // container) rather than trusting the `container` directive above.
    val(task.container), emit: container_id

    script:
    prefix = "${meta.id}.gs_panel"
    quality_args = genotypeQualityArguments(prefix)

    """
    build_gs_panel.py \
        --gs-pass-vcf ${gs_pass_vcf} \
        --cohort-id '${meta.id}' \
        --sample-ploidy ${params.sample_ploidy} \
        --matrix-output ${prefix}.genotype_matrix.tsv.gz \
        --sample-metadata-output ${prefix}.sample_metadata.tsv \
        --variant-metadata-output ${prefix}.variant_metadata.tsv \
        --genotype-accounting-output ${prefix}.genotype_encoding_accounting.tsv \
        --genotype-accounting-summary-output ${prefix}.genotype_encoding_accounting.summary.txt${quality_args}
    """
}
