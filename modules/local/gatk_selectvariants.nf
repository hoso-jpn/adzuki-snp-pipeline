// Issue #17: `meta.variant_type` is this pipeline's single canonical
// representation of "which variant type is this lineage", and it is
// deliberately lowercase (`snp`/`indel`) because that same value is
// interpolated into published artifact names (cohort.snp.vcf.gz,
// cohort.indel.filtered.vcf.gz), into process tags, and into the
// stage/type keys the QC accounting files are grouped by. GATK's own
// `--select-type-to-include` wants the uppercase VariantContext type
// name instead, so this table -- not a second `val(variant_type)`
// process input carrying `SNP`/`INDEL` alongside the meta tag -- is
// where that one translation happens. Passing both made the pair
// `[variant_type: 'snp'] + 'INDEL'` a representable (and silently
// wrong: snp-named output holding indel records) state; deriving the
// GATK argument from the canonical tag makes it unrepresentable.
def gatkSelectTypeFor(meta) {
    def select_types = [
        snp: 'SNP',
        indel: 'INDEL',
    ]
    def expected = select_types.keySet().join(', ')

    if (!(meta instanceof Map) || !meta.containsKey('variant_type')) {
        error(
            "GATK_SELECTVARIANTS requires meta.variant_type; " +
                "expected one of: ${expected}"
        )
    }

    def variant_type = meta['variant_type']

    // Deliberately not `variant_type.toUpperCase()` with no allowlist:
    // that would accept `mnp`/`mixed`/`foo` and hand GATK an argument
    // it either rejects with a much less specific message or -- for
    // `MNP`/`MIXED`, which are real VariantContext types -- silently
    // accepts, producing a cohort.mnp.vcf.gz that no downstream stage
    // in this pipeline knows about. Case is also not normalized here:
    // accepting `SNP` would publish cohort.SNP.vcf.gz and split the QC
    // accounting keys, so a non-canonical spelling is an error rather
    // than something quietly repaired.
    def select_type = variant_type instanceof CharSequence
        ? select_types[variant_type.toString()]
        : null

    if (!select_type) {
        error(
            "unsupported variant type '${variant_type}' for " +
                "${meta['id'] ?: 'unknown'}; expected one of: ${expected}"
        )
    }

    return select_type
}

process GATK_SELECTVARIANTS {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple(
        val(meta),
        path(vcf),
        path(vcf_index)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.variant_type}.vcf.gz"),
        path("${meta.id}.${meta.variant_type}.vcf.gz.tbi"),
        emit: vcf
    )

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    // Validated at this process boundary rather than only in
    // workflows/adzuki_snp_pipeline.nf, so the module stays safe when
    // used standalone (as tests/modules/gatk_selectvariants.nf.test
    // does) and so the failure is a diagnosable Nextflow error naming
    // the rejected value, not a confusing GATK usage error or a
    // mis-typed output file.
    select_type = gatkSelectTypeFor(meta)

    """
    gatk --java-options "-Xmx${memory_gb}g" SelectVariants \
        --variant ${vcf} \
        --select-type-to-include ${select_type} \
        --output ${meta.id}.${meta.variant_type}.vcf.gz \
        --create-output-variant-index true
    """
}
