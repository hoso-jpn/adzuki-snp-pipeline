// Issue #30: dedicated resource label (see nextflow.config for the
// real-data benchmark this was sized from). This process reads the
// entire cohort's normalized VCF into memory at once; that footprint
// exceeded process_low's previous 4 GiB first-attempt allocation on
// Issue #26's real 5-sample cohort.
process CLASSIFY_NORMALIZED_VARIANTS {
    tag "${meta.id}"
    label 'process_variant_classification'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required. This image also has no
    // bgzip/tabix, which is why this process writes a plain-text VCF;
    // modules/local/gs_index_classified_variants.nf compresses and
    // indexes it in a separate, bcftools-based process.
    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(meta), path(normalized_vcf), path(normalized_vcf_index)

    output:
    tuple(val(meta), path("${meta.id}.classified.vcf"), emit: vcf)
    path("${meta.id}.classification_accounting.tsv"), emit: accounting
    path("${meta.id}.classification_accounting.summary.txt"), emit: summary

    script:
    """
    classify_normalized_variants.py \
        --normalized-vcf ${normalized_vcf} \
        --cohort-id '${meta.id}' \
        --output ${meta.id}.classified.vcf \
        --accounting-output ${meta.id}.classification_accounting.tsv \
        --summary-output ${meta.id}.classification_accounting.summary.txt
    """
}
