// Issue #30 introduced this dedicated resource label after the original
// all-record implementation exceeded process_low on a real 5-sample cohort.
// Issue #35 removed that O(N-records) representation: classification and
// output now stream by (CHROM, POS), retaining only the current locus while
// preserving all-occurrence duplicate-key exclusion. Formal 10/20-sample
// replays at 70ae4d6 measured the Python process at about 21 MiB at both
// scales with zero swap growth and production-equivalent output. See
// nextflow.config for the separate page-cache-aware resource contract.
process CLASSIFY_NORMALIZED_VARIANTS {
    tag "${meta.id}"
    label 'process_variant_classification'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required. This image also has no
    // bgzip/tabix, which is why this process writes a plain-text VCF;
    // modules/local/gs_index_classified_variants.nf compresses and
    // indexes it in a separate, bcftools-based process.
    container params.containers.python

    input:
    tuple val(meta), path(normalized_vcf), path(normalized_vcf_index)

    output:
    tuple(val(meta), path("${meta.id}.classified.vcf"), emit: vcf)
    path("${meta.id}.classification_accounting.tsv"), emit: accounting
    path("${meta.id}.classification_accounting.summary.txt"), emit: summary
    // Issue #52: see gs_normalize_variants.nf for why this records
    // task.container (the resolved, post-override effective container)
    // rather than trusting the `container` directive above.
    val(task.container), emit: container_id

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
