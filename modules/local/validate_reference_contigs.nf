// Issue #11: neither the pipeline-generated reference-bundle path
// (SAMTOOLS_FAIDX + GATK_CREATE_SEQUENCE_DICTIONARY) nor the prebuilt path
// (--reference_fai/--reference_dict) previously confirmed that the .fai and
// .dict actually describe the same reference, in the same contig order,
// with the same declared lengths, before GATK_HAPLOTYPECALLER/
// GATK_GENOMICSDBIMPORT/GATK_GENOTYPEGVCFS started consuming both. A
// mismatched pair (most concretely: a correct prebuilt .fai/.dict whose
// contigs are the same *set* but listed in a different order -- something
// no downstream GATK tool is guaranteed to reject up front) would otherwise
// only surface as whatever error message, at whatever point in a run, the
// first GATK process that happens to notice chooses to produce.
//
// This process runs once, right after reference_fai_ch/reference_dict_ch
// are established in workflows/adzuki_snp_pipeline.nf -- for *both* the
// generated and the prebuilt path, since both converge to the same
// (meta, file) channel shape at that point -- and every downstream
// reference-dependent process depends on its output rather than on the
// pre-validation channels directly, so a mismatch fails the whole run
// before any GATK process starts, with an actionable message identifying
// the first point of disagreement (see bin/validate_reference_contigs.py).
process VALIDATE_REFERENCE_CONTIGS {
    tag "${fai_meta.id}"
    label 'process_low'

    container 'python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7'

    input:
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dict)

    output:
    tuple val(fai_meta), path(fai), emit: fai
    tuple val(dict_meta), path(dict), emit: dict

    script:
    """
    validate_reference_contigs.py \
        --fai ${fai} \
        --dict ${dict}
    """
}
