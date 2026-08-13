process BCFTOOLS_STATS {
    tag "${meta.id}:${meta.qc_stage}:${meta.variant_type}"
    label 'process_low'

    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.bcftools.stats.tsv"),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.variant_qc.tsv"),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.sample_qc.tsv"),
        path("${meta.id}.${meta.qc_stage}.${meta.variant_type}.summary.txt"),
        emit: qc
    )

    script:
    prefix = "${meta.id}.${meta.qc_stage}.${meta.variant_type}"

    """
    bcftools stats \
        --threads ${task.cpus} \
        --samples - \
        ${vcf} \
        > ${prefix}.bcftools.stats.tsv

    awk \
        -F '\\t' \
        -v OFS='\\t' \
        -v cohort_id='${meta.id}' \
        -v stage='${meta.qc_stage}' \
        -v variant_type='${meta.variant_type}' \
        -v variant_out='${prefix}.variant_qc.tsv' \
        -v sample_out='${prefix}.sample_qc.tsv' \
        -v summary_out='${prefix}.summary.txt' \
        '
        \$1 == "SN" {
            if (\$3 == "number of samples:") {
                number_of_samples = \$4
            } else if (\$3 == "number of records:") {
                number_of_records = \$4
            } else if (\$3 == "number of SNPs:") {
                number_of_snps = \$4
            } else if (\$3 == "number of MNPs:") {
                number_of_mnps = \$4
            } else if (\$3 == "number of indels:") {
                number_of_indels = \$4
            } else if (\$3 == "number of others:") {
                number_of_others = \$4
            } else if (\$3 == "number of multiallelic sites:") {
                number_of_multiallelic_sites = \$4
            }
        }

        \$1 == "TSTV" {
            transitions = \$3
            transversions = \$4
            transition_transversion_ratio = \$5
        }

        \$1 == "PSC" {
            sample_count++
            sample_name[sample_count] = \$3
            reference_homozygous[sample_count] = \$4
            non_reference_homozygous[sample_count] = \$5
            heterozygous[sample_count] = \$6
            average_depth[sample_count] = \$10
            singleton_count[sample_count] = \$11
            missing_count[sample_count] = \$14
        }

        END {
            number_of_samples += 0
            number_of_records += 0
            number_of_snps += 0
            number_of_mnps += 0
            number_of_indels += 0
            number_of_others += 0
            number_of_multiallelic_sites += 0
            transitions += 0
            transversions += 0

            if (transition_transversion_ratio == "") {
                transition_transversion_ratio = "NA"
            }

            sample_names = ""
            cohort_missing = 0

            print \
                "cohort_id", \
                "stage", \
                "variant_type", \
                "metric", \
                "value" \
                > variant_out

            print \
                cohort_id, stage, variant_type, \
                "number_of_samples", number_of_samples \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_records", number_of_records \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_snps", number_of_snps \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_mnps", number_of_mnps \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_indels", number_of_indels \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_others", number_of_others \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "number_of_multiallelic_sites", \
                number_of_multiallelic_sites \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "transitions", transitions \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "transversions", transversions \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "transition_transversion_ratio", \
                transition_transversion_ratio \
                >> variant_out

            print \
                "cohort_id", \
                "stage", \
                "variant_type", \
                "sample", \
                "reference_homozygous", \
                "non_reference_homozygous", \
                "heterozygous", \
                "missing", \
                "missingness_rate", \
                "average_depth", \
                "singletons" \
                > sample_out

            for (sample_index = 1; sample_index <= sample_count; sample_index++) {
                cohort_missing += missing_count[sample_index]

                if (sample_names == "") {
                    sample_names = sample_name[sample_index]
                } else {
                    sample_names = \
                        sample_names "," sample_name[sample_index]
                }

                if (number_of_records > 0) {
                    sample_missingness = sprintf("%.6f", missing_count[sample_index] / number_of_records)
                } else {
                    sample_missingness = "NA"
                }

                print \
                    cohort_id, \
                    stage, \
                    variant_type, \
                    sample_name[sample_index], \
                    reference_homozygous[sample_index], \
                    non_reference_homozygous[sample_index], \
                    heterozygous[sample_index], \
                    missing_count[sample_index], \
                    sample_missingness, \
                    average_depth[sample_index], \
                    singleton_count[sample_index] \
                    >> sample_out
            }

            cohort_total = \
                number_of_records * number_of_samples

            if (cohort_total > 0) {
                cohort_missingness = sprintf("%.6f", cohort_missing / cohort_total)
            } else {
                cohort_missingness = "NA"
            }

            print \
                cohort_id, stage, variant_type, \
                "cohort_missing_genotypes", cohort_missing \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "cohort_total_genotypes", cohort_total \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "cohort_missingness_rate", cohort_missingness \
                >> variant_out
            print \
                cohort_id, stage, variant_type, \
                "sample_names", sample_names \
                >> variant_out

            print "Variant QC summary" > summary_out
            print "Cohort ID: " cohort_id >> summary_out
            print "Stage: " stage >> summary_out
            print "Variant type: " variant_type >> summary_out
            print \
                "Samples (" number_of_samples "): " sample_names \
                >> summary_out
            print \
                "Records: " number_of_records \
                >> summary_out
            print \
                "SNPs: " number_of_snps \
                >> summary_out
            print \
                "MNPs: " number_of_mnps \
                >> summary_out
            print \
                "Indels: " number_of_indels \
                >> summary_out
            print \
                "Other variants: " number_of_others \
                >> summary_out
            print \
                "Multiallelic sites: " \
                number_of_multiallelic_sites \
                >> summary_out
            print \
                "Transitions: " transitions \
                >> summary_out
            print \
                "Transversions: " transversions \
                >> summary_out
            print \
                "Ti/Tv ratio: " transition_transversion_ratio \
                >> summary_out
            print \
                "Missing genotypes: " \
                cohort_missing "/" cohort_total \
                " (" cohort_missingness ")" \
                >> summary_out
        }
        ' \
        ${prefix}.bcftools.stats.tsv
    """
}
