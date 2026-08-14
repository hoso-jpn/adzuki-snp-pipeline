#!/usr/bin/env nextflow

include {
    validateParameters
    paramsSummaryLog
    samplesheetToList
} from 'plugin/nf-schema'

include {
    ADZUKI_SNP_PIPELINE
} from './workflows/adzuki_snp_pipeline'

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // `.toString().toBoolean()` (not a bare `if (params.enable_gs_panel)`)
    // because Nextflow resolves a CLI-provided `--enable_gs_panel false`
    // to the *String* "false", not the Groovy boolean `false` -- and a
    // non-empty String, including the literal text "false", is truthy.
    // Round-tripping through `.toString()` first makes this correct
    // whether the value came from the CLI (String) or from
    // nextflow.config's own `true` default (already a real Boolean).
    def gsPanelEnabled = params.enable_gs_panel.toString().toBoolean()

    // The same class of bug as the boolean coercion above, in the
    // opposite direction: a CLI-provided `--sample_ploidy 2` resolves
    // to the *String* "2", not the Groovy Integer 2. Groovy's `!=`
    // does not coerce a String and an Integer to compare by numeric
    // value, so `"2" != 2` evaluates to `true` -- meaning a user who
    // explicitly, correctly passes `--sample_ploidy 2` would otherwise
    // hit the fail-fast below by mistake. `.toString().toInteger()`
    // first is correct whether the value came from the CLI (String) or
    // from nextflow.config's own `2` default (already a real Integer).
    def samplePloidy = params.sample_ploidy.toString().toInteger()

    if (gsPanelEnabled && samplePloidy != 2) {
        error(
            "params.sample_ploidy is ${params.sample_ploidy}, but params.enable_gs_panel " +
            'is true (the default). The GS panel schema (v1) is diploid-only and would ' +
            'fail only after variant calling has already run, wasting that work. Set ' +
            '--enable_gs_panel false to run non-diploid variant calling without the GS ' +
            'panel, or leave sample_ploidy at its default (2) to keep the GS panel enabled.'
        )
    }

    sample_rows = samplesheetToList(
        params.input,
        'assets/schema_input.json'
    )

    invalid_read_pairs = sample_rows
        .findAll { row -> row[1] == row[2] }
        .collect { row -> row[0].read_group_id }

    if (!invalid_read_pairs.isEmpty()) {
        error(
            'fastq_1 and fastq_2 must reference different files for ' +
            "read groups: ${invalid_read_pairs.join(', ')}"
        )
    }

    reused_fastqs = sample_rows
        .collectMany { row -> [row[1], row[2]] }
        .countBy { fastq -> fastq }
        .findAll { _fastq, count -> count > 1 }
        .keySet()
        .collect { fastq -> fastq.toString() }
        .sort()

    if (!reused_fastqs.isEmpty()) {
        error(
            'the same FASTQ file is referenced by multiple read groups: ' +
            reused_fastqs.join(', ')
        )
    }

    samples_ch = channel.fromList(sample_rows)

    // Issue #8: computed once, synchronously, directly from the fully
    // materialized samplesheet list -- before any channel operation
    // runs -- so that groupTuple() downstream can be told exactly how
    // many read groups to expect per sample via groupKey() and emit
    // each sample's merged BAM the moment its own read groups are all
    // mapped, rather than waiting for every sample's mapping to finish.
    read_group_counts_by_sample = sample_rows
        .collect { row -> row[0].id }
        .countBy { sample_id -> sample_id }

    reference_meta = [
        id       : params.reference_id,
        name     : params.reference_name,
        accession: params.reference_accession,
        species  : params.reference_species,
        cultivar : params.reference_cultivar
    ]

    reference_fasta = file(
        params.reference_fasta,
        checkIfExists: true
    )

    if (params.reference_fai) {
        reference_fai = file(
            params.reference_fai,
            checkIfExists: true
        )
        expected_fai_name = "${reference_fasta.name}.fai"

        if (reference_fai.name != expected_fai_name) {
            error(
                "reference_fai must be named ${expected_fai_name}; " +
                "found ${reference_fai.name}"
            )
        }
    }

    if (params.reference_dict) {
        reference_dict = file(
            params.reference_dict,
            checkIfExists: true
        )
        expected_dict_name = "${reference_fasta.baseName}.dict"

        if (reference_dict.name != expected_dict_name) {
            error(
                "reference_dict must be named ${expected_dict_name}; " +
                "found ${reference_dict.name}"
            )
        }
    }

    if (params.bwa_index_prefix) {
        bwa_index_prefix = file(
            params.bwa_index_prefix,
            checkIfExists: false
        )

        if (bwa_index_prefix.name != reference_fasta.name) {
            error(
                'bwa_index_prefix basename must match reference_fasta; ' +
                "expected ${reference_fasta.name}, " +
                "found ${bwa_index_prefix.name}"
            )
        }

        bwa_index_suffixes = [
            '.0123',
            '.amb',
            '.ann',
            '.bwt.2bit.64',
            '.pac'
        ]

        bwa_index_suffixes.each { suffix ->
            file(
                "${params.bwa_index_prefix}${suffix}",
                checkIfExists: true
            )
        }
    }

    reference_ch = channel.value(
        tuple(
            reference_meta,
            reference_fasta
        )
    )

    ADZUKI_SNP_PIPELINE(
        samples_ch,
        reference_ch,
        read_group_counts_by_sample
    )
}
