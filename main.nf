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

    reference_meta = [
        id       : params.reference_id,
        name     : params.reference_name,
        accession: params.reference_accession,
        species  : params.reference_species,
        cultivar : params.reference_cultivar
    ]

    reference_ch = channel.value(
        tuple(
            reference_meta,
            file(params.reference_fasta, checkIfExists: true)
        )
    )

    ADZUKI_SNP_PIPELINE(
        samples_ch,
        reference_ch
    )
}
