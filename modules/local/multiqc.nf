process MULTIQC {
    tag 'cohort'
    label 'process_low'

    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1@sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac'

    publishDir(
        "${params.outdir}/qc/multiqc",
        mode: 'copy'
    )

    input:
    path raw_fastqc_zips
    path trimmed_fastqc_zips
    path fastp_jsons
    path markduplicates_metrics
    path samtools_flagstats
    path samtools_stats
    path samtools_idxstats
    path multiqc_config, name: 'multiqc_config.yaml'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data
    path 'multiqc_config.yaml', emit: config
    path 'multiqc_version.txt', emit: version

    script:
    def asFileList = { value ->
        value == null
            ? []
            : (value instanceof List ? value : [value])
    }
    def categories = [
        raw_fastqc_zips: asFileList.call(raw_fastqc_zips),
        trimmed_fastqc_zips: asFileList.call(trimmed_fastqc_zips),
        fastp_jsons: asFileList.call(fastp_jsons),
        markduplicates_metrics: asFileList.call(markduplicates_metrics),
        samtools_flagstats: asFileList.call(samtools_flagstats),
        samtools_stats: asFileList.call(samtools_stats),
        samtools_idxstats: asFileList.call(samtools_idxstats),
    ]
    def empty_categories = categories
        .findAll { _name, files -> files == null || files.isEmpty() }
        .keySet()

    if (!empty_categories.isEmpty()) {
        error(
            'MULTIQC requires at least one artifact in every QC category; empty: ' +
            empty_categories.join(', ')
        )
    }

    def input_files = categories.values()
        .collectMany { files -> files }
        .collect { report -> report.getName() }
    def quoted_input_files = input_files
        .collect { report -> "'/multiqc/input/${report}'" }
        .join(' ')
    def stage_links = input_files
        .collect { report ->
            "ln -s \"\${PWD}/${report}\" '/multiqc/input/${report}'"
        }
        .join('\n')

    """
    mkdir -p /multiqc/input /multiqc/output /multiqc/scratch
    ${stage_links}
    cp multiqc_config.yaml /multiqc/multiqc_config.yaml
    printf '%s\\n' ${quoted_input_files} > /multiqc/multiqc_inputs.txt

    multiqc --version > multiqc_version.txt

    (
        cd /multiqc
        TMPDIR=/multiqc/scratch multiqc \
            --config multiqc_config.yaml \
            --filename multiqc_report.html \
            --outdir /multiqc/output \
            --force \
            --file-list multiqc_inputs.txt
    )

    cp /multiqc/output/multiqc_report.html .
    cp -R /multiqc/output/multiqc_data .
    """
}
