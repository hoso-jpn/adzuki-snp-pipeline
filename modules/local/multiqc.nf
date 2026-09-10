process MULTIQC {
    tag 'cohort'
    label 'process_low'

    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1@sha256:b65e3fe879df27b92334dda0fd987a6e21bdee09a2848551d4f287099a93b7ac'

    // Issue #51: container tasks now run as the launching host user
    // (nextflow.config's `docker.runOptions`), and an unprivileged user
    // cannot create `/multiqc` at the container filesystem root the way
    // the container's root previously could -- the task fails at its very
    // first line with `mkdir: cannot create directory '/multiqc'`.
    //
    // The fixed container-absolute staging root is not incidental: Issue
    // #38 chose it so that the source paths MultiQC records in
    // multiqc_data.json and multiqc_sources.json are container paths and
    // never this host's work directory (see docs/multiqc.md, and the
    // `/multiqc/input/` assertions in tests/pipeline). Relocating it would
    // change the published provenance those artifacts carry, so it is kept
    // exactly as it is and simply given a writable mount. A tmpfs rather
    // than a bind of the task's own work directory: binding the work
    // directory onto `/multiqc` makes the staging root and the task
    // directory the same directory, which turns the script's own
    // `cp multiqc_config.yaml /multiqc/multiqc_config.yaml` into a
    // same-file copy error. What lands here is MultiQC's staging and
    // report output, which the script then copies back into the task
    // directory; it is therefore charged to this task's memory budget
    // rather than to disk, and its size has been observed only at the
    // synthetic fixture's scale.
    containerOptions '--tmpfs /multiqc:rw'

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

    // Issue #42: this process's *effective* container -- Nextflow's own
    // task.container, resolved after any withName/alias/fully-qualified-
    // selector/profile override on top of the `container` directive above
    // -- so the run-level provenance manifest records what this task
    // actually ran in rather than a default the pipeline assumed. See
    // workflows/adzuki_snp_pipeline.nf for the canonical process key this
    // invocation is recorded under, and docs/run_manifest_data_contract.md
    // for the schema v2 contract.
    val(task.container), emit: container_id

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
