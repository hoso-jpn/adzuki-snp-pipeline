def sampleMapQuote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process BUILD_SAMPLE_NAME_MAP {
    tag 'cohort'
    label 'process_low'
    container params.containers.python

    input:
    val(entries)
    path(gvcfs)
    path(indexes)

    output:
    path('cohort.sample_name_map.tsv'), emit: sample_map
    path('cohort.joint_inputs.json'), emit: summary
    path(gvcfs), emit: gvcfs
    path(indexes), emit: indexes
    val(task.container), emit: container_id

    script:
    def entry_json = groovy.json.JsonOutput.toJson(entries)
    """
    printf '%s' ${sampleMapQuote(entry_json)} > joint_entries.json
    python3 ${sampleMapQuote("${projectDir}/bin/build_sample_name_map.py")} --entries joint_entries.json \
        --output cohort.sample_name_map.tsv --summary cohort.joint_inputs.json
    """
}
