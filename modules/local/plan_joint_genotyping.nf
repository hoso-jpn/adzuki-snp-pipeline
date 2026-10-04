def planQuote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process PLAN_JOINT_GENOTYPING {
    tag 'cohort'
    label 'process_low'
    container params.containers.python

    input:
    tuple val(meta), path(fai)
    path(samples)
    path(input_provenance)
    tuple val(reference_meta), path(reference_manifest)
    val(pipeline_sha)
    tuple val(review_enabled), path(review_file)
    val(output_free_bytes)

    output:
    path('cohort.joint_plan.json'), emit: plan
    path('cohort.joint_intervals.tsv'), emit: intervals
    val(task.container), emit: container_id

    script:
    def review_arg = review_enabled ? "--review ${planQuote(review_file)}" : ''
    def provenance_paths = input_provenance instanceof List ? input_provenance : [input_provenance]
    def provenance_args = provenance_paths.collect { path -> planQuote(path) }.join(' ')
    """
    python3 ${planQuote("${projectDir}/bin/plan_joint_genotyping.py")} \
        --fai ${planQuote(fai)} --samples ${planQuote(samples)} \
        --input-provenance ${provenance_args} \
        --reference-manifest ${planQuote(reference_manifest)} \
        --pipeline-sha ${planQuote(pipeline_sha)} --ploidy ${params.sample_ploidy} \
        --window ${params.joint_interval_size_bp} \
        --task-memory-gib ${params.joint_genotype_memory_gib} \
        --concurrency ${params.joint_max_concurrency} \
        --memory-budget-gib ${params.joint_launch_memory_gib} \
        --output-free-bytes ${output_free_bytes} ${review_arg} \
        --output cohort.joint_plan.json --intervals cohort.joint_intervals.tsv \
        --require-eligible
    """
}
