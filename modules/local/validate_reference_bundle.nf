def referenceShellQuote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process VALIDATE_REFERENCE_BUNDLE {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.python

    input:
    tuple val(meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(dictionary)
    tuple val(index_meta), path(indexes)
    tuple val(origin), path(proof, stageAs: 'input_reference_proof'), path(tool_version)

    output:
    tuple val(meta), path("${meta.id}.reference_bundle.json"), emit: manifest
    tuple val(fai_meta), path(fai), emit: fai
    tuple val(dict_meta), path(dictionary), emit: dict
    tuple val(index_meta), path(indexes), emit: indexes
    val(task.container), emit: container_id

    script:
    def index_args = ([indexes].flatten())
        .collect { index -> "--bwa-index ${referenceShellQuote(index)}" }.join(' ')
    def origin_args = origin.mode == 'generated'
        ? "--receipt ${referenceShellQuote(proof)} --tool-version ${referenceShellQuote(tool_version)} --container ${referenceShellQuote(origin.container)}"
        : "--prebuilt-manifest ${referenceShellQuote(proof)}"
    """
    python3 ${projectDir}/bin/validate_reference_bundle.py \
        --fasta ${referenceShellQuote(fasta)} \
        --fai ${referenceShellQuote(fai)} \
        --dictionary ${referenceShellQuote(dictionary)} \
        ${index_args} ${origin_args} \
        --output ${meta.id}.reference_bundle.json
    """
}
