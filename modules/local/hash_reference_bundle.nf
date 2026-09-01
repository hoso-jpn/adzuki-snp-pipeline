process HASH_REFERENCE_BUNDLE {
    tag "${meta.id}"
    label 'process_low'

    // See modules/local/summarize_variant_qc.nf for why the full
    // (non-"-slim") Python image is required.
    container params.containers.python

    input:
    tuple val(meta), path(fasta)
    tuple val(fai_meta), path(fai)
    tuple val(dict_meta), path(sequence_dictionary)
    tuple val(index_meta), path(bwa_indexes)

    output:
    tuple(
        val(meta),
        path("${meta.id}.reference_provenance.tsv"),
        emit: provenance
    )
    // Issue #42: see modules/local/hash_input_fastqs.nf for why every
    // containerized process emits its own effective task.container.
    val(task.container), emit: container_id

    script:
    // Issue #42: the BWA-MEM2 index files are hashed here alongside the
    // FASTA/FAI/dict because they are a real mapping input, and the
    // pipeline accepts a *prebuilt* index (--bwa_index_prefix) just as
    // readily as one it builds itself. A stale or mismatched prebuilt
    // index changes alignment while leaving reference_id and the FASTA
    // checksum identical, so a manifest that recorded only the FASTA
    // could not tell the two runs apart. `bwa_indexes` is whatever
    // actually reached BWA_MEM2_MEM_SORT, generated or prebuilt alike.
    def index_args = bwa_indexes
        .collect { index -> "--bwa-index ${index}" }
        .join(' \\\n        ')

    """
    hash_reference_bundle.py \
        --fasta ${fasta} \
        --fai ${fai} \
        --dict ${sequence_dictionary} \
        ${index_args} \
        --output ${meta.id}.reference_provenance.tsv
    """
}
