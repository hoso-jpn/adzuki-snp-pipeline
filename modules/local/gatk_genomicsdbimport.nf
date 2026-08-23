// Issue #11: two independent fixes, both found from the same real-data
// evidence (Issue #8 Phase 5's Seedcore-01 run of a single SRR29909135
// gVCF against the real Longxiaodou 4 reference, which failed here).
//
// 1. List-safety: Nextflow's `path` input qualifier silently unwraps a
//    single-element List into a bare scalar Path/File when exactly one
//    file resolves for that channel slot -- this is documented Nextflow
//    behavior, not specific to this pipeline. `gvcfs`/`gvcf_indexes` are
//    always Lists in every synthetic fixture (>= 2 samples), so this
//    never surfaced until a real single-accession run hit it: calling
//    `.size()` on the resulting scalar File returned its *byte length*,
//    not "1" -- the real run's error message ("1176577035 gVCFs, 642455
//    indexes") is literally the real gVCF's and its .tbi index's file
//    sizes. Both inputs are coerced to Lists before anything else in the
//    script block touches them.
// 2. Native-memory headroom: GenomicsDBImport's actual memory footprint
//    is JVM heap (-Xmx) *plus* a separate native/C++ TileDB storage
//    layer (buffers, mmap regions, thread stacks) that the JVM heap
//    setting does not account for at all. The previous
//    `task.memory.toGiga() - 1` formula reserved a fixed 1 GiB
//    regardless of task.memory's actual size -- at the 16 GiB
//    'process_high' allocation this leaves the native layer only
//    1/16 = 6.25% of the task's memory (Xmx = 93.75%), not the >= 20%
//    headroom this Issue requires. Xmx is now capped at a flat 80% of
//    task.memory, computed in MiB (matching the precedent set in
//    BWA_MEM2_MEM_SORT for Issue #8) rather than rounded to whole GiB,
//    so small task.memory values cannot round in a way that overshoots
//    80%. The remaining >= 20% for the native layer falls out of the
//    same ceiling automatically; it is not asserted as a separate
//    number.
//
// `process_genomicsdb` is a dedicated resource label so this process's
// memory can be tuned independently of the rest of 'process_high' as
// real multi-sample cohort data (Issue #26) reveals GenomicsDBImport's
// actual footprint; its cpus/memory/time are carried over from
// 'process_high' unchanged -- see nextflow.config -- since no real
// multi-sample measurement yet justifies different values.
process GATK_GENOMICSDBIMPORT {
    tag "${interval_meta.id}"
    label 'process_genomicsdb'

    container 'broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb'

    input:
    tuple val(interval_meta), val(interval)
    path(gvcfs)
    path(gvcf_indexes)

    output:
    tuple(
        val(interval_meta),
        val(interval),
        path("${interval_meta.id}.genomicsdb"),
        emit: genomicsdb
    )

    script:
    gvcf_list = gvcfs instanceof List ? gvcfs : [gvcfs]
    gvcf_index_list = gvcf_indexes instanceof List ? gvcf_indexes : [gvcf_indexes]

    if (gvcf_list.size() != gvcf_index_list.size()) {
        error(
            'the number of gVCFs and indexes must match: ' +
            "${gvcf_list.size()} gVCFs, " +
            "${gvcf_index_list.size()} indexes"
        )
    }

    variant_arguments = gvcf_list
        .collect { gvcf -> "--variant ${gvcf}" }
        .join(" \\\n        ")

    // See BWA_MEM2_MEM_SORT (Issue #8) for the same toMega()-based,
    // fail-fast-rather-than-clamp precedent. task.memory is re-read
    // fresh on every OOM-retry attempt (errorStrategy scales it as
    // `16.GB * task.attempt` via 'process_genomicsdb'), so this ratio
    // is recomputed correctly on every attempt rather than holding a
    // stale value from the first.
    total_memory_mib = task.memory.toMega()
    xmx_mib = Math.floor(total_memory_mib * 0.8).intValue()

    if (xmx_mib < 1) {
        error(
            "GATK_GENOMICSDBIMPORT (${interval_meta.id}) needs more memory: " +
            "task.memory=${task.memory} yields an 80% JVM heap ceiling of " +
            "${xmx_mib} MiB, below the 1 MiB floor. Increase the " +
            "'process_genomicsdb' resource label's memory in " +
            'nextflow.config or the active profile config.'
        )
    }

    """
    gatk --java-options "-Xmx${xmx_mib}m" GenomicsDBImport \
        ${variant_arguments} \
        --genomicsdb-workspace-path ${interval_meta.id}.genomicsdb \
        --intervals '${interval}' \
        --reader-threads ${task.cpus} \
        --batch-size ${params.genomicsdb_batch_size} \
        --tmp-dir .
    """
}
