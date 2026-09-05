// Issue #8: BWA-MEM2's alignment stream is piped directly into
// `samtools sort` instead of being written to an uncompressed
// intermediate `.sam` file first. For real WGS data an uncompressed SAM
// can be many times the size of the final coordinate-sorted BAM, which
// risked filling the task's work-directory disk before sorting even
// began. `set -o pipefail` (bash's default `-e` alone does not check
// exit codes of earlier commands in a pipe) makes a failure in either
// tool -- a segfault, truncated/corrupt input, disk full -- fail the
// whole task, rather than silently succeeding on `samtools sort`'s own
// exit code over a truncated alignment stream.
//
// Combined image, not two separate bwa-mem2/samtools containers: the
// two tools run concurrently inside one pipe and must share a single
// process's CPU/memory budget, which a container boundary would hide
// from Nextflow's resource accounting. This is a Seqera Wave-built
// multi-package image, verified directly (`docker run ... which`,
// `bwa-mem2 version`, `samtools --version`) to bundle bwa-mem2 2.2.1 --
// the same real binary version already running under this repository's
// previously separate, pinned `bwa-mem2:2.3` tag; see the BWA-MEM2
// section of README for the documented tag/binary version mismatch --
// and samtools 1.22.1. No existing mulled/Wave image combining
// bwa-mem2 with samtools 1.24 (the version pinned everywhere else in
// this pipeline: SAMTOOLS_MERGE, SAMTOOLS_INDEX, SAMTOOLS_QC,
// SAMTOOLS_FAIDX) was found; this 1.22.1-vs-1.24 gap is a deliberate,
// documented limitation of this specific choice, not an oversight --
// see README. The BAM produced here has been confirmed, for this
// pipeline's synthetic fixture, to carry the expected coordinate-sort
// header, and the generated-index and prebuilt-index paths -- both
// sorted by this same samtools 1.22.1 -- produce a byte-identical raw
// cohort VCF once fed through the 1.24-pinned downstream stages; this
// is not a comparison of samtools 1.22.1's own output against
// samtools 1.24's. That is not evidence that 1.22.1 and 1.24 are
// interchangeable in general (BGZF interoperability does not imply
// identical sort-tie-break ordering, compression, or bug fixes
// between versions; the SAM spec leaves same-RNAME/POS record order
// unspecified). Any behavioral difference on real data between these
// two samtools versions themselves is unverified and left to Phase 5
// (real-reference profiling) to surface.
process BWA_MEM2_MEM_SORT {
    tag "${meta.id}:${meta.read_group_id}"
    label 'process_mapping'

    container 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113@sha256:5ebd1290d9680195817ce75915b79ae2e608834c017824b7e2bc7b141509b242'

    publishDir(
        "${params.outdir}/logs/mapping",
        mode: 'copy',
        pattern: '*.bwa-mem2.log'
    )

    input:
    tuple val(meta), path(read1), path(read2)
    tuple val(reference_meta), path(fasta)
    tuple val(index_meta), path(indexes)

    output:
    tuple(
        val(meta),
        path("${meta.read_group_id}.sorted.bam"),
        emit: bam
    )
    tuple(
        val(meta),
        path("${meta.read_group_id}.bwa-mem2.log"),
        emit: log
    )

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
    platform_unit = meta.platform_unit
        ? "\\tPU:${meta.platform_unit}"
        : ''
    read_group = (
        "@RG\\tID:${meta.read_group_id}" +
        "\\tSM:${meta.id}" +
        "\\tLB:${meta.library_id}" +
        "\\tPL:${meta.platform}" +
        platform_unit
    )

    // Neither tool receives task.cpus directly: they run concurrently
    // in the same pipe and would otherwise both claim the task's full
    // CPU share, oversubscribing it. bwa-mem2's alignment/scoring work
    // dominates the pair's total CPU cost and scales close to linearly
    // with threads, so it gets the majority (80%); samtools sort's
    // speedup from additional threads is smaller past a handful, so it
    // gets the remainder. This 80/20 split is a documented starting
    // assumption, not a measurement -- Issue #8 Phase 5 (real-reference
    // profiling) is expected to confirm or retune it once real hardware
    // and a real reference are available.
    //
    // PR #25 review (P1-1): a previous version of this script clamped
    // both thread counts and the per-thread memory share to a floor of
    // 1 instead of validating task.cpus/task.memory up front. On a
    // small enough resource override, those floors silently violated
    // the very budget this formula exists to respect -- e.g.
    // task.cpus=1 gave bwa_threads=1 *and* sort_threads=1 (2 threads
    // sharing a 1-CPU task), and task.memory=1024 MiB gave
    // bwa_share_mib=512 + os_overhead_mib=512 + a memory-clamped
    // sort_share_mib=1, summing to 1025 MiB against a 1024 MiB budget.
    // Failing fast here, with a diagnosable message, is safer than a
    // clamp that hides an already-broken allocation.
    if (task.cpus < 2) {
        error(
            "BWA_MEM2_MEM_SORT (${meta.id}:${meta.read_group_id}) needs " +
            "task.cpus >= 2 -- bwa-mem2 and samtools sort run " +
            "concurrently in the same pipe and each need at least one " +
            "thread of their own. Got task.cpus=${task.cpus} from the " +
            "'process_mapping' resource label; increase its cpus in " +
            'nextflow.config or the active profile config.'
        )
    }

    sort_threads = Math.max(1, (task.cpus * 0.2).intValue())
    bwa_threads = task.cpus - sort_threads

    // `samtools sort -m` is memory *per thread*. Omitting it defaults
    // to a fixed 768 MiB/thread regardless of what the task actually
    // has available, which on a small allocation risks `sort_threads`
    // threads collectively claiming more than task.memory and
    // triggering exactly the OOM-kill behavior this Issue hardens
    // against. The budget is split with documented, non-measured
    // assumptions: half of task.memory is reserved for bwa-mem2 itself
    // (its footprint is dominated by the reference index, which scales
    // with genome size in a way this formula deliberately does not try
    // to predict ahead of Phase 5's real measurement), and a fixed
    // 512 MiB is reserved for OS/tool overhead (page cache pressure,
    // allocator fragmentation, samtools' own non-buffer bookkeeping).
    // Whatever remains is divided evenly across sort's own threads; if
    // that would be below 1 MiB/thread, task.memory is too small for
    // this task.cpus and the task fails fast rather than silently
    // running sort under-provisioned.
    total_memory_mib = task.memory.toMega()
    os_overhead_mib = 512
    bwa_share_mib = Math.max(1, (total_memory_mib * 0.5).intValue())
    sort_share_mib = total_memory_mib - bwa_share_mib - os_overhead_mib

    if (sort_share_mib < sort_threads) {
        error(
            "BWA_MEM2_MEM_SORT (${meta.id}:${meta.read_group_id}) needs " +
            "more memory: task.memory=${task.memory} splits into " +
            "bwa_share_mib=${bwa_share_mib} + os_overhead_mib=${os_overhead_mib}, " +
            "leaving only ${sort_share_mib} MiB for ${sort_threads} " +
            "samtools-sort thread(s) -- below the 1 MiB/thread floor. " +
            "Increase the 'process_mapping' resource label's memory in " +
            'nextflow.config or the active profile config.'
        )
    }

    sort_mem_per_thread_mib = (sort_share_mib / sort_threads).intValue()

    """
    set -o pipefail

    bwa-mem2 mem \
        -t ${bwa_threads} \
        -R '${read_group}' \
        ${fasta} \
        ${read1} \
        ${read2} \
        2> ${meta.read_group_id}.bwa-mem2.log \
        | samtools sort \
            -@ ${sort_threads} \
            -m ${sort_mem_per_thread_mib}M \
            -o ${meta.read_group_id}.sorted.bam \
            -
    """
}
