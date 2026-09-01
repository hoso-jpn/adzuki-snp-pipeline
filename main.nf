#!/usr/bin/env nextflow

include {
    validateParameters
    paramsSummaryLog
    samplesheetToList
} from 'plugin/nf-schema'

include {
    ADZUKI_SNP_PIPELINE
} from './workflows/adzuki_snp_pipeline'

// Issue #42: the full 40-character commit this run executed from, and
// proof that the working tree it ran from is that commit.
//
// The run manifest is a provenance artifact the pipeline produces for
// itself, so "which code produced this result" is not allowed to be
// unknown -- or, worse, confidently wrong. Two things have to hold, and
// both are checked before any analysis process starts:
//
//   * a full 40-character commit can be resolved at all. Nextflow's own
//     `workflow.commitId` is populated only when Nextflow pulled the
//     pipeline from a git host (`nextflow run owner/repo`); for this
//     repository's own documented `nextflow run .` invocation it is null
//     (measured on Nextflow 26.04.6, and the same under nf-test), so the
//     commit is read from the project directory's git checkout instead.
//   * the working tree is *clean*. This is the part a HEAD lookup alone
//     cannot give: with uncommitted changes the code that actually ran is
//     not the code at HEAD, and a manifest recording that SHA would be
//     false provenance -- the most damaging kind, because it looks
//     authoritative and is checkable by nobody.
//
// `git status --porcelain` covers all three ways a tree can differ from
// its commit -- unstaged changes to tracked files, staged changes, and
// untracked files -- in one check, and honors .gitignore, so a run's own
// `work/`, `results/` and `.nf-test/` output does not count as a
// modification of the pipeline.
//
// Nothing host-specific reaches the error messages: no project path, no
// list of dirty files, no user name, and none of git's own stderr. A
// provenance guard that leaked the layout of the machine it guarded would
// be self-defeating.
def gitCommandOutput(String projectPath, List<String> arguments) {
    def process = (['git', '-C', projectPath] + arguments).execute()
    def standardOutput = new StringBuilder()
    def errorOutput = new StringBuilder()
    // Consume both streams while waiting: a `git status` in a repository
    // with many changes can fill a pipe buffer and deadlock a bare
    // waitFor().
    process.waitForProcessOutput(standardOutput, errorOutput)

    return [
        exit_code: process.exitValue(),
        output: standardOutput.toString().trim(),
    ]
}

def resolvePipelineCommit(projectDirectory, declaredCommitId) {
    def projectPath = projectDirectory.toString()

    def head = gitCommandOutput(projectPath, ['rev-parse', 'HEAD'])
    if (head.exit_code != 0 || !(head.output ==~ /^[0-9a-f]{40}$/)) {
        error(
            'could not resolve the git commit this pipeline is running from. ' +
            'The run-level provenance manifest (Issue #42) records the exact ' +
            'commit a result came from and will not record an unknown or ' +
            'abbreviated one, so this run stops before doing any analysis. Run ' +
            'the pipeline from a git checkout of this repository, or from a ' +
            'git-hosted revision.'
        )
    }

    // `git -C <dir>` walks *up* to the nearest enclosing repository, so a
    // project directory that is not itself under version control -- an
    // unpacked tarball that happens to sit inside some unrelated checkout,
    // say -- would otherwise resolve that repository's HEAD and record it
    // as the code that ran. Requiring the pipeline's own entry point to be
    // a tracked file of the repository whose commit is about to be
    // recorded rules that out, while still allowing the legitimate case
    // where the pipeline lives in a subdirectory of a larger repository
    // that genuinely versions it. (`mainScript` in nextflow.config.)
    def tracked = gitCommandOutput(projectPath, ['ls-files', '--error-unmatch', 'main.nf'])
    if (tracked.exit_code != 0) {
        error(
            "the pipeline's own source is not tracked by the git repository " +
            'its commit would be read from, so that commit does not identify ' +
            'the code being run. The run-level provenance manifest (Issue #42) ' +
            'records the exact commit a result came from, so this run stops ' +
            'before doing any analysis. Run the pipeline from a git checkout ' +
            'of this repository, or from a git-hosted revision.'
        )
    }

    if (declaredCommitId && declaredCommitId.toString() != head.output) {
        error(
            'the revision Nextflow reports for this run does not match the ' +
            "commit checked out in the pipeline's project directory. Refusing " +
            'to record either one as the code that produced this run.'
        )
    }

    def status = gitCommandOutput(projectPath, ['status', '--porcelain'])
    if (status.exit_code != 0) {
        error(
            "could not determine whether the pipeline's working tree matches " +
            'its git commit. The run-level provenance manifest (Issue #42) ' +
            'records that commit as the code that produced the run, so this ' +
            'run stops rather than recording a claim it cannot support.'
        )
    }

    if (!status.output.isEmpty()) {
        error(
            "the pipeline's working tree has uncommitted changes (modified, " +
            'staged, or untracked files). The run-level provenance manifest ' +
            '(Issue #42) would record the HEAD commit as the code that ' +
            'produced this run, which is not true while the working tree ' +
            'differs from it. Commit, stash, or remove the changes before a ' +
            'formal provenance run.'
        )
    }

    return head.output
}

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
            'is true (the default). The GS panel genotype encoding ' +
            '(diploid_additive_dosage_v1) is diploid-only and would ' +
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

    // Issue #42: a provenance-only view of the same rows, each tagged
    // with its zero-based samplesheet position. Built from the already
    // materialized `sample_rows` list rather than by transforming
    // samples_ch, so the scientific channel's own cardinality and
    // ordering contract is untouched: Nextflow makes no ordering promise
    // across parallel tasks, and this rank is how the run manifest
    // restores samplesheet order afterwards. The rank is a transport
    // detail -- bin/build_run_manifest.py sorts on it and drops it.
    input_provenance_rows_ch = channel.fromList(
        sample_rows.withIndex().collect { row, index ->
            tuple(row[0] + [rank: index], row[1], row[2])
        }
    )

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
        read_group_counts_by_sample,
        input_provenance_rows_ch,
        resolvePipelineCommit(workflow.projectDir, workflow.commitId)
    )
}
