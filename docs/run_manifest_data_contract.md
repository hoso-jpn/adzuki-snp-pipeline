# Run manifest data contract

`bin/build_run_manifest.py` produces this pipeline's whole-run
provenance record: which samples went in, which reference they were
aligned to, which container each process actually ran in, which
parameters were in force, what the accounting came out to, and the
checksums of the deliverables that came out.

There are two schema versions, and they are two different documents
produced two different ways. This file describes both, because both
exist in this repository's published record.

| | schema v1 | schema v2 |
| --- | --- | --- |
| Produced by | `build_run_manifest.py --mode legacy-v1`, run by hand after a run | `BUILD_RUN_MANIFEST`, an ordinary process in the Nextflow DAG |
| Introduced in | Issue #26 | Issue #42 |
| `containers` | four tool-level entries, hand-passed | one entry per *process actually executed*, read from that process's own effective `task.container` |
| `git_commit` | nullable | required, full 40-character SHA |
| `sample_accounting` | absent | present |
| Where it lives | `docs/real_cohort_*_manifest.json` (historical runs) | `<outdir>/provenance/<cohort_id>.run_manifest.json` (every run) |

Existing schema v1 documents are **not** regenerated or rewritten. They
remain valid records of the runs that produced them, and legacy mode
remains the default so the invocation that produced them still works.

## Why the DAG-generated manifest is schema v2, not an extension of v1

Schema v1's `containers` had four keys — `bwa_mem2`, `samtools`,
`gatk`, `python` — each a container reference the caller passed on the
command line, copied by hand from a `container` directive. That shape
cannot describe this pipeline: FastQC, fastp, MultiQC and bcftools each
run in their own image, and any process can be overridden individually
by a `withName` selector, a process alias, a fully-qualified selector, a
profile, or an external config. Two aliases of one module
(`GATK_VARIANTFILTRATION` and `GATK_VARIANTFILTRATION_GS`) can
legitimately run in different images in the same run.

Keeping the field name while changing what it means would silently
change what every already-published v1 manifest appears to claim. So the
field changed shape and the version changed with it.

## Schema v2

```json
{
  "schema_version": 2,
  "run_id": "20260901T001535Z-29f8fdc4",
  "generated_at": "2026-09-01T00:15:35Z",
  "cohort_id": "cohort",
  "pipeline_version": "0.2.0",
  "git_commit": "402f437b4a4f00023e0571bbe6d73dab3a83fadf",
  "nextflow_version": "26.04.6",
  "containers": { "fastqc_raw": "...", "multiqc": "...", "...": "..." },
  "reference": {
    "reference_id": "synthetic-adzuki-v1",
    "reference_name": "Synthetic adzuki test reference",
    "reference_species": "Vigna angularis",
    "reference_cultivar": "synthetic",
    "reference_accession": "",
    "fasta": { "filename": "synthetic.fa", "checksum": "sha256:..." },
    "fai":   { "filename": "synthetic.fa.fai", "checksum": "sha256:..." },
    "dict":  { "filename": "synthetic.dict", "checksum": "sha256:..." },
    "bwa_index": [ { "filename": "synthetic.fa.0123", "checksum": "sha256:..." }, "..." ]
  },
  "parameters": { "sample_ploidy": 2, "enable_gs_panel": true, "...": "..." },
  "samples": [
    {
      "sample_id": "sample_a",
      "read_group_id": "sample_a_L001",
      "library_id": "lib_a",
      "platform": "ILLUMINA",
      "platform_unit": "flowcell1.L001.ATCACG",
      "fastq_1": { "filename": "sample_a_L001_R1.fastq.gz", "checksum": "sha256:..." },
      "fastq_2": { "filename": "sample_a_L001_R2.fastq.gz", "checksum": "sha256:..." }
    }
  ],
  "cohort_accounting": { "number_of_samples": "2", "sample_names": "sample_a,sample_b", "...": "..." },
  "sample_accounting": [ { "sample": "sample_a", "missingness_rate": "0.000000", "...": "..." } ],
  "variant_type_accounting": { "raw_all_records": "2", "...": "..." },
  "gs_panel": { "schema_version": 2, "run_id": "...", "cohort_id": "cohort", "panel_status": "empty", "manifest_hash": "sha256:..." },
  "checksums": { "cohort.raw.vcf.gz": "sha256:...", "...": "..." },
  "manifest_hash": "sha256:..."
}
```

### Output path

`<outdir>/provenance/<cohort_id>.run_manifest.json` — for the default
cohort, `results/provenance/cohort.run_manifest.json`. Its own directory
rather than inside a lineage's output tree, because it describes the
whole run, including which lineages ran at all.

### `containers`: one entry per executed process

Each key is the canonical, lowercase process key the workflow assigns
that invocation; each value is Nextflow's own `task.container` for that
process, captured from inside the task after every
`withName`/alias/fully-qualified-selector/profile override has already
been resolved. This is the same mechanism Issue #52 established for the
GS panel manifest, applied to the whole run.

Three properties follow, and all three are deliberate:

* **Aliases stay separate.** `fastqc_raw` and `fastqc_trimmed` are two
  entries, as are `gatk_variantfiltration` / `gatk_variantfiltration_gs`
  and `gatk_selectpassvariants` / `gatk_selectpassvariants_gs`. They are
  independently overridable, so collapsing them by module would record a
  value that is false for one of them.
* **One process, one container.** A process that ran 300 tasks
  contributes 300 identical rows, which collapse to one entry. If two
  tasks of the same process report *different* images, the manifest is
  not written and the run fails: no single value would be a true
  statement about that process, and picking one silently would be worse
  than stopping. (A future need for per-task dynamic containers would be
  a schema change, not a quiet relaxation of this.)
* **Only what ran.** A process that did not execute has no entry.
  `samtools_faidx`, `gatk_create_sequence_dictionary` and
  `bwa_mem2_index` are absent when the run was given a prebuilt
  reference bundle; the eight GS-lineage keys are absent when
  `enable_gs_panel = false`. The field lists processes that touched this
  data, not processes the pipeline could have run.

The canonical keys are assigned explicitly in
`workflows/adzuki_snp_pipeline.nf` rather than derived from a runtime
value such as `task.process`. A published schema should not depend on
how a given Nextflow version happens to spell a qualified process name,
and a module-derived name could not distinguish the aliases above.

`build_run_manifest` is itself one of the keys: the process that writes
the manifest runs in a container like any other, and reads its own
`task.container` from inside its own script.

### `git_commit` and `nextflow_version`

`git_commit` is a full 40-character SHA and is required. Nextflow's
`workflow.commitId` is populated only when Nextflow itself pulled a
git-hosted pipeline (`nextflow run owner/repo`); for this repository's
own documented `nextflow run .` invocation it is null (measured on
Nextflow 26.04.6, and the same under nf-test). `main.nf`'s
`resolvePipelineCommit()` therefore falls back to reading the commit
from the project directory's git checkout, and the run fails if it
cannot resolve a full SHA. Only the SHA is used — not the project
directory, not the git command's output.

`nextflow_version` is the engine version that actually ran, recorded
separately from any container: Nextflow runs on the host, outside every
pinned image, and its own runtime semantics have previously been the
root cause of a real production failure here (Issue #11).

### `samples`: input FASTQ provenance

One entry per read group, in samplesheet order, carrying the
samplesheet's own metadata plus each FASTQ's basename and SHA-256.

The hashing happens in `HASH_INPUT_FASTQS`, one task per read group,
not in the manifest process. Handing every raw FASTQ in the cohort to a
single terminal task would re-stage the entire input dataset — at this
pipeline's 327-sample target, hundreds of gigabytes — purely to write a
JSON file.

Nextflow makes no ordering promise across parallel tasks, so each row
carries its zero-based samplesheet position, and the manifest builder
sorts on it and then drops it. That rank never appears in the document:
it means nothing outside the run that produced it, whereas "the order
the samplesheet listed them in" is something a reader can check against
their own input.

### `reference`

Identity (`reference_id`, `reference_name`, `reference_species`,
`reference_cultivar`, `reference_accession`) plus the checksum of every
reference file the run actually consumed: the FASTA, its FAI, its
sequence dictionary, and the five BWA-MEM2 index files — whether this
run generated the index or was given a prebuilt one. A stale or
mismatched prebuilt index changes alignment while leaving `reference_id`
and the FASTA checksum identical, so recording only the FASTA could not
tell two such runs apart.

### `cohort_accounting`, `sample_accounting`, `variant_type_accounting`

All three are read from artifacts the pipeline already published, never
recomputed:

* `cohort_accounting` — `SUMMARIZE_VARIANT_QC`'s raw/all
  `variant_qc.tsv`;
* `sample_accounting` — the `sample_qc.tsv` from that same invocation,
  one row per sample in that file's own order;
* `variant_type_accounting` — `RECONCILE_VARIANT_TYPE_COUNTS`'s TSV.

A second implementation of the same numbers is a second chance to
disagree with the QC the pipeline published, so there isn't one.

The builder cross-checks them before recording anything: the cohort's
`number_of_samples` and `sample_names` must match the per-sample rows
exactly, in order; every file must be for this cohort; and the QC files
must be the raw/all ones they claim to be. A provenance record that
embedded two contradictory statements about the same run without
complaint would be worse than none, so a mismatch fails the run.

### `gs_panel`

When the GS panel ran, this embeds the GS panel manifest's own summary —
`schema_version`, `run_id`, `cohort_id`, `panel_status`,
`manifest_hash` — as a pointer. The GS manifest already checksums its
own artifacts (see `docs/gs_panel_data_contract.md`); nothing is
re-derived or duplicated here. A GS manifest for a different cohort
fails the run rather than being attached to this one.

When `enable_gs_panel = false`, `gs_panel` is explicit `null` — a
present field with a null value, so a reader can distinguish "this run
had no GS panel" from "this manifest predates the field" — and none of
the eight GS process keys appears in `containers`.

### `checksums`

The run's own deliverables: each sample's gVCF, the cohort raw VCF, and
the primary SNP and indel PASS VCFs. This is the same artifact set the
historical schema v1 manifests recorded.

They are hashed by `HASH_RUN_ARTIFACTS`, invoked once per artifact group
(per sample for gVCFs, once per cohort-level VCF), so hashing fans out
with the rest of the pipeline instead of funnelling every deliverable
through the terminal task. Keys are basenames; two different artifacts
sharing a basename fails rather than silently recording one of them.

## Failure is a run failure

`BUILD_RUN_MANIFEST` is an ordinary required process in the DAG, not a
`workflow.onComplete` handler, an `afterScript`, or a post-run script.
An `onComplete` handler cannot fail the run it is reporting on, which is
exactly the outcome this contract exists to prevent: *"the analysis
succeeded but there is no provenance record"* is not a successful run.
The process has no `errorStrategy 'ignore'`, and the manifest is written
atomically — a failure leaves neither a partial document nor a stale
temp file.

Because the manifest consumes the aggregated container provenance of
every executed process, it is also a genuine dataflow barrier: it cannot
start until every process that contributes to it has finished, including
the MultiQC side branch, which produces nothing the manifest's *content*
depends on and would otherwise be free to still be running.

## Privacy contract

A published manifest records basenames, checksums, scientifically
meaningful identifiers and numbers. It must never record:

* an absolute host path, a home directory, or a `file://` URI;
* a user name;
* `workDir`, `launchDir`, `projectDir`, or the command line;
* a private or loopback network location;
* a credential, token or password;
* raw sequence, `.nextflow.log`, or raw trace output.

This is achieved primarily by never passing host-specific values into
the payload at all — every file arrives as a basename plus a checksum,
and no working directory or command line is an input to the builder.
`manifest_utils.assert_no_host_metadata()` is the backstop: it walks the
finished document and refuses to write it if any such value reached it
through some field nobody thought about.

The check is deliberately narrow. It does not try to detect user names
by guessing which strings look like personal names: a cohort that
legitimately contains a sample called `alice` is science, and a privacy
check that rejected it would push callers toward renaming their data to
satisfy a linter.

## Legacy schema v1 (standalone)

`build_run_manifest.py` with no `--mode` (or `--mode legacy-v1`) behaves
exactly as it did before Issue #42: run by hand against a completed
run's published artifacts, with `--samplesheet`, the four
`--*-container` flags, `--reference-fasta/fai/dict` and repeatable
`--checksum-file`. It writes `schema_version: 1`.

This mode exists because it produced the historical real-cohort
manifests referenced elsewhere in this repository's documentation
(`docs/real_cohort_e2e_run_manifest.json`,
`docs/real_cohort_scale_validation_10sample_manifest.json`,
`docs/real_cohort_scale_validation_20sample_manifest.json`). Those
documents are unchanged by Issue #42 and remain valid schema v1 records
of the runs that produced them.
