# Joint Genotyping execution and recovery (Issue #69)

The production DAG now consumes the adopted #45 design: header-checked
sample-name-map, disjoint reference-ordered windows, separate small contigs,
variant-start ownership at boundaries, and reference-rank gathering. Reblock,
small-scaffold grouping and consolidation remain off. A dedicated GenotypeGVCFs
label reserves 20% of its allocation outside the JVM and limits concurrency;
it does not double memory on an OOM retry beyond the reviewed launch budget.

This is an implementation change, not a completed real-scale validation.
The public 20-sample E2E and 51-sample targeted evidence in #45 are historical
observations. The new DAG's 51-sample replay, production resource measurements
and commercial SLA are **not established** by those earlier observations.
The 327-sample decision remains **NO-GO**.

## Before reads are processed

Reference-only preparation and the #68 bundle gate run first. The planner
then records input byte counts, measured FASTQ SHA-256s bound to read-group
metadata, sample IDs, reference fingerprint, pipeline
SHA, interval plan, resource settings and a free-space snapshot. FASTQ QC,
trimming and mapping cannot start until that plan is eligible.

- Up to 20 samples: eligible for *research execution* if the other gates pass.
  Unreviewed references are limited to the bundled synthetic FASTA and the
  exact Longxiaodou 4 FASTA SHA recorded in
  `docs/evidence/issue45/old20_lineage_audit.json`. Other references require
  an input-bound validation review even for a small cohort.
- 21–51 samples: require `--joint_review_file`, a JSON review of public-scale
  validation bound to the exact code, reference bundle, FASTQ content and
  read-group metadata summary, and budgets. Equal input size is not equal input
  identity; swapping or changing FASTQ content invalidates an earlier review.
  Fill the generated `review_template` with reviewer/rationale after assessing
  the plan. This is not an accuracy or customer-release approval.
- Input FASTQs above 207,358,941,216 bytes also require that review, even if
  storage is available. This is the observed 51-sample input size in
  `docs/evidence/issue45/host_and_generation_observed.json`, not a guarantee
  that smaller inputs fit a resource allocation.
- More than 51 samples: refused, including 327, even if a review file is supplied.
- References above 500 Mb and non-diploid real references above 1 Mb are outside
  this initial envelope. These are conservative scope boundaries, not biological
  thresholds or claims of performance on every smaller dataset.

Defaults are `joint_interval_size_bp=20000000`,
`joint_genotype_memory_gib=16`, `joint_max_concurrency=3`, and
`joint_launch_memory_gib=110`. The last value also bounds the local executor's
scheduled memory. Set it for the actual host; it does not reserve RAM against
unrelated applications. The 20 Mb value is a maximum window span, not a safety
guarantee across cohorts. Higher variant density can require smaller windows.
Plans above 10,000 intervals are rejected before allocating the interval list;
this is an operational planning limit, not a measured capacity endorsement.

Storage uses observed input bytes and a clearly labeled **6x planning
assumption**, compared with free space on both the task and output filesystem.
It is not a proven upper bound and does not replace monitoring during the run.
The profile currently targets local filesystems and a local executor. Remote
object-store capacity and multi-host scheduling require separate support.

The existing per-read-group `HASH_INPUT_FASTQS` stage now precedes approval;
the planner reuses its measured checksums instead of hashing reads itself.
This stage and the reference-bundle validator use Nextflow `cache 'deep'` so a
same-size input change with unchanged mtime cannot reuse stale provenance.
Deep caching reads content to form cache keys and the task may then read it
again to validate/hash it. Budget this additional sequential I/O, including for
large reference indexes, before approval and on resume. Keep inputs immutable
while a run is active; this is not an atomic filesystem snapshot.

Successful runs publish `provenance/cohort.joint_plan.json` and
`provenance/cohort.joint_inputs.json`; the run manifest binds both with SHA-256.
Input tuples preserve sample/gVCF/index associations and headers are checked.
The sample map declares lexicographic sample order explicitly; downstream GS
metadata continues to follow the actual VCF header, never an assumed sheet
position. The data matrix is not relabeled by column position. The verified
input-summary and reference-bundle fingerprints also travel into the scientific
task metadata: a changed binding invalidates cached preprocessing/mapping even
when file size and mtime were preserved. A changed cohort input currently
invalidates preprocessing for the whole cohort, trading reuse for safety.

## Resume and retention

Retain Nextflow's launch-directory cache and `work/` until acceptance and any
agreed rerun window expire. A final report's retention differs from reusable
intermediates. No automatic cleanup is added by this change.

After Nextflow has stopped, before files can be changed, inventory successful
tasks' regular output files (input symlinks and task-root hidden command files
excluded; nested hidden database files are included):

```
python bin/resume_cache_inventory.py snapshot --work-dir /run/work --inventory /run/cache.json
python bin/resume_cache_inventory.py verify --work-dir /run/work --inventory /run/cache.json
```

Only then invoke Nextflow with the same launch directory and `-resume`. An
inventory detects same-size corruption that ordinary output-existence cache
checks can miss. Changed files and unexpected added directory members are
refused; use a fresh work directory or
explicitly remove/rebuild the affected task after reviewing the report. For
intentionally removed top-level file outputs, `verify --permit-missing` reports
the missing files and permits a resume that must re-execute those tasks.
Missing members inside directory outputs (such as GenomicsDB fragments) always
block: the parent directory can still satisfy Nextflow's cache existence check.
Use a fresh work directory or explicitly rebuild the whole affected task after
reviewing those failures. The inventory
does not delete anything and must not be regenerated over corruption to hide
it. Hashing large BAM/GenomicsDB caches has real I/O cost. Keep inventories
private with the run; they are operational artifacts, not public fixtures.

The bounded synthetic regression performs a controlled failure at GS matrix
generation, verifies stopped cache, resumes and checks biological equivalence,
then detects a same-length corrupted cache file, removes that one disposable
test output and confirms missing-output re-execution. It also changes the gzip
header of an isolated input copy while preserving size and filesystem mtime,
then requires a new input binding, non-cached preprocessing, and unchanged
biological output. No repository or real input is mutated. It refuses a non-empty
output directory:

```
python tests/scripts/check_issue69_resume.py --output-dir /tmp/issue69-new-run
```

Each invocation keeps its own numbered trace and log. Nextflow can complete
while refusing to replace an existing trace, so a later resume must never
reuse an earlier trace path for its cache assertions. CI runs this regression
in a parallel job and retains its synthetic traces, logs and completed evidence
as the `issue69-synthetic-resume` artifact, including diagnostics after failure.

This is wired into CI alongside a real GATK boundary-spanning indel fixture
and a nonlexical sample-ID / split-interval E2E case. Availability of these
tests is not evidence that they passed in an environment without Nextflow
and Docker; use the PR's recorded checks to assess validation status.

## Remaining real-data acceptance

Use the existing, authorized public 51-sample gVCFs and fixed reference rather
than downloading/re-calling the cohort. Freeze input checksums and the new DAG
SHA, first replay one window with the measured budget, then compare the
adopted interval set against #45's baseline: sample order, CHROM/POS/REF/ALT,
GT/AC/AN, boundary accounting and newly appearing QD differences. Record wall
time, CPU time, RSS, concurrency and retained bytes. Classify differences
using the existing `benchmarks/issue45` helpers; do not assume all QD changes
are acceptable because some were explained previously. Keep this acceptance
open until evidence exists. No 327-sample full run is authorized by this PR.
