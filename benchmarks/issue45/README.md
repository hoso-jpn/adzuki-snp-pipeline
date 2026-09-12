# Issue45 isolated execution helpers

These helpers operate outside production and benchmark checkouts. They do not
modify the frozen production modules or pipeline defaults. Public evidence belongs
in `docs/evidence/issue45`; FASTQ, BAM, gVCF, workspace and private operational logs
remain outside Git.

`stage_generation.py` consumes the fixed cohort selection and checksum-verified
download inventory. It rechecks every FASTQ SHA256 and reference checksum before
creating a new run directory. All 51 BioSamples must be unique. Its manifest
records the production SHA, wrapper/helper hashes, reference identity, input
checksums, library identity and scientific parameters before execution. The
benchmark-ready and lineage-verified flags remain false until output validation.
The production `bin` files are copied from the git-tracked file list, preserving
executable modes. Their SHA256 values are recorded and checked before launch and
after generation. An external `bin` symlink is insufficient: a task container may
mount the launch directory without mounting the symlink target. The regression
test mounts only the launch directory and invokes the real reference validator.
The old40 FASTQs are copied locally from their audited archive into the dedicated
input directory, with SHA256 verification before replacing the new run's own
symlinks. The archived originals are preserved. This adds 56.7 GB to the initial
1.9 TB conservative peak-storage estimate and avoids reliance on trash retention.
The first failed generation's 46,032,388 KiB (about 47.1 GB) are also retained;
the conservative total estimate therefore increases to about 2.004 TB. The new
launch still requires the original 2.95 TB free-space gate.

The generation template imports the frozen main modules for FASTP, reference
validation/indexing, mapping/sorting, per-sample merge, duplicate marking, BAM
index/QC and HaplotypeCaller. It deliberately ends at gVCF. The same modules and
production resource labels render the scientific commands. Each input run is one
read group, with the unchanged run accession as sample ID, ENA library name as
library ID, and no inferred physical platform unit. A first FASTQ header is not
evidence that every read came from one flowcell/lane. The old20 ENA library IDs
agree with the prior samplesheet.

Only publication policy changes in the isolated run: FASTP reports are copied,
trimmed reads and BAMs remain in the preserved work directory, and gVCFs are hard
linked into the output directory. This avoids large redundant copies. No FASTQ,
BAM, gVCF or work directory is removed.

Example (all paths below are caller-supplied local paths):

```bash
python3 benchmarks/issue45/stage_generation.py \
  --run-root "$RUN_ROOT" \
  --production-checkout "$PRODUCTION_CHECKOUT" \
  --cohort-selection docs/evidence/issue45/cohort_selection.json \
  --run-dir "$GENERATION_RUN"
```

`run_generation.py` is specific to the inventoried local trial container. Use it
only with authorization to temporarily stop that workload. It verifies the exact
container identity, localhost-only exposure, absence of active connections and
the model endpoint before stopping. The 51-sample generation receives 24 logical
CPUs and a 96 GiB aggregate Nextflow allocation. It requires 110 GiB available RAM
and 2.95 TB free storage at launch. Host metrics are sampled every five seconds.
Persistent available RAM below 8 GiB or swap growth over 512 MiB, or free storage
below 1.2 TB, aborts the owned run. These are conservative operational guards, not
measured 327-sample resource recommendations.

The controller restores the same trial container after success, Nextflow failure,
headroom failure or handled interruption, and verifies configuration fingerprints
and the model endpoint. It stops only orphan containers labeled for this exact
generation before restoring the large model. Operational state records exclude
environment values but still remain private. A failure to stop owned tasks or
restore the trial is an explicit error requiring follow-up.

```bash
python3 benchmarks/issue45/run_generation.py \
  --run-dir "$GENERATION_RUN" \
  --expected-trial-id "$INVENTORIED_TRIAL_ID"
```

`measure_process.py` runs **inside** a measurement container. Linux
`RUSAGE_CHILDREN` reports the tool and waited-for descendants, including the JVM
behind the GATK Python launcher; measuring the host Docker client would not report
that JVM's RSS. It preserves the subprocess exit status and refuses to overwrite
an existing measurement. Host available memory/swap and workspace size are
separate measurements.

The generation wrapper has been exercised with the repository's synthetic
multi-read-group fixture and compared to the full frozen-main workflow. Both
samples' complete gVCF record streams were byte-identical. This verifies the
wrapper path only; the real 51-sample lineage and E0–E4 gates are separate.

`validate_generated_cohort.py` checks all 51 gVCFs after generation and service
restoration finish. It binds the published hard links to successful task outputs,
checks the actual launcher/container and expanded HaplotypeCaller parameters,
matches ordered contigs and shared INFO/FORMAT/FILTER/ALT definitions, reads every
BGZF member to verify CRC/EOF, and compares the complete sequential and indexed
bcftools record streams. Reference FASTA, FAI, dictionary, launch files and each
gVCF/index are checksum verified. Only its successful final manifest opens the
lineage gate; a partially validated cohort leaves no successful final manifest.
The bcftools traversal containers use `--log-driver none`: their complete record
stream goes directly to the validator, avoiding a second sequencing-data copy in
Docker's root-filesystem logs. stderr remains in the validation directory. This
is a per-container option and does not change Docker daemon configuration.

```bash
python3 benchmarks/issue45/validate_generated_cohort.py \
  --run-dir "$GENERATION_RUN" \
  --production-checkout "$PRODUCTION_CHECKOUT" \
  --output-dir "$NEW_VALIDATION_DIRECTORY"
python3 benchmarks/issue45/run_benchmarks.py \
  --validated-cohort "$NEW_VALIDATION_DIRECTORY/validated_cohort.private.json" \
  --output-dir "$BENCHMARK_DIRECTORY" \
  --no-trial-pause
```

The benchmark controller requires 110 GiB available RAM and 2.0 TB free storage
before execution. Runtime host guards match generation. It pauses and restores the
authorized trial only when `--expected-trial-id` is given; the pause exists purely
to free RAM, so `--no-trial-pause` measures without touching any user workload once
the host already has that headroom. That mode still enforces the same launch gate
and runtime host monitor, and it records that nothing was stopped.

Three interval tasks run concurrently with 8 CPU each. GenomicsDBImport gets
16 GiB with a 13,107 MiB heap; its measured 51-sample peak was 1.55 GiB, so that
allocation is already ample and stays unchanged. GenotypeGVCFs has its own,
larger tier because the first real E0 was OOM-killed there (see below). Gather
uses 4 CPU / 8 GiB, and Reblock 4 CPU / 16 GiB. Docker memory+swap limits equal
the memory allocation for these targeted tasks. This fixed comparison policy
differs from unconstrained host scheduling and is recorded explicitly. No cache
dropping, host tuning, or production resource change is performed.

Every experiment shares one GenotypeGVCFs ceiling, so an allocation difference
is never a confound between a compared pair. The published evidence records that
uniformity explicitly and names the values that could otherwise confound it. The
16 GiB production-baseline E0 failure below is a distinct measurement under a
different allocation: it is retained as negative evidence about the current
production resource contract, and the common-ceiling comparison run does not
replace it. That experiment was not retried upward until it passed.

The first real E0 attempt OOM-killed GenotypeGVCFs on the longest contig
(65.4 Mb) at 51 samples: exit 247, `OOMKilled=true`, peak RSS 15.98 GiB against
its own 16 GiB ceiling, after the progress meter collapsed from about 1.2M to 3
records per minute. A collapse of that shape is a GC death spiral, so the live
set genuinely approached the 15 GiB heap rather than the container merely
clipping native allocations. This repository's standing methodology
(`nextflow.config`, Issues #30/#33) treats a reading at its own cgroup ceiling as
untrue, so GenotypeGVCFs was re-measured alone at a deliberately generous
ceiling, and the resulting tier is applied identically to every experiment.
Raising one allocation is a resource change only: no scientific parameter,
threshold, ploidy, reference or interval semantic differs, and the per-contig
E0/E1 plans are still the same plans. The completed E0 tasks also show
GenotypeGVCFs peak RSS growing with interval length at fixed sample count, which
is measured input to the interval-strategy decision rather than an assumption
carried into it. That relationship is an observation at 51 samples under this
reference, GATK version, ploidy and HaplotypeCaller condition, and the two
longest contigs' own readings were capped by the ceiling being tested; it is not
evidence of how memory scales with sample count, and nothing here extrapolates
it linearly to 327 samples.

The campaign is tens of hours long, so the controller is launched detached:
`setsid` puts it in its own session and process group, so closing the operator's
SSH connection delivers no SIGHUP to it, and `nohup` with a `/dev/null` stdin
means it never blocks on or dies with a terminal. Its own SIGTERM/SIGINT
handlers still stop the suite deliberately and stop the containers it owns. A
recorded controller PID makes a second concurrent launch into the same run
directory refuse rather than interleave.

A suite that stops part-way keeps its finished experiments. Re-running the same
command against the same output directory reuses every experiment whose result is
already `COMPLETED_AWAITING_COMPARATIVE_REVIEW`, and refuses to start when an
experiment directory holds retained failure evidence or no result at all, so a
new run directory is required to re-attempt one. A resumed run re-verifies the
validated cohort checksum, every helper checksum, the resource policy and the
interval plan, and stops if any of them changed since the directory was
started.
Task wall-time sums and complete experiment elapsed time are distinct metrics;
the latter also includes index/integrity verification and scheduling overhead.

The executor uses the preparation helper's plans, then independently checks exact
tiling and dictionary order. Each import must log both the 50-sample first batch
and one-sample second batch, plus both completions. Sample-name-map uses explicit
index paths and enables GATK's map validation in addition to the stricter external
header checks. Multiple intervals trigger GATK's documented serial reader
initialization fallback even when eight reader threads were requested; this is
recorded as part of the grouping result, without changing the requested flags.

All GenotypeGVCFs calls use `--only-output-calls-starting-in-intervals true` to give
each window ownership of variant starts. For whole-contig tasks this does not
exclude any reference coordinate. Synthetic whole-contig/split tests check the
resulting sample, variant, genotype and accounting contracts. The flag exists in
the pinned 4.6.2.0 executable, although its help marks it deprecated.

E3 invokes the actual tool name **ReblockGVCF** (singular), with explicit GQ bands
20/100, `keep-all-alts=true`, `floor-blocks=false`, and `drop-low-quals=false`.
Keeping alternate alleles avoids introducing optional allele dropping into this
first compression experiment. Reblocking still changes reference-confidence
representation and can change genotyping annotations/qualities; only measured
comparisons can justify adoption. Its outputs are a separate checksum-tracked
input set and never replace the production gVCFs.

The comparison preserves full-record hashes and separately audits every variant
key, INFO field and per-sample FORMAT change. The pinned
[QualByDepth implementation](https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/QualByDepth.java)
replaces raw QD at least 35 with a random draw around 30 (standard deviation 3).
Different task partitions can therefore change QD while QUAL, AD, GT and all
other fields remain identical. Such differences are classified only if the raw
QD reconstructed from the unchanged genotypes/AD exceeds 35 with a rounding
margin, every other record field agrees, and membership of the frozen main's
QD<2 filter is unchanged. All other differences remain unexplained until review;
neither QD nor any other annotation is silently discarded from evidence.

`summarize_evidence.py` turns a finished private run into the sanitized JSON
committed under `docs/evidence/issue45`. It selects, renames and aggregates what
the run already recorded and re-derives no measurement, so a disagreement with
the private record is a bug in the summarizer. Host directory layout, absolute
paths and usernames are never published; sample identities are public accessions
and stay as they are. An experiment that failed is reported with its status and
error type and without invented results.

```bash
python3 benchmarks/issue45/summarize_evidence.py \
  --run-dir "$BENCHMARK_DIRECTORY" \
  --validated-cohort "$NEW_VALIDATION_DIRECTORY/validated_cohort.private.json" \
  --output docs/evidence/issue45/benchmark_results.json
```

Each task retains command, log, allocation, exit/OOM state and measurements.
Only finished, explicitly named benchmark containers are removed; bind-mounted
data and logs remain. Tool failures retain their experiment evidence. E3 failure
does not prevent the independent E4 comparison; failures in their common baseline
stop the suite. Results always await comparative scientific review and do not
automatically select thresholds, architecture decisions or the 327-sample gate.
