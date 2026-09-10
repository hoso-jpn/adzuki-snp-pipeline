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
The old40 FASTQs are copied locally from their audited archive into the dedicated
input directory, with SHA256 verification before replacing the new run's own
symlinks. The archived originals are preserved. This adds 56.7 GB to the initial
1.9 TB conservative peak-storage estimate and avoids reliance on trash retention.

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
