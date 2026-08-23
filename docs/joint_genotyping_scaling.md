# Joint Genotyping scale hardening (Issue #11)

This document records Issue #11: hardening the existing Joint Genotyping architecture
(`GATK_HAPLOTYPECALLER` -> per-contig `GATK_GENOMICSDBIMPORT` -> `GATK_GENOTYPEGVCFS` ->
`GATK_GATHERVCFS`, first implemented around PR #10) ahead of real multi-sample cohort
ingestion, plus the real single-accession observation that motivated it. It is not a
rewrite of that architecture, and it does not validate cohort scale: the target is a
327-sample cohort (Issue #26), and everything measured here is either synthetic-fixture
evidence or a single real accession on a single machine.

## Real-data observation that motivated this Issue

Issue #8 Phase 5 left a real pipeline run going on `seedcore-01`
(`RUN_ROOT=/home/yusuke-hosokawa/scratch/adzuki-snp-pipeline-profile/phase5-2026-08-23`)
against the real Longxiaodou 4 reference (`GCF_016808095.1`) and a real WGS accession
(`SRR29909135`). `GATK_HAPLOTYPECALLER` completed successfully (40m46s, peak RSS 10.7 GB),
producing a real ~1.1 GB gVCF and its `.tbi` index. The run then **failed** at
`GATK_GENOMICSDBIMPORT`, before any task was even dispatched (no `trace.txt` row for it at
all), with:

```
the number of gVCFs and indexes must match: 1176577035 gVCFs, 642455 indexes
```

Those two numbers are not gVCF/index counts -- they are literally the real gVCF's and its
`.tbi` index's file sizes in bytes. The root cause: Nextflow's `path` input qualifier
silently unwraps a single-element `List` into a bare scalar `Path`/`File` object when
exactly one file resolves for that channel slot. This is documented Nextflow behavior, not
a bug specific to this pipeline. `GATK_GENOMICSDBIMPORT`'s (and, identically,
`GATK_GATHERVCFS`'s) Groovy `script:` block called `.size()` on `gvcfs`/`gvcf_indexes`
assuming they were always `List`s; on a scalar `File`, Groovy's `.size()` returns the
file's byte length instead. Every synthetic fixture in this repository uses two samples,
so this never triggered before a real single-accession run hit it.

**This is one accession's evidence, not a multi-sample or cohort-scale validation.** It
identified a real bug and confirmed the pre-existing Xmx-ratio problem below is real; it
says nothing about `GenomicsDBImport`'s behavior with tens or hundreds of samples
(Issue #26).

## Batch-size contract

`params.genomicsdb_batch_size` (default `50`, minimum `1`, validated as an integer via
`nextflow_schema.json`) is passed to `GATK_GENOMICSDBIMPORT`'s `--batch-size`, controlling
how many samples GenomicsDBImport opens concurrently while building each per-contig
workspace.

`50` is a **current initial operational value**, carried over from GATK's own general
tooling conventions -- it is **not** validated as optimal for this pipeline's target
327-sample cohort, or for any intermediate 20-30 sample scale (Issue #26). No measurement
in this repository informs what value would actually be best at those scales; this
parameter exists so that value can be tuned later without a code change, not because `50`
has been shown to be correct here.

Because this repository's synthetic fixtures use only 1-2 samples, the raw cohort VCF's
actual content (`REF`/`ALT`/`GT`/`AC`/`AN`/`AF`/`DP`) is identical regardless of
`--batch-size` -- output equivalence cannot prove this parameter is actually wired up.
`tests/pipeline/adzuki_snp_pipeline.nf.test` instead reads each `GATK_GENOMICSDBIMPORT`
task's own real, fully rendered `.command.sh` from disk and asserts `--batch-size`
literally appears with the expected value, for both the default (`50`) and a
non-default (`7`) run.

## Memory contract

`GATK_GENOMICSDBIMPORT`'s actual memory footprint is JVM heap (`-Xmx`) plus a separate
native/C++ TileDB storage layer (buffers, mmap regions, thread stacks) that the JVM heap
setting does not account for at all. The previous formula --
`memory_gb = task.memory.toGiga().intValue() - 1` -- reserved a fixed 1 GiB regardless of
`task.memory`'s actual size. At the existing 16 GiB `process_high` allocation this is
Xmx = 15 GiB = 93.75% of `task.memory`, leaving the native layer only 6.25% headroom.

Xmx is now capped at a flat 80% of `task.memory`, computed in MiB (matching the precedent
set for `BWA_MEM2_MEM_SORT` in Issue #8) rather than rounded to whole GiB, so a small
`task.memory` cannot round in a way that overshoots 80%:

```groovy
total_memory_mib = task.memory.toMega()
xmx_mib = Math.floor(total_memory_mib * 0.8).intValue()
```

The required >= 20% native-layer headroom falls out of this same 80% ceiling
automatically -- it is not a separately asserted number. If `xmx_mib` would fall below 1,
the task fails fast with a diagnosable message naming the `process_genomicsdb` label,
mirroring `BWA_MEM2_MEM_SORT`'s existing fail-fast style, rather than silently running
GenomicsDBImport with an unusable heap.

`GATK_GENOMICSDBIMPORT` moves to its own `process_genomicsdb` resource label (in
`nextflow.config` and `conf/test.config`), independent of `process_high`, so its memory can
be tuned separately once real multi-sample data (Issue #26) reveals its actual footprint.
Its `cpus`/`memory`/`time` are carried over from `process_high` unchanged (`cpus=8`,
`memory = 16.GB * task.attempt`, `time='24h'`) -- no new numbers are introduced without
real justification.

`task.attempt`-based OOM-retry scaling (established in Issue #8) is preserved: `memory` is
still `{ 16.GB * task.attempt }`, and because `xmx_mib` is computed fresh from
`task.memory` inside the `script:` block on every execution, the 80% ratio is recomputed
correctly on every retry attempt rather than holding a stale value from the first.

`--reader-threads ${task.cpus}` is unchanged -- there is no evidence from either the
synthetic fixture or the real single-accession run that it needs to be.

`GATK_GATHERVCFS` received the identical List-safety fix (found via the same real-data
evidence -- it has the exact same `path(vcfs)`/`path(vcf_indexes)` + `.size()` pattern),
but its Xmx formula and `process_medium` label are intentionally unchanged: Issue #11's
memory-contract scope is GenomicsDBImport specifically.

## Reference contract

Neither the pipeline-generated reference-bundle path (`SAMTOOLS_FAIDX` +
`GATK_CREATE_SEQUENCE_DICTIONARY`) nor the prebuilt path (`--reference_fai`/
`--reference_dict`) previously confirmed that the `.fai` and `.dict` actually describe the
same reference, in the same contig order, with the same declared lengths, before
`GATK_HAPLOTYPECALLER`/`GATK_GENOMICSDBIMPORT`/`GATK_GENOTYPEGVCFS` started consuming both.
A mismatched pair -- most concretely, a correct `.fai`/`.dict` whose contigs are the same
*set* but listed in a different order, which no downstream GATK tool is guaranteed to
reject up front -- would otherwise only surface as whatever error message, at whatever
point in a run, the first GATK process that happens to notice chooses to produce.

`bin/validate_reference_contigs.py` parses both files and compares them **positionally**
(name and length, in file order) -- not as sets -- and reports the first point of
disagreement:

```
validate_reference_contigs.py: error: reference FAI and sequence dictionary are
inconsistent: first mismatch at index 12: FAI: Chr03 length=... DICT: Chr04 length=...
```

`modules/local/validate_reference_contigs.nf` wraps this as a pass-through gate, wired
into `workflows/adzuki_snp_pipeline.nf` right after `reference_fai_ch`/`reference_dict_ch`
are established for *either* reference path -- both converge to the same channel shape at
that point -- so every downstream reference-dependent process now depends on this gate's
output rather than on the pre-validation channels directly. A mismatch fails the whole run
before any GATK process starts.

Test coverage (`tests/modules/validate_reference_contigs.nf.test`,
`tests/pipeline/adzuki_snp_pipeline.nf.test`) includes the specific case a naive
set-equality comparison would miss: the same contig set, same lengths, different order.

## Consolidation policy: `--consolidate` is not enabled

GenomicsDBImport's `--consolidate` flag merges a workspace's incremental TileDB fragments
into fewer, larger ones after import, at the cost of extra I/O and time during the import
itself. It is not enabled here, and no new `genomicsdb_consolidate` parameter is added by
this Issue.

Re-evaluation conditions (i.e., when this should be revisited, not implemented now):

- A real multi-sample cohort run (Issue #26, 20-30 samples or the full 327) shows many
  small per-batch imports into the same per-contig workspace, producing enough fragments
  that downstream `GenotypeGVCFs` reads measurably slow down.
- GenomicsDBImport is ever run incrementally against an *existing* workspace (adding
  samples to a cohort already imported) rather than building each per-contig workspace
  fresh in one run, which is the only mode this pipeline currently uses.

Neither condition has been observed; this is a documented non-use decision, not an
oversight.

## Real reference contig/scaffold distribution (context only, not acted on)

The real Longxiaodou 4 reference (`GCF_016808095.1_ASM1680809v1_genomic.fna.fai`, the same
file profiled in Issue #8 Phase 5) has **36 contigs total**, summing to 448,362,642 bp:

| Length bucket | Contigs | Notes |
| --- | ---: | --- |
| >= 10 Mb | 11 | The 11 chromosome-scale pseudomolecules (`NC_068970.1` .. `NC_068980.1`, 65.4 Mb down to 27.7 Mb) |
| 100 kb - 1 Mb | 2 | `NC_021092.1` (404,466 bp) and `NC_021091.1` (151,683 bp) -- organellar-scale sequences by size, not independently confirmed here as chloroplast/mitochondrial |
| 10 kb - 100 kb | 5 | Unplaced scaffolds (`NW_...`) |
| < 10 kb | 18 | Unplaced scaffolds, smallest 1,000 bp |

Because this workflow creates one `GenomicsDBImport`/`GenotypeGVCFs` task pair per FAI row
(`intervals_ch` in `workflows/adzuki_snp_pipeline.nf`), a real run against this reference
launches 36 such task pairs regardless of sample count. The 23 small unplaced scaffolds
(most under 10 kb) are the most obvious candidates for a future interval-splitting or
scaffold-grouping strategy, so that a 327-sample cohort does not pay 23 separate
per-interval task-scheduling overheads for a combined ~0.15% of the genome. **This Issue
does not implement that.** It is recorded here as input for whoever scopes that work
(tracked as a follow-up below and in Issue #26), not as a decision already made.

## Targeted real-data smoke test

Reused the exact real gVCF/index already produced by `GATK_HAPLOTYPECALLER` in the Issue #8
Phase 5 run above -- the FASTQ -> mapping -> HaplotypeCaller sequence (which took over 40
minutes for HaplotypeCaller alone) was **not** re-run. Provenance:

- Source: `seedcore-01`,
  `.../phase5-2026-08-23/work/f6/3e0f78230d6c69394ad78d714824af/SRR29909135.g.vcf.gz{,.tbi}`
- Copied to a dedicated, checksummed location:
  `.../issue11-genomicsdb-smoke-2026-08-23/gvcf/` and `.../reference/`, with
  `PROVENANCE.sha256` recording every file's checksum. The gVCF and index checksums match
  the originals exactly (`a7935d71...` and `4989d027...` respectively) -- byte-identical
  copies, not regenerated.

Using the new `GATK_GENOMICSDBIMPORT` code directly (an ad-hoc, uncommitted nf-test module
test -- not part of this repository's permanent suite, since it depends on a real-data path
that exists only on `seedcore-01`), this single real gVCF/index pair (exactly the scenario
that broke the original run: one sample resolving to a scalar, not a `List`) was imported
successfully for two real contigs:

| Contig | Length | Real command (`.command.sh`) | Result |
| --- | ---: | --- | --- |
| `NC_068970.1` (largest chromosome-scale contig) | 65,407,200 bp | `gatk --java-options "-Xmx13107m" GenomicsDBImport --variant SRR29909135.g.vcf.gz --genomicsdb-workspace-path interval_000001_NC_068970.1.genomicsdb --intervals 'NC_068970.1' --reader-threads 8 --batch-size 50 --tmp-dir .` | Success, 67.5s (`process_genomicsdb` overridden to cpus=8/memory=16 GB for this smoke test; 16384 MiB * 0.8 = 13107 MiB, matching `-Xmx13107m` exactly) |
| `NW_026294847.1` (small unplaced scaffold) | 2,353 bp | Same shape, `--intervals 'NW_026294847.1'` | Success, 3.2s |

This confirms, on real production-scale gVCF content (not a synthetic fixture): the
List-unwrap fix works against a genuine single-sample Nextflow channel (not just reasoned
about in Groovy), `--batch-size` and the 80%-Xmx-ceiling formula both reach the real
command line correctly, and GenomicsDBImport itself succeeds where it previously failed.

**This remains one accession's evidence.** It does not exercise `GenomicsDBImport` with
more than one sample, does not exercise `--batch-size` at a value where it would actually
change import behavior (batching only matters once sample count exceeds the batch size),
and says nothing about memory behavior once GenomicsDB workspaces hold tens or hundreds of
samples' data. `-Xmx` was overridden to a generous 16 GB for this smoke test specifically
so that a real chromosome-scale import would not be memory-constrained by
`conf/test.config`'s tiny CI-sized values; this is not evidence that 16 GB (or any other
value) is right for a 327-sample cohort.

`-Xmx` was NOT verified via Nextflow's `-resume` against the original Issue #8 Phase 5 run
directory. Introducing `VALIDATE_REFERENCE_CONTIGS` upstream of `GATK_HAPLOTYPECALLER`
changes what physically flows into HaplotypeCaller's `path(fasta)`/`path(fai)`/`path(dict)`
inputs (re-staged through a new task's work directory, even though the content is
byte-identical), which risks invalidating Nextflow's task-level resume cache for
HaplotypeCaller and forcing a ~40-minute re-run this Issue explicitly should not require.
The direct real-gVCF smoke test above was chosen instead, specifically to avoid that risk.

## Explicitly out of scope for this Issue

Documented here as follow-up items only -- not implemented, and no additional accessions
were downloaded to explore them:

- `--sample-name-map` (a GenomicsDBImport alternative to repeated `--variant` flags,
  relevant once sample counts make the current command-line-argument approach unwieldy)
- Interval splitting or scaffold-grouping (see the contig-distribution analysis above)
- `ReblockGVCFs`
- Issue #26 in full: validating this contract at 20-30 samples, then 327 samples

## Baseline and full regression

Re-established on `seedcore-01` (docker/nextflow/nf-test are not available on the
development machine used for this Issue) against `main` at commit `92743ec` (PR #27's
merge commit) before any Issue #11 change:

```
$ python3 -m unittest discover -s tests/bin
Ran 197 tests in 0.213s — OK

$ nextflow lint .
✅ 34 files had no errors

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test --ci
SUCCESS: Executed 14 tests in 82.289s
```

Final, on this Issue's branch:

```
$ python3 -m unittest discover -s tests/bin
Ran 220 tests in 0.223s — OK

$ nextflow lint .
✅ 35 files had no errors

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test --ci
SUCCESS: Executed 24 tests in 135.057s
```

All pre-existing tests -- including the `REF`/`ALT`/`GT`/`AC`/`AN`/`AF`/`DP` genotype
contract from Issue #12/#20 -- pass unchanged.
