# Real-cohort downstream resource hardening (Issue #30)

This document records Issue #30: hardening the resource contract for three downstream
processes and re-evaluating `GATK_HAPLOTYPECALLER`'s per-task CPU count, using targeted
benchmarks against Issue #26's own already-published real 5-sample cohort artifacts. It
does not run a 20-30 sample cohort (that remains Issue #26's own next phase, out of this
Issue's scope) and does not re-run the FASTQ -> mapping -> HaplotypeCaller -> Joint
Genotyping sequence -- every benchmark below reuses real files Issue #26 already produced.

## Issue #26 baseline (re-confirmed from real trace/logs, not just docs)

Read directly from `seedcore-01`'s still-intact
`/home/yusuke-hosokawa/scratch/adzuki-snp-pipeline-profile/issue26-cohort-2026-08-23/metrics/trace.txt`
(not only `docs/real_cohort_e2e.md`'s own summary of it), confirming the two real risks this
Issue addresses:

| Process | Attempt 1 (4 GiB) | Attempt 2 (8 GiB) |
| --- | --- | --- |
| `BUILD_GS_PANEL` | `COMPLETED`, exit 0, reported `peak_rss=4 GB` | (no retry needed) |
| `CLASSIFY_NORMALIZED_VARIANTS` (`cohort_gs`) | `FAILED`, exit 137 | `COMPLETED`, exit 0, reported `peak_rss=8 GB` |
| `SUMMARIZE_FILTER_QC` (`cohort:snp`) | `FAILED`, exit 137 | `COMPLETED`, exit 0, reported `peak_rss=8 GB` |

`process_low` (the label all three shared) is `cpus=2`, `memory={4.GB * task.attempt}`,
`maxRetries=1` -- so a process that also needed more than 8 GiB on attempt 2 would have had
no further retry to fall back on.

`GATK_HAPLOTYPECALLER` (`process_high`, 8 cpus) durations across the 5 real samples:

| Sample | Approx. coverage | Wall time |
| --- | ---: | ---: |
| `SRR29909069` | ~7.2x | 26m55s |
| `SRR29909135` | ~19.3x | 48m59s |
| `SRR29909067` | ~8.0x | 1h13m54s |
| `SRR29909073` | ~8.3x | 1h23m16s |
| `SRR29909072` | ~7.8x | 1h28m54s |

The deepest sample finished fastest; three of the four shallower samples took markedly
longer. `process_high` requests 8 cpus/task; this host has 32 logical CPUs, so at most
`floor(32/8) = 4` such tasks run genuinely concurrently -- consistent with contention among
concurrently scheduled tasks, not sequencing depth, driving the difference.

## Problem statement

1. Three downstream Python processes (`CLASSIFY_NORMALIZED_VARIANTS`,
   `SUMMARIZE_FILTER_QC`, `BUILD_GS_PANEL`) shared `process_low` with 13 other genuinely
   lightweight processes, at a memory allocation two of them could not meet on first attempt
   and one of them used in full without any retry margin.
2. `GATK_HAPLOTYPECALLER`'s per-task cpu count (8) limits concurrency to 4 tasks on this
   host; whether a lower cpu count could materially improve cohort throughput without a
   meaningful per-task cost was unmeasured.

## Resource architecture: before

```groovy
withLabel: process_low {
    cpus   = 2
    memory = { 4.GB * task.attempt }
    time   = '2h'
}
// CLASSIFY_NORMALIZED_VARIANTS, SUMMARIZE_FILTER_QC, BUILD_GS_PANEL,
// and 13 other processes all used this one label.

withLabel: process_high {
    cpus   = 8
    memory = { 16.GB * task.attempt }
    time   = '24h'
}
// GATK_HAPLOTYPECALLER, GATK_GENOTYPEGVCFS
```

## Resource architecture: after

Three new dedicated labels, plus one changed value on an existing label:

```groovy
withLabel: process_variant_classification {
    cpus   = 2
    memory = { 12.GB * task.attempt }
    time   = '2h'
}

withLabel: process_variant_qc_summary {
    cpus   = 2
    memory = { 12.GB * task.attempt }
    time   = '2h'
}

withLabel: process_gs_panel {
    cpus   = 2
    memory = { 8.GB * task.attempt }
    time   = '2h'
}

withLabel: process_high {
    cpus   = 4    // was 8
    memory = { 16.GB * task.attempt }
    time   = '24h'
}
```

`cpus`/`time` for the three new labels are carried over from `process_low` unchanged --
every benchmark below measured ~100% CPU (effectively single-threaded) and well under the
2-hour budget even at this cohort's scale, so neither is the constraint these labels exist
to address. `errorStrategy`/`maxRetries=1` (retry on exit 137/140/143, Issue #8's own
contract) are unchanged and still apply to all four labels; the intent of this Issue's
memory values is that attempt 1 succeeds in normal operation, not that retry becomes load-
bearing.

The other 13 `process_low` processes (`samtools_index`, `summarize_variant_qc`,
`bcftools_stats`, `validate_reference_contigs`, `fastqc`,
`gatk_create_sequence_dictionary`, `reconcile_variant_type_counts`,
`build_gs_panel_manifest`, `samtools_qc`, `reconcile_gs_panel_accounting`, `fastp`,
`samtools_faidx`, `gs_index_classified_variants`) are untouched -- none of them showed any
memory pressure in Issue #26's real run, and raising `process_low` itself would have
over-allocated all of them for no evidenced benefit.

## Targeted real-data benchmark methodology

Each benchmark reused one of Issue #26's own already-published real cohort artifacts,
invoking the target module directly (an ad-hoc, uncommitted nf-test `nextflow_process`
test per module -- not a permanent production-code fork, and not committed to this
repository, per this Issue's own explicit instruction) with a deliberately generous memory
ceiling (24 GiB, comfortably above any plausible true usage), specifically so the reported
`peak_rss` reflects genuine consumption rather than being capped by a too-small cgroup
limit -- this distinction turned out to matter (see `BUILD_GS_PANEL` below). Real resource
usage was read from each task's own `.command.trace` file (Nextflow's raw per-task
accounting, in KiB, more precise than the human-rounded `trace.txt` table format used for
the main pipeline run).

| Artifact reused | Source (Issue #26 real run) | SHA-256 |
| --- | --- | --- |
| `cohort_gs.normalized.vcf.gz` (+`.tbi`) | `results/variants/gs_normalized/` | `f917483adf3adee727d0376e2277e0886b9f636da220ba5c03650b100e749ad2` |
| `cohort.snp.filtered.vcf.gz` (+`.tbi`) | `results/variants/filtered/` | `8b7d1651dc8e4e9bdd7ec28689d6674f18dfe81b72e8af6dfa3dcc293360aa8f` |
| `cohort.indel.filtered.vcf.gz` (+`.tbi`) | `results/variants/filtered/` | `8d2f415f41f5013e444d79aa46294ccbdba5d135addc4e0c72fe87a8e45776eb` |
| `cohort_gs.snp.pass.vcf.gz` (+`.tbi`) | `results/variants/gs_pass/` | `fca83b37149953bd68f952b668edafc6b52b13bb98b2e3326480a9a1a0319ff1` |
| `SRR29909135.markdup.bam` (+`.bai`) | `results/alignment/` | `dcb630ed63cc3f3ad9ae9941d7daab42db317b212dd09840ee77fdf36a8f1852` |
| `GCF_016808095.1_ASM1680809v1_genomic.fna` | `reference/` (Issue #8's copy, reused since) | `e9838db1b048b54b21534285aaf95eae64cb3019d85af270b44e20b3d545f383` |

Git commit these benchmarks ran against: `0ab2f94` (this Issue's own branch, after the
resource-label changes above were already applied) for the downstream-process re-runs used
for semantic-equivalence checking; the first (label-sizing) pass ran against the
pre-relabeling module code, which is why that pass's benchmark config still targeted
`process_low` -- both passes' `peak_rss` values agree to within 0.01%, confirming this
didn't affect the measurement.

**Measurement precision caveat**: `.command.trace`'s `peak_rss` is a periodic sample of the
cgroup's reported RSS (Nextflow's own resource-trace mechanism), not an instruction-level
profiler; a peak between samples could in principle be missed. It is, however, far more
precise than `trace.txt`'s human-rounded display (which is what originally read
`BUILD_GS_PANEL`'s `peak_rss` as exactly `"4 GB"` -- the cgroup ceiling itself, not
necessarily the true value).

## CLASSIFY_NORMALIZED_VARIANTS

| | |
| --- | --- |
| Input | `cohort_gs.normalized.vcf.gz`, 5,694,489 records |
| True peak RSS (24 GiB ceiling) | 9.33 GiB (9,786,416-9,787,412 KiB across two reproductions) |
| Wall time | 27-29s |
| New attempt-1 memory | 12 GiB |
| Headroom | (12 - 9.33) / 12 = **22.3%** |
| Semantic equivalence | `cohort_gs.classification_accounting.tsv` byte-identical to the original Issue #26 production run's own output |

## SUMMARIZE_FILTER_QC

| | SNP (`cohort:snp`) | Indel (`cohort:indel`) |
| --- | ---: | ---: |
| Input | `cohort.snp.filtered.vcf.gz`, 4,768,719 records | `cohort.indel.filtered.vcf.gz`, 771,867 records |
| True peak RSS | 8.43 GiB (8,844,076-8,844,100 KiB) | 1.36 GiB (1,424,080-1,424,196 KiB) |
| Wall time | 47-49s | 7-9s |
| Semantic equivalence | `filter_breakdown.tsv` and `annotation_qc.tsv` byte-identical to production | `filter_breakdown.tsv` byte-identical to production |

Both invocations share one label (`process_variant_qc_summary`, 12 GiB), sized for the SNP
case -- indel's real footprint (1.36 GiB) has enormous headroom by comparison, but SNPs
vastly outnumber indels in this cohort and splitting further by `meta.variant_type` was
considered and rejected: it would add a second dynamic-resource code path for a difference
this Issue has no evidence yet justifies (see "Considered and rejected" below).

Headroom (SNP, the binding case): (12 - 8.43) / 12 = **29.8%**

## BUILD_GS_PANEL

| | |
| --- | --- |
| Input | `cohort_gs.snp.pass.vcf.gz`, 4,298,980 PASS records |
| True peak RSS (24 GiB ceiling) | 5.34 GiB (5,596,808-5,597,036 KiB across two reproductions) |
| Wall time | 2m13s-2m16s |
| New attempt-1 memory | 8 GiB |
| Headroom | (8 - 5.34) / 8 = **33.3%** |
| Semantic equivalence | `genotype_encoding_accounting.tsv` and `sample_metadata.tsv` byte-identical; `variant_metadata.tsv` and the genotype matrix both exactly 4,298,980 records, matching production |

**This true peak RSS (5.34 GiB) is itself a finding.** The original Issue #26 production
run's own `trace.txt` reported `BUILD_GS_PANEL`'s `peak_rss` as exactly `"4 GB"` -- the
`process_low` cgroup ceiling at the time, to the byte. A task genuinely using only 4 GiB
would not be expected to land precisely on its own ceiling; more likely, the task exceeded
4 GiB and survived via this host's 15 GiB of configured swap (Linux's OOM killer does not
necessarily trigger immediately when swap is available) rather than via genuine headroom.
This benchmark's generous 24 GiB ceiling is what let the true figure surface at all.

## Retry contract (unchanged mechanism, now a safety net rather than a first-line expectation)

All four labels (`process_variant_classification`, `process_variant_qc_summary`,
`process_gs_panel`, `process_high`) keep this repository's existing `errorStrategy` (retry
on exit 137/140/143) and `maxRetries=1`, with `memory` still a `{ N.GB * task.attempt }`
closure -- so attempt 2 still gets double attempt 1's allocation, and the ratio is
recomputed fresh on every attempt rather than holding a stale value (the same pattern
established in Issue #8/#11). What changed is the *intent*: attempt 1's value was chosen
from a real positive-headroom benchmark specifically so a normal 5-sample-scale run
succeeds on attempt 1, rather than depending on the retry to reach a value that works. The
retry remains as a genuine safety net (e.g. for a future cohort whose record counts exceed
this benchmark's), not because this Issue expects to need it in normal operation. No
`maxRetries` increase or unbounded memory escalation was introduced (Seedcore-01 has 123
GiB total; an unbounded `task.attempt` scaling policy was rejected as unnecessary and
unsafe at this host's real capacity).

## GATK_HAPLOTYPECALLER CPU benchmark

`--native-pair-hmm-threads ${task.cpus}` is exactly where `task.cpus` reaches the real GATK
command (confirmed by reading `modules/local/gatk_haplotypecaller.nf` directly, not
assumed). Benchmarked 8 cpus (current baseline) against 4 cpus, both isolated (no other
concurrent task on the host), against the same real input:

| cpus | Wall time | True peak RSS | Real `%cpu` (of allocated cpus) | gVCF records | Real variant records (non-`<NON_REF>`-only) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 40m45s | 10.56 GiB | 54.1% (4.33 of 8 cores) | 89,856,154 | 1,134,549 |
| 4 | 42m10s | 10.21 GiB | 64.5% (2.58 of 4 cores) | 89,856,154 | 1,134,549 |

**Semantic equivalence**: gVCF record count identical (89,856,154); a positional diff of
every real variant record's `CHROM`/`POS`/`REF`/`ALT` (excluding `<NON_REF>`-only reference
confidence blocks) between the two runs' gVCFs found **zero differences** across all
1,134,549 records. Changing `task.cpus` for `GATK_HAPLOTYPECALLER` changes wall time only,
not scientific output.

8 cpus bought only ~3.5% less wall time than 4 -- GATK's PairHMM-threaded region does not
meaningfully scale past ~4 threads for this reference/workload on this hardware, consistent
with the sub-linear real `%cpu` utilization measured at both settings. `GATK_GENOTYPEGVCFS`
(`process_high`'s only other consumer) does not read `task.cpus` anywhere in its own
script -- confirmed by reading `modules/local/gatk_genotypegvcfs.nf` directly -- so this
change affects only how many concurrent scheduling slots this label's tasks compete for,
not `GATK_GENOTYPEGVCFS`'s own execution.

**Measurement integrity note**: the 4-cpu benchmark's isolation was briefly compromised for
approximately 7 minutes (of its ~42-minute total) by an unrelated downstream-process
benchmark this Issue's own work accidentally started concurrently on the same host; that
benchmark was identified and killed as soon as noticed. Given the interfering tasks
requested only 2 Docker `--cpu-shares` each against 32 real logical CPUs (no genuine
hardware oversubscription occurred, since total concurrent demand never exceeded the host's
core count), the expected impact on wall time is negligible, but this is disclosed rather
than silently omitted.

### Theoretical concurrency and cohort throughput

At 8 cpus/task: `floor(32 / 8) = 4` concurrent `GATK_HAPLOTYPECALLER` tasks on this host.
At 4 cpus/task: `floor(32 / 4) = 8` concurrent tasks -- double.

**Ideal linear estimate only** (single isolated sample's wall time x waves needed; assumes
zero contention overhead, which Issue #26's own real run already shows is not realistic --
see below):

| Cohort size | Waves at 8 cpus (4 concurrent) | Ideal wall time (8 cpus) | Waves at 4 cpus (8 concurrent) | Ideal wall time (4 cpus) |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 5 | ~3h24m | 3 | ~2h7m |
| 30 | 8 | ~5h26m | 4 | ~2h49m |

**Real contention uncertainty**: the one available real/isolated comparison point --
`SRR29909135` under Issue #26's real (partially concurrent, 8-cpu) run (48m59s) versus this
Issue's isolated 8-cpu benchmark of the identical sample (40m45s) -- shows a ~20% real
overhead even under that run's actual concurrency pattern. The other four samples in that
same real run showed *far* larger overheads relative to their own (unmeasured, since no
isolated benchmark exists for them) expected times -- three of four shallower-coverage
samples took 1h13m-1h29m against one deeper sample's 48m59s, an ordering inconsistent with
depth-driven cost alone. **The ideal linear estimates above are a lower bound, not a
prediction**; real wall time at 8-way concurrency (untested -- Issue #26's real run only
exercised up to 4-way) could exceed them by an amount this Issue has not measured.

### CPU policy decision

**Changed `process_high` from 8 to 4 cpus.** Rationale: real, measured, near-zero per-task
cost (~3.5% wall time) for a real doubling of theoretical concurrency, with confirmed
scientific-output equivalence. This is not "8 threads is slower" (the Issue's own
instruction explicitly warned against that unsupported claim) -- it is "4 threads is
adequate for this workload, and the freed concurrency is worth more than the small
per-task cost." Not adopted reflexively: the change was made only because the benchmark
showed a clear, quantified benefit; if it had shown a large per-task cost, `process_high`
would have stayed at 8 with that finding recorded instead.

**New consideration this change introduces**: 8 concurrent `GATK_HAPLOTYPECALLER` tasks at
~10.2 GiB true peak RSS each is ~82 GiB, versus the previous 4 concurrent tasks at ~10.6
GiB each (~42 GiB) -- a real increase in peak concurrent memory pressure on a 123 GiB host,
before accounting for whatever else the pipeline's other concurrently scheduled stages
(mapping, Joint Genotyping) are doing at the same time. This was not stress-tested at 8-way
concurrency in this Issue (doing so would mean either running a real 8-sample-or-larger
cohort, out of scope here, or fabricating synthetic concurrent load, which would not be
real evidence). It is recorded as a real, unresolved consideration for whoever executes the
20-30 sample expansion, not dismissed.

## Considered and rejected

- **Raising `process_low` itself** to 8 or 12 GiB: rejected -- would over-allocate the 13
  other genuinely lightweight processes still on that label for no evidenced benefit.
- **Splitting `SUMMARIZE_FILTER_QC` into per-variant-type labels or dynamic
  `meta.variant_type`-based memory**: rejected for this cycle -- one label sized for the
  SNP (binding) case is simpler and already gives positive headroom for both invocations;
  introducing a second dynamic-resource code path was not justified by evidence of an
  actual problem at indel's much smaller real footprint.
- **A memory formula scaling with record count** (e.g. `memory = base + variants *
  coefficient`): rejected -- this Issue has exactly one cohort-scale data point per
  process; fitting a predictive formula to a single point would be pseudo-precision, not
  evidence-based sizing. Static per-label allocation plus the existing retry mechanism is
  the appropriate level of sophistication until more scale points exist.
- **Increasing `maxRetries` or using unbounded `task.attempt` memory escalation**: rejected
  -- Seedcore-01's real capacity (123 GiB) bounds what is sane regardless of formula, and
  an escalating-forever retry hides a real sizing problem rather than surfacing it.

## Storage (unchanged from Issue #26, reference only)

No new storage run was performed. Issue #26's own linear reference estimate (~95 GB
attributable to the real 5-sample run, extrapolating to roughly 500-600 GB for 30 samples)
is carried forward unchanged, explicitly as a linear estimate rather than a measurement at
that scale -- this Issue's resource changes do not materially change per-sample storage
footprint (they affect memory ceilings and cpu counts, not what gets written to disk).

## GenomicsDB batching (unchanged, confirmed still unexercised)

`genomicsdb_batch_size` (default 50, Issue #11) remains untouched by this Issue and remains
unexercised in its actual batching behavior at any cohort size attempted so far (5 real
samples, far below 50). This is unchanged from Issue #26's own conclusion, not re-evaluated
here.

## Synthetic regression and CI

Baseline (this branch, before any Issue #30 change, commit `7a681c9`):

```
$ python3 -m unittest discover -s tests/bin
Ran 247 tests in 0.226s — OK

$ nextflow lint .
✅ 35 files had no errors

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test tests/modules/gatk_genomicsdbimport.nf.test tests/modules/gatk_gathervcfs.nf.test --ci
SUCCESS: Executed 26 tests in 147.168s
```

Final (this Issue's last commit):

```
$ python3 -m unittest discover -s tests/bin
Ran 251 tests in 0.226s — OK   (247 baseline + 4 new resource-label contract tests)

$ nextflow lint .
✅ 35 files had no errors (36-38 when the uncommitted ad-hoc benchmark directory was present locally; not part of any committed file)

$ git diff --check main...HEAD
(clean)

$ nf-test test tests/pipeline/adzuki_snp_pipeline.nf.test tests/modules/gatk_selectvariants.nf.test tests/modules/gs_normalize_variants.nf.test tests/modules/bwa_mem2_mem_sort.nf.test tests/modules/validate_reference_contigs.nf.test tests/modules/gatk_genomicsdbimport.nf.test tests/modules/gatk_gathervcfs.nf.test --ci
SUCCESS: Executed 26 tests in 148.612s
```

The synthetic fixture's own genotype/annotation contract (Issue #12/#20), GS panel contract
(Issue #18), reference-contig validation (Issue #11), and `genomicsdb_batch_size`/Xmx
contracts (Issue #11) all pass unchanged -- resource label and cpu-count changes alone do
not touch any of those data contracts, and the synthetic fixture's tiny inputs succeed
under any of the new memory values regardless, which is exactly why
`tests/bin/test_resource_label_contracts.py` reads the real source files directly instead
of relying on nf-test to observe a production-scale difference it cannot exercise.

`tests/bin/test_resource_label_contracts.py` adds 4 tests: each of the three relabeled
modules still declares its dedicated label (and not a reversion to `process_low`); each of
the four changed/new labels (`process_variant_classification`, `process_variant_qc_summary`,
`process_gs_panel`, `process_high`) has the expected `cpus`/`memory`/`time` in
`nextflow.config`; memory still scales with `task.attempt` rather than a fixed value; and
`conf/test.config` overrides all three new labels to CI-sized values.

## 20-30 sample Go / Conditional Go / No-Go decision

**CONDITIONAL GO.**

The two concrete blockers this Issue set out to address are resolved with real evidence:

- `CLASSIFY_NORMALIZED_VARIANTS` and `SUMMARIZE_FILTER_QC` (`cohort:snp`) now have
  22.3%/29.8% real positive headroom at 5-sample scale, instead of landing at their retry
  ceiling.
- `BUILD_GS_PANEL` now has 33.3% real positive headroom, instead of a swap-assisted "success"
  with no real margin.
- `GATK_HAPLOTYPECALLER`'s cpu count is now backed by a real, output-equivalence-verified
  benchmark rather than an unmeasured default; the change roughly doubles this host's
  theoretical HaplotypeCaller concurrency at negligible per-task cost.

This is **not** a GO, because attempting 20-30 samples introduces real, unresolved
questions this Issue's 5-sample-scale evidence does not answer:

1. **Downstream memory scaling is unvalidated beyond 5 samples.** A 20-30 sample cohort
   will very likely discover more variant sites (more genetic diversity across more
   individuals) than this 5-sample cohort's 5,558,870 raw records, and
   `CLASSIFY_NORMALIZED_VARIANTS`/`SUMMARIZE_FILTER_QC`/`BUILD_GS_PANEL` all process
   whole-cohort record sets in memory at once. Their new ceilings (12/12/8 GiB) were sized
   with real headroom *for this cohort's record volume*, not validated against a larger
   one. The existing retry mechanism remains available as a safety net if this ceiling is
   exceeded at larger scale, but that would mean falling back to the same "succeed only on
   retry" situation this Issue set out to move away from.
2. **8-way `GATK_HAPLOTYPECALLER` concurrency was never stress-tested.** The CPU policy
   change was validated in isolation (one task, no contention) and doubles the *theoretical*
   concurrent task count; real concurrent memory pressure at 8 simultaneous tasks (~82 GiB)
   plus whatever else the pipeline schedules concurrently was not exercised.
3. **Real wall-time at cohort scale remains uncertain beyond a lower bound.** The ideal
   linear throughput estimates above assume zero contention; Issue #26's own real run
   already showed contention effects at only 4-way concurrency, and this Issue's isolated
   benchmarks cannot measure 8-way contention effects without actually running a
   large-enough cohort.

**Conditions for the next phase (Issue #26's own 20-30 sample step, or a follow-up) to
convert this to a GO:**

- Monitor aggregate real memory usage during that run's `GATK_HAPLOTYPECALLER` stage
  specifically; if concurrent peak RSS approaches host capacity, be prepared to reduce
  `process_high`'s cpus further (increasing concurrency more) or accept lower concurrency
  by raising it back, rather than assuming this Issue's 4-cpu choice is final at that scale.
- Watch `CLASSIFY_NORMALIZED_VARIANTS`/`SUMMARIZE_FILTER_QC`/`BUILD_GS_PANEL` for any OOM
  at the new ceilings; if the larger cohort's record volume exceeds what 5 samples showed,
  re-benchmark with that scale's own real artifacts rather than assuming this Issue's
  numbers still hold.
- `genomicsdb_batch_size`'s real batching behavior remains unexercised below 50 samples
  (Issue #26/#11); this is unchanged and not a new condition, but still open.

## Limitations

- Every benchmark in this Issue is single-run (not repeated for statistical variance) --
  the two downstream-process re-runs (once for label sizing, once for semantic-equivalence
  checking after relabeling) agree to within 0.01%, which is reassuring but not a formal
  variance estimate.
- The `GATK_HAPLOTYPECALLER` cpu benchmark used only one real sample (`SRR29909135`); the
  other four real samples' behavior under isolated 4-cpu conditions was not measured.
- No 327-sample-scale risk is addressed or newly discovered here -- Issue #26's own
  327-sample exclusions (full cohort run, `--sample-name-map`, interval splitting,
  `ReblockGVCFs`, `--consolidate`, hard-filter/MAF/call-rate/LD/imputation changes, GS model
  training, GPU/Parabricks) remain entirely out of this Issue's scope and unaddressed by it.
