# FILTER QC streaming benchmark (Issue #43)

## Decision and review

The production path retains counters for distinct FILTER values/tags and seven
annotations instead of every VCF record. FILTER, annotation, missing/NA and
multi-tag semantics are unchanged on valid inputs. Resource allocations and
scientific thresholds are unchanged.

The review additionally reproduced missing/duplicate/misplaced VCF headers and
partial final outputs after a late output-write failure. The CLI now validates
header/row structure, stages all three outputs, rejects output/input aliases and
rolls back ordinary publication failures while preserving earlier results.
This does not claim a multi-file transaction across process/host crashes.
The failing-before regression cases now pass, including a failed second rename.

## Input identity

Public BioProject **PRJNA1138464**, reference **GCF_016808095.1 / Longxiaodou 4**.
Only the existing Issue #33 artifacts were reused. No FASTQ download, E2E rerun,
reference mixing or sample renaming was performed.

| Samples | Source pipeline SHA | SNP records | SNP VCF SHA256 |
| --- | --- | ---: | --- |
| 10 | `f129a01ada403a138b01edcfaea4c261b25f8eb1` | 8,748,064 | `a9633638711ead06e5cb4792870243006db54d66a78fe31246165e03d3fbedb7` |
| 20 | `556f38fd93008e0e4093a2fbd836e8d141afd8a4` | 10,296,980 | `180aefa9a242cb189849582e15636979d2b1219e9a4b80bc1f5c5725f103ed3c` |

Before use, the run manifests were checked byte-for-byte against the committed
10/20-sample manifests. FASTA/FAI/dict SHA256, ordered VCF contigs, sample count and
order, raw VCF manifest checksums, filtered record counts and GATK 4.6.2.0 headers
were independently verified. The source Nextflow version is 26.04.6.
Full sample accessions and reference/input checksums are in
[`seedcore_replay.json`](evidence/issue43/seedcore_replay.json).

## Measurement design

seedcore-01 had only about 3 GiB available RAM at preflight. Its unrelated workload
was left intact. All four verified SNP/indel VCFs were replayed through the new
script there, under a 512 MiB address-space limit; the actual RSS was about 16 MiB.
The old implementation was not launched on this memory-constrained host.

The two SNP VCFs were copied read-only into a local temporary directory and their
SHA256 verified again. Old and new scripts were then measured sequentially on the
same local review host with Python 3.12.3 and `/usr/bin/time -v`. The host has
46 GiB RAM; each run required at least 32 GiB MemAvailable at launch and used a
28 GiB address-space ceiling. The measured peaks are below that ceiling.
Old code is `0aebc73efe298b82ace100152d10fbaf6d216283:bin/summarize_filter_qc.py`;
new code is the reviewed streaming implementation at `b221ba8` (full SHA and script
checksums are recorded in the JSON). This measures Python RSS, not Docker CLI RSS.

Local nf-test runs were also active. Wall times are observations under that load,
not a controlled throughput claim; do not compare local and seedcore wall times
as a speed ratio. Host swap-used and swap-in/out counters are recorded before and
after each run; these are host counters, not exclusively attributable to the CLI.

## Paired local measurements

| Samples | Implementation | Peak RSS (MiB) | Wall (s) | Swap-used delta (bytes) |
| --- | --- | ---: | ---: | ---: |
| 10 | old | 16634.99 | 321.53 | 4096 |
| 10 | new | 15.70 | 203.17 | 0 |
| 20 | old | 20349.39 | 387.31 | 0 |
| 20 | new | 16.20 | 300.96 | 0 |

Both sample scales produced **byte-identical FILTER breakdown, annotation QC and
summary files** between old/new and against the corresponding historical
production outputs. No semantic normalization was needed. New RSS remains
about 16 MiB as SNP record count rises from 8.75 million to 10.30 million; the old
implementation grows from 16.25 to 19.87 GiB. The bounded-counter code and the
4,000→64,000-record regression test support the same conclusion. This is not a
claim that arbitrary malformed FILTER alphabets or record widths use constant RAM.

The old 10-sample run recorded 4 KiB host swap-used growth (and 4 KiB swap-out).
At completion about 36 GiB was available; the measured process was below its
28 GiB ceiling. We report this small host delta rather than calling it zero or
claiming that it proves swap-free isolation. All other paired swap-used deltas
were zero. Exact start/end counters and output checksums are in
[`local_comparison.json`](evidence/issue43/local_comparison.json).

## seedcore-01 confirmation

| Samples | Type | Records | New RSS (MiB) | Wall (s) | Swap-used delta (bytes) |
| --- | --- | ---: | ---: | ---: | ---: |
| 10 | snp | 8,748,064 | 16.00 | 34.28 | 0 |
| 10 | indel | 1,422,224 | 16.25 | 5.98 | 0 |
| 20 | snp | 10,296,980 | 16.00 | 49.47 | 0 |
| 20 | indel | 1,649,635 | 16.00 | 8.19 | 0 |

All twelve output comparisons (four inputs × three QC files) were byte-identical
to the historical outputs. Real input stays outside Git. Only sanitized metrics,
checksums and public accessions are committed.

## Reproduction

Use the exact checksummed source VCFs above, retain the recorded source manifests,
and set `IMPLEMENTATION_SHA`, `INPUT_VCF` and a fresh `OUTPUT_DIR` for each replay:

```bash
git show "${IMPLEMENTATION_SHA}:bin/summarize_filter_qc.py" > summarize_filter_qc.py
sha256sum "$INPUT_VCF" summarize_filter_qc.py
mkdir "$OUTPUT_DIR"
# Capture MemAvailable, SwapTotal/SwapFree and pswpin/pswpout before and after.
/usr/bin/time -v python3 summarize_filter_qc.py \
  --filtered-vcf "$INPUT_VCF" --cohort-id cohort --stage filtered --variant-type snp \
  --filter-breakdown-output "$OUTPUT_DIR/filter_breakdown.tsv" \
  --annotation-qc-output "$OUTPUT_DIR/annotation_qc.tsv" \
  --summary-output "$OUTPUT_DIR/filter_qc.summary.txt"
```

Run old/new separately with sufficient host memory and fresh output directories;
compare all three files with `cmp` or SHA256. Do not use a memory-capped completion
as an unconstrained RSS measurement. The committed JSON preserves the actual
script hashes, memory ceiling, counters and checksums for this audit.
