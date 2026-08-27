# CLASSIFY_NORMALIZED_VARIANTS streaming benchmark (Issue #35)

## Scope and controls

Formal targeted replays ran on `seedcore-01` on 2026-08-27 against the detached,
clean benchmark worktree at implementation SHA
`70ae4d638574043e858a6f3b94082a6d5fb767be`. Production's exact Python image was used:
`python:3.12@sha256:dd4fe98ab39f91e936f8e7e7a65a3ce59ecfb11e32f9a125b3132779920ba7f7`.
The pre-run `docker stats`/top-RSS snapshot found no heavy workload; Dify, Open WebUI,
PostgreSQL, and the other lightweight resident services were left running. vLLM was stopped.

Peak measurements deliberately distinguish:

- classifier Python RSS: `VmHWM` of the container's Python main PID;
- container peak: cgroup v2 `memory.peak`, including charged filesystem page cache;
- host pressure: minimum `/proc/meminfo` `MemAvailable`;
- host swap: used swap sampled before, throughout, and after the classifier.

An initial exit-127 attempt (the image has no `/usr/bin/time`) and two measurement/invocation
shakedowns were excluded. The final runs below poll all four measurements every 0.2 seconds
and use the production cohort ID `cohort_gs`.

## Inputs and provenance

Both inputs are the already-published `GS_NORMALIZE_VARIANTS` artifacts from Issue #33's
real Longxiaodou 4 cohort runs. Their source cohorts and upstream checksums are recorded in
`real_cohort_scale_validation_10sample_manifest.json` and
`real_cohort_scale_validation_20sample_manifest.json`.

| Scale | Input artifact | SHA256 | Input records | Samples |
| --- | --- | --- | ---: | ---: |
| 10 | `runs/10sample/results/variants/gs_normalized/cohort_gs.normalized.vcf.gz` | `421cd9f6ac6ce39dc09ce34314a791fe3e6c19956056ed60ed411f414f1e6394` | 10,787,857 | 10 |
| 20 | `runs/20sample/results/variants/gs_normalized/cohort_gs.normalized.vcf.gz` | `bff0d959554dbf5660ff208b01eb4e2e5b11ff6f9bd7b14b5bd37af1bbebf03a` | 12,878,591 | 20 |

## Formal results

| Measurement | 10 samples | 20 samples |
| --- | ---: | ---: |
| Start / end (JST) | 17:43:02 / 17:43:39 | 17:44:03 / 17:44:55 |
| Wall time | 37 s | 52 s |
| Classifier Python peak RSS | 21,328 KiB (20.83 MiB) | 21,400 KiB (20.90 MiB) |
| Container/cgroup peak | 4,867,792,896 B (4.53 GiB) | 7,976,562,688 B (7.43 GiB) |
| Host `MemAvailable` minimum | 121,754,228 kB | 122,611,816 kB |
| Swap before / peak / after | 1,191,972 / 1,191,972 / 1,191,972 kB | 1,191,972 / 1,191,972 / 1,191,972 kB |
| Swap delta | **0 kB** | **0 kB** |
| Exit / OOM | 0 / false | 0 / false |
| Output records | 8,920,725 | 10,556,952 |

The 72 KiB RSS increase from 10 to 20 samples is not cohort-proportional. The implementation
retains the current `(CHROM, POS)` locus only, so memory is bounded by local locus width.
The much larger cgroup readings track the 4.0/7.7 GB plain output files: Linux charges
reclaimable write page cache to the container cgroup. This is why the resource contract uses
the end-to-end cgroup peak plus headroom even though the Issue #35 Python-RSS target is met.

## Output hashes and equivalence

| Artifact | 10-sample SHA256 | 20-sample SHA256 |
| --- | --- | --- |
| Classified plain VCF | `52938c14ce303e325694a9b667e4050e0a487a3880da9e6ef08e6b9ea9848450` | `52ffccc85b749df4f23b8a07f8a067a39a3c4e66be0fbcc90db3dfd57e615642` |
| Accounting TSV | `cd8e803a5cdb64ebeb789ab9c4e620d4341052d799aa1dd2f43781222fccc758` | `db2fadec2c9c75a9c2deff17e8d966241fb1196b3e4036a0c8d9b974a7f846f2` |
| Summary | `b15a598a03baab201bb4bbcf0593f3161ad342e89df4f216f0bbe52a13426e65` | `598ed086fdd3ecae27252ea72c6ed381614cbdd69fc1eacb5fdb8ca415868232` |

At both scales, the formal plain VCF is byte-identical to the historical pre-index targeted
replay. Accounting and summary are byte-identical to the published production files. The
published classified VCF is bgzip-compressed by `GS_INDEX_CLASSIFIED_VARIANTS`; after decode,
bcftools has added only the standard PASS header plus `bcftools_viewVersion` and
`bcftools_viewCommand`. Removing those three compression provenance headers makes each
decoded published VCF byte-identical to its formal plain VCF, including every variant record.

## Resource and retry decision

`process_variant_classification` is set to 10 GiB: about 25% headroom over the larger 7.43
GiB cgroup peak and orders of magnitude over the 20.90 MiB classifier working set. This is
not a mechanical 72-to-4 GiB change. The fixed fail-fast/no-retry policy remains: an OOM at
a new scale or pathological locus should trigger evidence-preserving re-measurement rather
than silently changing the contract on retry.
