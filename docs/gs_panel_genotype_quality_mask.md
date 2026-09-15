# GS genotype quality mask: real-cohort evidence (Issue #64)

The contract for the optional genotype-level DP/GQ mask is in
[`gs_panel_data_contract.md`](gs_panel_data_contract.md#genotype-quality-mask-issue-64-optional-off-by-default).
This document records what it was measured to do on real public data. **It
is a sensitivity analysis, not an accuracy result.** There is no truth set
for this cohort, so nothing below says which threshold yields better
genotypes, and no threshold is adopted as a default.

## Input and lineage

- **Cohort**: Issue #45's 51 unique BioSamples (PRJNA1138464, WGS, Illumina
  HiSeq X), joint-genotyped with the production per-contig plan (suite E1,
  record-identical to the production baseline). Joint callset SHA256
  `e83b6531…`.
- **GS lineage**: rebuilt from that callset with the production GS modules
  and hard-filter definitions of a clean `main@7ef04cb`
  (`benchmarks/issue64/gs_lineage_from_callset.nf.template`). No FASTQ
  processing and no variant calling.
- **Scale**: 13,219,170 raw/all records → 9,746,661 GS-eligible PASS SNPs
  × 51 samples = **497,079,711 genotype cells**.
- **Branch code**: `2f5ca34`, run in the pinned `python:3.12` and
  `bcftools:1.24` containers on seedcore-01, one step at a time.
- **Inputs never modified**: the PASS VCF's SHA256 was `f1ab1c3a…` before
  and after all checks.

Machine-readable records:
[`evidence/issue64/real_cohort_checks.json`](evidence/issue64/real_cohort_checks.json)
and
[`evidence/issue64/real_cohort_genotype_quality_sensitivity.json`](evidence/issue64/real_cohort_genotype_quality_sensitivity.json).

## 1. Mask off: byte-identical to main

The branch builder with no quality option reproduced main's matrix, sample
metadata, variant metadata, genotype accounting and its summary **byte for
byte** (all five SHA256 equal). It ran in 700.7 s at 28.4 MiB peak RSS, the
same as main's 679.5 s at 26.6 MiB.

## 2. Manifest: only the intended v3 change

main's and the branch's GS manifest builders, given the same unmasked
panel and arguments, agree on every field and every checksum. The only
differences are the per-run identifiers and the intended change: v2 → v3
plus a `genotype_quality_mask` block with `enabled: false`.

## 3. The default `reject` policy stops on this cohort

With `min_dp=10`, `min_gq=20` and every unreadable-value policy left at
`reject`, the build failed after 5.3 s, publishing nothing:

```
line 26705: NC_068970.1:1993424 sample 'SRR29909068': DP is value_missing and the policy rejects it (missing_value=reject)
```

GATK wrote a called genotype with `FORMAT/DP = .`. Across the panel there are
4,590 such calls, and no missing GQ, truncated, absent or malformed value.
The fail-closed default does what it is meant to do: a run on this cohort
must state `gs_genotype_missing_value` explicitly rather than have those
calls silently kept or dropped.

## 4. One illustrative enabled build, verified

`min_dp=10`, `min_gq=20`, all unreadable-value policies `unevaluated`.
These values were chosen only to exercise the enabled path; they are not a
recommendation.

| | Value |
| --- | ---: |
| evaluated calls | 465,500,816 |
| passed / kept unevaluated / masked | 254,580,845 / 4,412 / **210,915,559 (45.3%)** |
| masked by original dosage (hom-ref / het / hom-alt) | 185,203,732 / 6,242,066 / 19,469,761 |
| masked: DP low only / GQ low only / both | 80,590,206 / 7,462,793 / 122,862,560 |
| cells treated as missing | 31,578,895 → 242,494,454 (6.35% → 48.78%) |
| rows with at least one masked call | 9,743,049 of 9,746,661 |
| verifier | `consistent`: all 497,079,711 cells checked |
| build / verify peak RSS | 29.8 MiB / 24.7 MiB |
| build / verify / index wall time | 1,907 s / 1,965 s / 12 s |
| masked VCF | 2.38 GB, BGZF, indexed by bcftools |

Memory stays flat, as the bounded-memory contract requires. Time does not:
the enabled build takes 2.7× the disabled one, and verification adds about
the same again (see Risks).

## 5. Sensitivity over a threshold grid

One streaming pass (1,851 s, 24.4 MiB) evaluated 6 DP × 3 GQ thresholds plus
"no threshold" on each field. Unreadable values are kept here, so masking
follows the thresholds alone. The `dp>=10,gq>=20` point reproduces the build
in step 4 exactly on all nine cross-checked counts: total and evaluated calls,
masked total, masked by dosage, and each reason.

Evaluated calls' DP distribution: DP≤2: 4.8%; DP≤4: 13.9%; DP≤7: 31.3%; DP≤9: 43.7%; DP≤14: 71.8%; DP≤19: 89.4%; DP≤29: 98.7%. The cohort's per-call
depth is low, which is the context for every number below.

Columns: masked calls (% of evaluated) · reason split · overall post-mask
missingness · per-sample missingness min / median / max · variants whose AC
changed · whose AN changed · that went from polymorphic to monomorphic ·
whose AN became 0 · mean AN (2 × 51 = 102 at most).

| min DP | min GQ | masked | DP only | GQ only | DP+GQ | missing | sample min/med/max | AC changed | AN changed | poly→mono | AN=0 | mean AN |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| — | — | 0 (0.0%) | 0 | 0 | 0 | 6.35% | 0.028 / 0.046 / 0.196 | 0 | 0 | 0 | 0 | 95.5 |
| — | 10 | 48,050,395 (10.3%) | 0 | 48,050,395 | 0 | 16.02% | 0.045 / 0.088 / 0.542 | 4,526,832 | 9,559,473 | 1,052,312 | 2,639 | 85.7 |
| — | 20 | 130,325,353 (28.0%) | 0 | 130,325,353 | 0 | 32.57% | 0.069 / 0.207 / 0.839 | 6,955,015 | 9,726,473 | 3,130,806 | 16,489 | 68.8 |
| — | 30 | 221,273,562 (47.5%) | 0 | 221,273,562 | 0 | 50.87% | 0.099 / 0.460 / 0.949 | 7,743,251 | 9,739,805 | 4,830,719 | 49,790 | 50.1 |
| 3 | — | 22,304,595 (4.8%) | 22,304,595 | 0 | 0 | 10.84% | 0.031 / 0.061 / 0.392 | 2,521,474 | 8,341,576 | 471,360 | 857 | 90.9 |
| 3 | 10 | 48,129,963 (10.3%) | 79,568 | 25,825,368 | 22,225,027 | 16.04% | 0.045 / 0.088 / 0.543 | 4,560,269 | 9,561,720 | 1,058,576 | 2,741 | 85.6 |
| 3 | 20 | 130,398,723 (28.0%) | 73,370 | 108,094,128 | 22,231,225 | 32.59% | 0.069 / 0.207 / 0.840 | 6,975,979 | 9,726,804 | 3,140,643 | 16,854 | 68.8 |
| 3 | 30 | 221,341,253 (47.5%) | 67,691 | 199,036,658 | 22,236,904 | 50.88% | 0.099 / 0.460 / 0.949 | 7,757,797 | 9,739,927 | 4,843,955 | 50,427 | 50.1 |
| 5 | — | 64,625,259 (13.9%) | 64,625,259 | 0 | 0 | 19.35% | 0.037 / 0.095 / 0.670 | 5,649,502 | 9,621,945 | 1,811,362 | 4,597 | 82.3 |
| 5 | 10 | 69,709,475 (15.0%) | 21,659,080 | 5,084,216 | 42,966,179 | 20.38% | 0.048 / 0.109 / 0.674 | 5,962,259 | 9,680,276 | 1,842,770 | 5,451 | 81.2 |
| 5 | 20 | 131,145,699 (28.2%) | 820,346 | 66,520,440 | 63,804,913 | 32.74% | 0.070 / 0.208 / 0.850 | 7,246,540 | 9,729,912 | 3,324,067 | 18,101 | 68.6 |
| 5 | 30 | 222,054,759 (47.7%) | 781,197 | 157,429,500 | 63,844,062 | 51.02% | 0.099 / 0.460 / 0.954 | 7,973,606 | 9,741,117 | 5,058,507 | 53,493 | 50.0 |
| 8 | — | 145,722,707 (31.3%) | 145,722,707 | 0 | 0 | 35.67% | 0.053 / 0.229 / 0.905 | 7,869,772 | 9,727,378 | 4,374,343 | 23,345 | 65.6 |
| 8 | 10 | 149,720,854 (32.2%) | 101,670,459 | 3,998,147 | 44,052,248 | 36.47% | 0.061 / 0.239 / 0.907 | 8,014,569 | 9,735,926 | 4,421,397 | 24,331 | 64.8 |
| 8 | 20 | 155,714,985 (33.5%) | 25,389,632 | 9,992,278 | 120,333,075 | 37.68% | 0.073 / 0.252 / 0.910 | 8,128,558 | 9,738,795 | 4,502,227 | 26,857 | 63.6 |
| 8 | 30 | 223,818,416 (48.1%) | 2,544,854 | 78,095,709 | 143,177,853 | 51.38% | 0.099 / 0.460 / 0.962 | 8,456,429 | 9,743,328 | 5,537,888 | 57,375 | 49.6 |
| 10 | — | 203,452,766 (43.7%) | 203,452,766 | 0 | 0 | 47.28% | 0.063 / 0.395 / 0.961 | 8,579,446 | 9,738,431 | 5,717,288 | 47,873 | 53.8 |
| 10 | 10 | 207,011,467 (44.5%) | 158,961,072 | 3,558,701 | 44,491,694 | 48.00% | 0.072 / 0.404 / 0.962 | 8,662,949 | 9,741,883 | 5,775,518 | 48,819 | 53.0 |
| 10 | 20 | 210,915,559 (45.3%) | 80,590,206 | 7,462,793 | 122,862,560 | 48.78% | 0.083 / 0.414 / 0.963 | 8,718,805 | 9,743,049 | 5,844,549 | 51,283 | 52.2 |
| 10 | 30 | 225,290,234 (48.4%) | 4,016,672 | 21,837,468 | 199,436,094 | 51.68% | 0.099 / 0.462 / 0.967 | 8,746,912 | 9,744,048 | 5,872,369 | 59,059 | 49.3 |
| 15 | — | 334,320,795 (71.8%) | 334,320,795 | 0 | 0 | 73.61% | 0.141 / 0.800 / 0.995 | 9,405,193 | 9,744,834 | 7,459,028 | 186,796 | 26.9 |
| 15 | 10 | 336,790,299 (72.4%) | 288,739,904 | 2,469,504 | 45,580,891 | 74.11% | 0.149 / 0.805 / 0.995 | 9,415,885 | 9,745,396 | 7,486,236 | 190,993 | 26.4 |
| 15 | 20 | 338,875,081 (72.8%) | 208,549,728 | 4,554,286 | 125,771,067 | 74.53% | 0.159 / 0.808 / 0.995 | 9,425,113 | 9,745,640 | 7,511,341 | 194,756 | 26.0 |
| 15 | 30 | 340,810,696 (73.2%) | 119,537,134 | 6,489,901 | 214,783,661 | 74.92% | 0.170 / 0.812 / 0.996 | 9,434,875 | 9,745,781 | 7,547,272 | 198,826 | 25.6 |
| 20 | — | 415,941,634 (89.4%) | 415,941,634 | 0 | 0 | 90.03% | 0.373 / 0.954 / 0.998 | 9,619,323 | 9,745,930 | 8,337,811 | 767,231 | 10.2 |
| 20 | 10 | 416,281,366 (89.4%) | 368,230,971 | 339,732 | 47,710,663 | 90.10% | 0.377 / 0.954 / 0.999 | 9,621,264 | 9,746,058 | 8,351,746 | 770,808 | 10.1 |
| 20 | 20 | 417,164,942 (89.6%) | 286,839,589 | 1,223,308 | 129,102,045 | 90.28% | 0.384 / 0.955 / 0.999 | 9,624,304 | 9,746,133 | 8,384,047 | 786,321 | 9.9 |
| 20 | 30 | 418,431,470 (89.9%) | 197,157,908 | 2,489,836 | 218,783,726 | 90.53% | 0.394 / 0.957 / 0.999 | 9,628,332 | 9,746,175 | 8,444,895 | 817,342 | 9.7 |

Originally missing (before any mask): 31,578,895 calls (6.35%); no
non-diploid or non-biallelic-index calls.

### Per-variant missingness (number of variants per bin)

| policy | ==0.0 | (0.0,0.05] | (0.05,0.1] | (0.1,0.2] | (0.2,0.3] | (0.3,0.5] | (0.5,0.8] | (0.8,1.0] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `dp>=none,gq>=none` | 2,885,873 | 3,644,699 | 1,634,372 | 855,730 | 327,828 | 245,058 | 131,502 | 21,599 |
| `dp>=3,gq>=none` | 397,654 | 2,556,828 | 3,838,467 | 1,847,538 | 518,787 | 355,354 | 184,825 | 47,208 |
| `dp>=5,gq>=10` | 31,338 | 178,848 | 1,323,894 | 5,012,880 | 1,858,053 | 920,519 | 323,845 | 97,284 |
| `dp>=none,gq>=20` | 10,635 | 50,189 | 153,133 | 1,107,620 | 3,924,293 | 3,526,259 | 793,526 | 181,006 |
| `dp>=10,gq>=20` | 2,430 | 10,730 | 39,693 | 166,685 | 375,288 | 5,312,431 | 3,403,860 | 435,544 |

### Allele frequency

|AF change| after masking, and variants crossing MAF 0.05:

| policy | Δ=0 | (0, 0.01] | (0.01, 0.05] | > 0.05 | MAF ≥0.05 → <0.05 | <0.05 → ≥0.05 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `dp>=3,gq>=none` | 1,411,304 | 5,722,248 | 2,469,320 | 142,932 | 455,452 | 87,951 |
| `dp>=5,gq>=10` | 86,371 | 4,030,593 | 4,798,628 | 825,618 | 1,294,772 | 142,662 |
| `dp>=none,gq>=20` | 45,455 | 2,730,995 | 5,257,388 | 1,696,334 | 1,680,050 | 246,085 |
| `dp>=10,gq>=20` | 39,615 | 1,440,477 | 5,663,749 | 2,551,537 | 2,285,393 | 169,961 |

### Per sample

Masking does not fall evenly across samples. The table lists the eight
samples with the highest missingness before masking, then the three least
affected by the illustrative policy. The first group includes the highly
divergent samples Issue #45 identified, and masking removes far more of their
calls.

| sample | no mask | DP≥5, GQ≥10 | DP≥10, GQ≥20 |
| --- | ---: | ---: | ---: |
| SRR29909074 | 0.196 | 0.670 | 0.955 |
| SRR29909070 | 0.195 | 0.674 | 0.958 |
| SRR29909072 | 0.181 | 0.657 | 0.957 |
| SRR29909073 | 0.166 | 0.599 | 0.943 |
| SRR29909068 | 0.128 | 0.540 | 0.932 |
| SRR29908888 | 0.118 | 0.321 | 0.780 |
| SRR29909423 | 0.106 | 0.338 | 0.831 |
| SRR29909067 | 0.101 | 0.614 | 0.963 |
| … | | | |
| SRR29908806 | 0.034 | 0.052 | 0.083 |
| SRR29908807 | 0.036 | 0.058 | 0.108 |
| SRR29908805 | 0.038 | 0.060 | 0.108 |

## Reading

- Thresholds often quoted for higher-depth data remove a very large share
  of this cohort's calls. DP≥10 alone masks 43.7% of evaluated calls; with
  GQ≥20 added, 45.3%. Most masked calls are hom-ref, and GQ masks mainly
  calls that DP already masks.
- Even DP≥3 changes AN at 8.3 M of 9.7 M variants. Moderate GQ thresholds
  move many variants across MAF 0.05, mostly downwards. A mask therefore
  changes downstream MAF selection, not only missingness.
- The effect is concentrated in specific samples, so a mask interacts with
  sample-level QC (out of scope here). Composition affects the outcome as
  much as the thresholds do.
- None of this measures genotype accuracy. A call masked here may have been
  correct, and a call kept may be wrong.

## Adoption conditions (recorded, not met)

No threshold is adopted, and the default stays `gs_genotype_quality_mask =
false` with no threshold. Before any threshold becomes a production
default or a recommended setting, all of these would be needed:

1. **An accuracy reference**: replicate concordance, pedigree/Mendelian
   consistency, or an orthogonal genotyping comparison, to show what a
   threshold does to genotype error, not just to call counts.
2. **Downstream impact**: the effect on the GS model the panel feeds
   (imputation, relationship matrices, prediction), measured with its
   consumer, because masking at these rates changes MAF and missingness
   materially.
3. **A decision on sample-level QC first**: low-depth and divergent samples
   dominate the masking, and sample exclusion is a separate stage.
4. **An explicit `missing_value` policy**: GATK emits called genotypes with
   `DP=.`, so the fail-closed default stops on real data by design.
5. **Re-measured resource limits at the target cohort size**: see Risks.

## Risks and limits

- **Time at larger cohorts.** At 51 samples the enabled build took 32 min
  and verification 33 min. `BUILD_GS_PANEL` (`process_gs_panel`) and
  `VERIFY_GS_GENOTYPE_QUALITY_MASK` (`process_low`) both allow 2 h. Both
  read every cell once, so assuming their time grows roughly with cell count
  (an assumption, not measured here), a cohort several times larger is
  expected to exceed that with the mask enabled. The limits were not raised without
  a measurement at that size; the mask is off by default and the 327-sample
  run is NO-GO (Issue #45).
- A 51-sample accession-order subset of one BioProject, at one depth
  profile. Other cohorts will differ.
- The masked VCF recomputes only AC, AN and AF. Other INFO annotations
  still describe the original calls.
- `value_malformed` never occurred and cannot come out of htslib for an
  Integer field; it is covered only by synthetic tests.
