# Published output contract (Issue #54)

The table records the pre-refactor process/alias → path/pattern/mode mapping.
All paths are relative to `params.outdir`; every policy now lives in the
`process` block of `nextflow.config`. No `saveAs` rule is introduced or removed.

| Process / alias | Relative path | Pattern | Mode | Policy location |
| --- | --- | --- | --- | --- |
| `BWA_MEM2_INDEX` | `reference` | `*.{0123,amb,ann,bwt.2bit.64,pac}` | copy | module → config |
| `BWA_MEM2_MEM_SORT` | `logs/mapping` | `*.bwa-mem2.log` | copy | module → config |
| `FASTP` | `reads/trimmed` | `*.fastq.gz` | copy | module → config |
| `FASTP` | `qc/fastp` | `*.{json,html}` | copy | module → config |
| `GATK_CREATE_SEQUENCE_DICTIONARY` | `reference` | `*.dict` | copy | module → config |
| `GATK_MARKDUPLICATES` | `qc/markduplicates` | `*.metrics.txt` | copy | module → config |
| `MULTIQC` | `qc/multiqc` | `(all emitted path outputs)` | copy | module → config |
| `SAMTOOLS_FAIDX` | `reference` | `*.fai` | copy | module → config |
| `SAMTOOLS_INDEX` | `alignment` | `*.markdup.bam*` | copy | module → config |
| `SAMTOOLS_QC` | `qc/samtools` | `*.txt` | copy | module → config |
| `FASTQC_RAW` | `qc/fastqc/raw` | `*_fastqc.{html,zip}` | copy | config unchanged |
| `FASTQC_TRIMMED` | `qc/fastqc/trimmed` | `*_fastqc.{html,zip}` | copy | config unchanged |
| `GATK_HAPLOTYPECALLER` | `variants/gvcf` | `*.g.vcf.gz*` | copy | config unchanged |
| `GATK_GATHERVCFS` | `variants/raw` | `*.raw.vcf.gz*` | copy | config unchanged |
| `GATK_SELECTVARIANTS` | `variants/by_type` | `*.vcf.gz*` | copy | config unchanged |
| `GATK_VARIANTFILTRATION` | `variants/filtered` | `*.filtered.vcf.gz*` | copy | config unchanged |
| `GATK_SELECTPASSVARIANTS` | `variants/pass` | `*.pass.vcf.gz*` | copy | config unchanged |
| `BCFTOOLS_STATS` | `qc/variants` | `*.{tsv,txt}` | copy | config unchanged |
| `SUMMARIZE_VARIANT_QC` | `qc/variants` | `*.{tsv,txt}` | copy | config unchanged |
| `SUMMARIZE_FILTER_QC` | `qc/variants` | `*.{tsv,txt}` | copy | config unchanged |
| `RECONCILE_VARIANT_TYPE_COUNTS` | `qc/variants` | `*.{tsv,txt}` | copy | config unchanged |
| `GS_NORMALIZE_VARIANTS` | `variants/gs_normalized` | `*.{normalized.vcf.gz,normalized.vcf.gz.tbi,normalize.report.txt}` | copy | config unchanged |
| `CLASSIFY_NORMALIZED_VARIANTS` | `qc/variants` | `*.{tsv,txt}` | copy | config unchanged |
| `GS_INDEX_CLASSIFIED_VARIANTS` | `variants/gs_classified` | `*.classified.vcf.gz*` | copy | config unchanged |
| `GATK_VARIANTFILTRATION_GS` | `variants/gs_filtered` | `*.filtered.vcf.gz*` | copy | config unchanged |
| `GATK_SELECTPASSVARIANTS_GS` | `variants/gs_pass` | `*.pass.vcf.gz*` | copy | config unchanged |
| `BUILD_GS_PANEL` | `gs_panel` | `*.{tsv,tsv.gz,txt}` | copy | config unchanged |
| `RECONCILE_GS_PANEL_ACCOUNTING` | `gs_panel` | `*.{tsv,txt}` | copy | config unchanged |
| `BUILD_GS_PANEL_MANIFEST` | `gs_panel` | `*.json` | copy | config unchanged |
| `BUILD_RUN_MANIFEST` | `provenance` | `*.json` | copy | config unchanged |

Only the nine module-local policies move. Removing the old `publishDir` statements
leaves every module body byte-identical to the baseline. Removing the newly added
selectors leaves `nextflow.config` byte-identical to the baseline. The workflow,
main entry point and container definitions are unchanged: task count, DAG,
resource labels, optional output declarations and scientific processing are preserved.

`PublishedOutputContract.groovy` checks real copied files in the existing synthetic
E2E, GS-disabled, prebuilt BWA and prebuilt FAI/dict cases. Its expected filenames
were first checked against the pre-refactor #62 outputs. These cases already use
an explicit non-default `outdir`, so the checks cover output relocation as well.
They check both FASTP destinations, raw/trimmed FastQC, four MultiQC outputs,
alignment/index files, reference outputs and the GS lineage. Existing tests also
check genotypes, accounting, metadata and MultiQC content; path checks do not
replace those scientific assertions. Module tests continue to use the repository's
explicit `tests/nextflow.config` and `test,docker` nf-test profile.
