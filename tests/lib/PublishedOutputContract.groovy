import java.nio.file.Files

// Public names observed on the pre-refactor synthetic E2E, independent of
// nextflow.config's selectors. Checking actual copied files catches selectors
// that parse correctly but never match a process or one of its aliases.
class PublishedOutputContract {
    static void names(File root, String directory, Collection expected) {
        def parent = new File(root, directory)
        def children = parent.listFiles() ?: []
        assert children.collect { it.name } as Set == expected as Set :
            "Unexpected published files in ${directory}"
        children.each { child ->
            assert !Files.isSymbolicLink(child.toPath()) : "Expected copy: ${child.name}"
            if (child.isFile()) assert child.length() > 0
        }
    }

    static void verify(String outputDir, boolean gsEnabled = true,
                       boolean generatedBwa = true, boolean generatedFaiDict = true) {
        def root = new File(outputDir)
        def groups = ['sample_a_L001', 'sample_a_L002', 'sample_b_L001']
        def samples = ['sample_a', 'sample_b']
        names(root, 'reads/trimmed', groups.collectMany { group ->
            ["${group}.trimmed_R1.fastq.gz".toString(), "${group}.trimmed_R2.fastq.gz".toString()]
        })
        names(root, 'qc/fastp', groups.collectMany { group ->
            ["${group}.fastp.json".toString(), "${group}.fastp.html".toString()]
        })
        ['raw', 'trimmed'].each { stage ->
            names(root, "qc/fastqc/${stage}", groups.collectMany { group ->
                ['R1', 'R2'].collectMany { mate ->
                    ['html', 'zip'].collect { extension ->
                        "${group}.${stage}.${mate}_fastqc.${extension}".toString()
                    }
                }
            })
        }
        names(root, 'alignment', samples.collectMany { sample ->
            ["${sample}.markdup.bam".toString(), "${sample}.markdup.bam.bai".toString()]
        })
        names(root, 'logs/mapping', groups.collect { "${it}.bwa-mem2.log".toString() })
        names(root, 'qc/markduplicates', samples.collect { "${it}.markduplicates.metrics.txt".toString() })
        names(root, 'qc/samtools', samples.collectMany { sample ->
            ['flagstat', 'idxstats', 'stats'].collect { "${sample}.${it}.txt".toString() }
        })
        names(root, 'qc/multiqc', [
            'multiqc_report.html', 'multiqc_data', 'multiqc_config.yaml', 'multiqc_version.txt'
        ])
        assert new File(root, 'qc/multiqc/multiqc_data/multiqc_data.json').isFile()
        def reference = []
        if (generatedBwa) reference.addAll(['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'].collect {
            "synthetic.fa.${it}".toString()
        })
        if (generatedFaiDict) reference.addAll(['synthetic.fa.fai', 'synthetic.dict'])
        names(root, 'reference', reference)
        def gsFiles = [
            'gs_normalized': ['cohort_gs.normalized.vcf.gz', 'cohort_gs.normalized.vcf.gz.tbi',
                              'cohort_gs.normalize.report.txt'],
            'gs_classified': ['cohort_gs.classified.vcf.gz', 'cohort_gs.classified.vcf.gz.tbi'],
            'gs_filtered': ['cohort_gs.snp.filtered.vcf.gz', 'cohort_gs.snp.filtered.vcf.gz.tbi'],
            'gs_pass': ['cohort_gs.snp.pass.vcf.gz', 'cohort_gs.snp.pass.vcf.gz.tbi']
        ]
        gsFiles.each { directory, expected ->
            names(root, "variants/${directory}", gsEnabled ? expected : [])
        }
        names(root, 'gs_panel', gsEnabled ? [
            'genotype_encoding_accounting.summary.txt', 'genotype_encoding_accounting.tsv',
            'genotype_matrix.tsv.gz', 'manifest.json', 'record_accounting.summary.txt',
            'record_accounting.tsv', 'sample_metadata.tsv', 'variant_metadata.tsv'
        ].collect { "cohort.gs_panel.${it}".toString() } : [])
    }
}
