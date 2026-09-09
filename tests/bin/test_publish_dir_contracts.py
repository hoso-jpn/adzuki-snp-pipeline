import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "nextflow.config"
PROCESS_SOURCES = (
    *(ROOT / "modules").rglob("*.nf"),
    *(ROOT / "workflows").rglob("*.nf"),
    *(ROOT / "subworkflows").rglob("*.nf"),
)

MOVED_POLICIES = {
    "BWA_MEM2_INDEX": ("${params.outdir}/reference", "*.{0123,amb,ann,bwt.2bit.64,pac}"),
    "BWA_MEM2_MEM_SORT": ("${params.outdir}/logs/mapping", "*.bwa-mem2.log"),
    "GATK_CREATE_SEQUENCE_DICTIONARY": ("${params.outdir}/reference", "*.dict"),
    "GATK_MARKDUPLICATES": ("${params.outdir}/qc/markduplicates", "*.metrics.txt"),
    "SAMTOOLS_FAIDX": ("${params.outdir}/reference", "*.fai"),
    "SAMTOOLS_INDEX": ("${params.outdir}/alignment", "*.markdup.bam*"),
    "SAMTOOLS_QC": ("${params.outdir}/qc/samtools", "*.txt"),
}


def extract_selector(config: str, process_name: str) -> str:
    marker = f"withName: {process_name} {{"
    start = config.index(marker)
    brace = config.index("{", start)
    depth = 0
    for index in range(brace, len(config)):
        if config[index] == "{":
            depth += 1
        elif config[index] == "}":
            depth -= 1
            if depth == 0:
                return config[start : index + 1]
    raise AssertionError(f"unbalanced config selector: {process_name}")


class PublishDirOwnershipTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = CONFIG.read_text(encoding="utf-8")

    def test_process_sources_do_not_own_publish_policy(self) -> None:
        directive = re.compile(r"(?m)^\s*publishDir\s*[=(]")
        offenders = [
            str(path.relative_to(ROOT))
            for path in PROCESS_SOURCES
            if directive.search(path.read_text(encoding="utf-8"))
        ]
        self.assertEqual([], offenders)

    def test_every_publish_assignment_is_centralized_in_nextflow_config(self) -> None:
        self.assertEqual(29, len(re.findall(r"(?m)^\s*publishDir\s*=", self.config)))

    def test_moved_single_destination_policies_are_unchanged(self) -> None:
        for process_name, (path, pattern) in MOVED_POLICIES.items():
            with self.subTest(process=process_name):
                block = extract_selector(self.config, process_name)
                self.assertIn(f'path   : "{path}"', block)
                self.assertIn("mode   : 'copy'", block)
                self.assertIn(f"pattern: '{pattern}'", block)

    def test_fastp_keeps_both_distinct_destinations(self) -> None:
        block = extract_selector(self.config, "FASTP")
        self.assertEqual(2, block.count("mode   : 'copy'"))
        self.assertIn('path   : "${params.outdir}/reads/trimmed"', block)
        self.assertIn("pattern: '*.fastq.gz'", block)
        self.assertIn('path   : "${params.outdir}/qc/fastp"', block)
        self.assertIn("pattern: '*.{json,html}'", block)

    def test_multiqc_keeps_unfiltered_copy_policy(self) -> None:
        block = extract_selector(self.config, "MULTIQC")
        self.assertIn('path: "${params.outdir}/qc/multiqc"', block)
        self.assertIn("mode: 'copy'", block)
        self.assertNotIn("pattern:", block)


if __name__ == "__main__":
    unittest.main()
