import gzip
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks/issue45"))
import validate_generated_cohort as validation  # noqa: E402


def header(sample="sample_a", extra=""):
    return (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller '
        "--input sample.bam --output sample.g.vcf.gz --reference ref.fa "
        "--sample-ploidy 2 --native-pair-hmm-threads 4 --emit-ref-confidence GVCF "
        '--create-output-variant-index true",Version="4.6.2.0",Date="test">\n'
        + extra
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + sample
        + "\n"
    )


class CohortValidationTests(unittest.TestCase):
    def read(self, text, expected="sample_a"):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.g.vcf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(text)
            return validation.read_header(path, expected, [("chr1", 100)])

    def test_sample_header_and_ploidy_mismatch_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "expected sample"):
            self.read(header("different_sample"))
        with self.assertRaisesRegex(ValueError, "sample-ploidy"):
            self.read(header().replace("--sample-ploidy 2", "--sample-ploidy 1"))

    def test_parameter_fingerprint_ignores_verified_filenames_but_not_defaults(self):
        first = self.read(header())["haplotypecaller_parameters_sha256"]
        renamed = header().replace("sample.bam", "another.bam")
        self.assertEqual(first, self.read(renamed)["haplotypecaller_parameters_sha256"])
        changed = header().replace(
            "HaplotypeCaller --input", "HaplotypeCaller --min-base-quality-score 30 --input"
        )
        self.assertNotEqual(first, self.read(changed)["haplotypecaller_parameters_sha256"])

    def test_shared_header_definition_change_is_detected(self):
        first = self.read(header())["shared_header_definitions_sha256"]
        second = self.read(header().replace("Type=String", "Type=Integer"))
        self.assertNotEqual(first, second["shared_header_definitions_sha256"])

    def test_record_validation_rejects_order_ploidy_and_allele_errors(self):
        row = b"chr1\t20\t.\tA\t<NON_REF>\t.\t.\tEND=30\tGT:DP\t0/0:10\n"
        result = validation.stream_records(io.BytesIO(row), [("chr1", 100)])
        self.assertEqual(result["record_count"], 1)
        for bad in (
            row.replace(b"0/0", b"0"),
            row.replace(b"0/0", b"0/2"),
            row.replace(b"END=30", b"END=101"),
            row + row.replace(b"\t20\t", b"\t10\t"),
        ):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                validation.stream_records(io.BytesIO(bad), [("chr1", 100)])

    def test_ordinary_gzip_is_not_accepted_as_complete_bgzf(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ordinary.gz"
            path.write_bytes(gzip.compress(header().encode()))
            with self.assertRaisesRegex(ValueError, "BGZF EOF"):
                validation.verify_bgzf(path)


if __name__ == "__main__":
    unittest.main()
