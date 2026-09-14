"""A representational change must never hide a semantic one in the E3 profile."""

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HELPERS = Path(__file__).resolve().parents[2] / "benchmarks/issue45"

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def row(pos="100", alt="C", qual="500.0", filt=".", info="AC=1;AN=4;QD=30.0", gts=("0/1", "0/0")):
    samples = "\t".join(f"{gt}:10,5:15" for gt in gts)
    return f"chr1\t{pos}\t.\tA\t{alt}\t{qual}\t{filt}\t{info}\tGT:AD:DP\t{samples}\n"


class DifferenceProfileTests(unittest.TestCase):
    def profile(self, left_rows, right_rows):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        base = Path(tmp.name)
        for name, rows in (("left.vcf.gz", left_rows), ("right.vcf.gz", right_rows)):
            (base / name).write_bytes(gzip.compress((HEADER + "".join(rows)).encode()))
        out = base / "profile.json"
        result = subprocess.run(
            [
                sys.executable,
                str(HELPERS / "profile_callset_differences.py"),
                "--left",
                str(base / "left.vcf.gz"),
                "--right",
                str(base / "right.vcf.gz"),
                "--left-label",
                "E2b",
                "--right-label",
                "E3",
                "--summary-json",
                str(out),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(0, result.returncode, result.stderr)
        return json.loads(out.read_text())

    def test_annotation_only_change_is_not_reported_as_semantic(self):
        p = self.profile([row()], [row(info="AC=1;AN=4;QD=27.5")])
        self.assertTrue(p["no_semantic_difference_observed"])
        self.assertEqual(1, p["totals"]["different_records"])
        self.assertIn("QD", p["annotation_differences"]["info_fields"])
        self.assertEqual(-2.5, p["info_numeric_deltas"]["QD"]["min"])

    def test_a_genotype_change_is_always_semantic(self):
        p = self.profile([row()], [row(gts=("1/1", "0/0"))])
        self.assertFalse(p["no_semantic_difference_observed"])
        self.assertEqual(1, p["scientific_semantic_categories"]["records_with_genotype_change"])
        self.assertIn("0/1->1/1", p["genotype_transitions"])

    def test_a_filter_change_is_always_semantic(self):
        p = self.profile([row()], [row(filt="LowQD")])
        self.assertFalse(p["no_semantic_difference_observed"])
        self.assertEqual(1, p["scientific_semantic_categories"]["records_with_filter_change"])

    def test_a_variant_appearing_or_disappearing_is_always_semantic(self):
        gained = self.profile([row()], [row(), row(pos="200")])
        self.assertFalse(gained["no_semantic_difference_observed"])
        self.assertEqual(1, gained["totals"]["right_only_variants"])
        lost = self.profile([row(), row(pos="200")], [row()])
        self.assertFalse(lost["no_semantic_difference_observed"])
        self.assertEqual(1, lost["totals"]["left_only_variants"])

    def test_an_added_or_removed_annotation_is_reported_distinctly(self):
        p = self.profile([row()], [row(info="AC=1;AN=4")])
        self.assertIn("QD (removed)", p["annotation_differences"]["info_fields"])

    def test_identical_callsets_report_no_difference_at_all(self):
        p = self.profile([row(), row(pos="200")], [row(), row(pos="200")])
        self.assertTrue(p["no_semantic_difference_observed"])
        self.assertEqual(0, p["totals"]["different_records"])
        self.assertEqual(2, p["totals"]["identical_records"])


if __name__ == "__main__":
    unittest.main()
