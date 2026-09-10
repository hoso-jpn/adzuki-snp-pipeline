import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks/issue45"))
from benchmark_tools import batches_from_log, compare_callsets  # noqa: E402
from execute_experiments import validate_tiling  # noqa: E402


class BenchmarkExecutionTests(unittest.TestCase):
    def test_qd_jitter_exception_does_not_hide_other_semantic_changes(self):
        row = "chr1\t1\t.\tA\tC\t300\t.\tAC=2;AN=2;QD=30;MQ=60\tGT:AD:DP\t1/1:0,6:6\n"
        with tempfile.TemporaryDirectory() as tmp:
            left, right = Path(tmp) / "left.gz", Path(tmp) / "right.gz"
            for source, changed, expected in (
                (row, row.replace("QD=30", "QD=27"), True),
                (row, row.replace("QD=30", "QD=1"), False),
                (row, row.replace("MQ=60", "MQ=50"), False),
                (row, row.replace("1/1", "0/1"), False),
                (
                    row.replace("\t300\t", "\t60\t"),
                    row.replace("\t300\t", "\t60\t").replace("QD=30", "QD=27"),
                    False,
                ),
            ):
                left.write_bytes(gzip.compress(source.encode()))
                right.write_bytes(gzip.compress(changed.encode()))
                audit = compare_callsets(left, right, [("chr1", 100)])
                self.assertEqual(
                    expected,
                    audit[
                        "all_differences_explained_by_high_qd_jitter_with_unchanged_current_filter"
                    ],
                )

    def test_51_sample_log_requires_both_real_batches_and_completion(self):
        complete = (
            "Importing batch 1 with 50 samples\nDone importing batch 1/2\n"
            "Importing batch 2 with 1 samples\nDone importing batch 2/2\n"
        )
        self.assertEqual(batches_from_log(complete, 51, 50)["observed"], [(1, 50), (2, 1)])
        for text in (
            "sample_count=51 batch_size=50",
            complete.replace("Done importing batch 2/2", ""),
            complete.replace("batch 2 with 1 samples", "batch 2 with 50 samples"),
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                batches_from_log(text, 51, 50)

    def test_plan_rejects_gap_overlap_order_and_incomplete_reference(self):
        valid = [{"intervals": ["chr1:1-5"]}, {"intervals": ["chr1:6-10", "small1", "small2"]}]
        contigs = [("chr1", 10), ("small1", 2), ("small2", 3)]
        validate_tiling(valid, contigs)
        for original, replacement in (
            ("6-10", "7-10"),
            ("6-10", "5-10"),
            ('"small1", "small2"', '"small2", "small1"'),
            (', "small2"', ""),
        ):
            invalid = json.loads(json.dumps(valid).replace(original, replacement))
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                validate_tiling(invalid, contigs)


if __name__ == "__main__":
    unittest.main()
