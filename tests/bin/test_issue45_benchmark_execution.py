import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from threading import Event

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks/issue45"))
import execute_experiments  # noqa: E402
from benchmark_tools import (  # noqa: E402
    ToolRunner,
    batches_from_log,
    compare_callsets,
)
from execute_experiments import (  # noqa: E402
    COMPLETED,
    completed_experiment,
    validate_tiling,
)


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
                (row.replace(";QD=30", ""), row, False),
                (row.replace("QD=30", "QD=NaN"), row, False),
                (
                    row.replace("0,6", ".,."),
                    row.replace("0,6", ".,.").replace("QD=30", "QD=27"),
                    False,
                ),
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


class ResumeContractTests(unittest.TestCase):
    """A multi-hour suite must reuse finished work without ever overwriting failure evidence."""

    def test_completed_experiment_is_reused_and_unfinished_evidence_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.assertIsNone(completed_experiment(root, "E0"))

            done = root / "E1"
            done.mkdir()
            payload = {"status": COMPLETED, "integrity": {"variant_count": 7}}
            (done / "experiment_result.private.json").write_text(json.dumps(payload))
            self.assertEqual(payload, completed_experiment(root, "E1"))

            failed = root / "E2a"
            failed.mkdir()
            (failed / "experiment_result.private.json").write_text(
                json.dumps({"status": "FAILED", "error_type": "RuntimeError"})
            )
            with self.assertRaisesRegex(ValueError, "retained FAILED evidence"):
                completed_experiment(root, "E2a")

            (root / "E2b").mkdir()
            with self.assertRaisesRegex(ValueError, "without a result"):
                completed_experiment(root, "E2b")


class ResourcePolicyTests(unittest.TestCase):
    """Issue #45: a 15 GiB heap inside a 16 GiB container OOM-killed the first real E0."""

    def runner(self, memory_limit=None):
        return ToolRunner(
            Path("/bench"), Path("/input"), Path("/reference"), Event(), memory_limit=memory_limit
        )

    def test_every_java_heap_stays_below_its_own_container_limit(self):
        runner = self.runner()
        for memory, reserve in (
            (
                execute_experiments.GENOTYPE_MEMORY_GIB,
                execute_experiments.GENOTYPE_NATIVE_RESERVE_GIB,
            ),
            (execute_experiments.GATHER_MEMORY_GIB, execute_experiments.GATHER_NATIVE_RESERVE_GIB),
            (
                execute_experiments.REBLOCK_MEMORY_GIB,
                execute_experiments.REBLOCK_NATIVE_RESERVE_GIB,
            ),
        ):
            with self.subTest(memory=memory):
                self.assertLess(runner.java_heap_gib(memory, reserve), memory)

    def test_a_capped_synthetic_allocation_also_shrinks_the_heap_it_contains(self):
        capped = self.runner(memory_limit=4)
        self.assertEqual(4, capped.effective_memory_gib(execute_experiments.GENOTYPE_MEMORY_GIB))
        heap = capped.java_heap_gib(
            execute_experiments.GENOTYPE_MEMORY_GIB, execute_experiments.GENOTYPE_NATIVE_RESERVE_GIB
        )
        self.assertGreaterEqual(heap, 1)
        self.assertLessEqual(heap, 4)

    def test_concurrent_tasks_fit_host_memory_at_the_declared_tiers(self):
        peak = (
            execute_experiments.MAXIMUM_CONCURRENT_TASKS * execute_experiments.GENOTYPE_MEMORY_GIB
        )
        self.assertLessEqual(peak, 110, "the launch gate requires 110 GiB available at start")


if __name__ == "__main__":
    unittest.main()
