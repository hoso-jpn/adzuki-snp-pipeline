"""The public Issue45 evidence must carry the required metrics and no host layout."""

import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

HELPERS = Path(__file__).resolve().parents[2] / "benchmarks/issue45"
sys.path.insert(0, str(HELPERS))
import summarize_evidence  # noqa: E402

COMPLETED = "COMPLETED_AWAITING_COMPARATIVE_REVIEW"
SECRET_PATH = "/data/home/someone-private/adzuki-run"


def metrics(wall, rss):
    return {
        "wall_seconds": wall,
        "peak_rss_bytes": rss,
        "user_cpu_seconds": wall / 2,
        "system_cpu_seconds": wall / 4,
        "attempt": 1,
        "retry_count": 0,
    }


def summary(wall, rss):
    return {
        "task_count": 1,
        "wall_seconds_sum": wall,
        "peak_rss_bytes_max": rss,
        "user_cpu_seconds_sum": wall / 2,
        "system_cpu_seconds_sum": wall / 4,
        "retry_count": 0,
    }


def experiment_result(experiment, compare_to, status=COMPLETED):
    integrity = {
        "sample_count": 51,
        "sample_order": ["SRR1"],
        "sample_order_sha256": "a" * 64,
        "ordered_contigs": [["NC_1", 100]],
        "variant_count": 1234,
        "per_contig_variant_counts": {"NC_1": 1234},
        "called_alleles": 99,
        "record_sha256": "b" * 64,
        "variant_gt_sha256": "c" * 64,
        "accounting_sha256": "d" * 64,
        "shared_header_sha256": "e" * 64,
        "vcf_sha256": "f" * 64,
        "vcf_index_sha256": "0" * 64,
        "sample_order_valid": True,
        "contig_order_valid": True,
        "allele_accounting_valid": True,
    }
    result = {
        "configuration": {
            "experiment_id": experiment,
            "compare_to": compare_to,
            "plan": "baseline_per_contig",
            "sample_count": 51,
            "batch_size": 50,
            "interval_count": 1,
            "sample_name_map": False,
            "reblock": False,
            "consolidate": False,
            "genomicsdb_cpus": 8,
            "genomicsdb_memory_bytes": 16 * 1024**3,
            "genotype_memory_bytes": 32 * 1024**3,
            "genotype_java_heap_gib": 26,
            "maximum_concurrent_tasks": 3,
            "only_output_calls_starting_in_intervals": True,
            # A private absolute path of the kind that must never be published.
            "workspace_path": SECRET_PATH + "/E0/intervals/interval_0001/workspace",
        },
        "status": status,
        "wall_seconds_including_validation": 4242.5,
        "tasks": [
            {
                "group": {"id": "interval_0001", "intervals": ["NC_1"], "total_bp": 100},
                "import": metrics(10.0, 2 * 1024**3),
                "genotype": metrics(20.0, 17 * 1024**3),
                "workspace": {"bytes": 5000, "file_count": 7},
                "batching": {
                    "expected": [[1, 50], [2, 1]],
                    "observed": [[1, 50], [2, 1]],
                    "completed": [[1, 2], [2, 2]],
                    "reader_initialization_fell_back_to_serial": False,
                },
            }
        ],
        "import": summary(10.0, 2 * 1024**3),
        "genotype": summary(20.0, 17 * 1024**3),
        "gather": metrics(5.0, 1024**3),
        "gather_index": metrics(2.0, 1024**3),
        "workspace": {"bytes": 5000, "file_count": 7},
        "integrity": integrity,
    }
    if compare_to:
        result["comparison"] = dict.fromkeys(summarize_evidence.COMPARISON_KEYS, True)
        result["record_difference_audit"] = {"identical_records": 1234}
    return result


class EvidenceSummaryTests(unittest.TestCase):
    def build(self, root, statuses):
        (root / "execution_lineage.json").write_text(
            json.dumps({"production_sha": "9" * 40, "genotype_memory_gib": 32})
        )
        (root / "host.tsv").write_text(
            "timestamp\tmem_available_bytes\tswap_bytes\tstorage_free_bytes\n"
            f"t0\t{100 * 1024**3}\t{3 * 1024**3}\t{2 * 10**12}\n"
            f"t1\t{80 * 1024**3}\t{4 * 1024**3}\t{1 * 10**12}\n"
        )
        previous = None
        for experiment, status in statuses.items():
            directory = root / experiment
            directory.mkdir()
            (directory / "experiment_result.private.json").write_text(
                json.dumps(experiment_result(experiment, previous, status))
            )
            previous = experiment
        cohort = root / "cohort.json"
        cohort.write_text(
            json.dumps(
                {
                    "evidence": {
                        "cohort_id": "issue45-51gvcf",
                        "sample_count": 51,
                        "unique_biosamples": 51,
                        "reference_accession": "GCF_016808095.1",
                        "production_sha": "9" * 40,
                        "nextflow_version": "26.04.6",
                        "gatk_version": "4.6.2.0",
                        "gatk_container": "broadinstitute/gatk:4.6.2.0@sha256:1",
                        "sample_ploidy": 2,
                        "lineage_verified": True,
                        "benchmark_ready": True,
                        "samples": [
                            {"bioproject": "PRJNA1138464", "gvcf_bytes": 10, "record_count": 5}
                        ],
                    }
                }
            )
        )
        return cohort

    def summarize(self, statuses):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name) / "issue45-e0-e4-run"
        root.mkdir()
        cohort = self.build(root, statuses)
        output = Path(tmp.name) / "evidence.json"
        with unittest.mock.patch.object(
            sys,
            "argv",
            [
                "summarize_evidence.py",
                "--run-dir",
                str(root),
                "--validated-cohort",
                str(cohort),
                "--output",
                str(output),
            ],
        ):
            summarize_evidence.main()
        return output

    def test_public_evidence_never_carries_a_host_absolute_path(self):
        output = self.summarize({"E0": COMPLETED, "E1": COMPLETED})
        text = output.read_text()
        self.assertNotIn(SECRET_PATH, text)
        self.assertNotIn("/data/", text)
        self.assertNotIn("someone-private", text)

    def test_required_metrics_and_comparison_are_all_present(self):
        evidence = json.loads(self.summarize({"E0": COMPLETED, "E1": COMPLETED}).read_text())
        self.assertEqual(51, evidence["cohort"]["sample_count"])
        self.assertEqual(80.0, evidence["host_during_run"]["mem_available_min_gib"])
        self.assertEqual(1.0, evidence["host_during_run"]["swap_delta_gib"])
        baseline = evidence["experiments"]["E0"]
        for key in (
            "import_summary",
            "genotype_summary",
            "gather_wall_seconds",
            "workspace_bytes",
            "workspace_file_count",
            "output_sample_count",
            "output_variant_count",
            "output_accounting_sha256",
            "genotype_memory_gib",
            "tasks",
        ):
            with self.subTest(key=key):
                self.assertIn(key, baseline)
        self.assertEqual([[1, 50], [2, 1]], baseline["tasks"][0]["batches_observed"])
        self.assertEqual(17.0, baseline["tasks"][0]["genotype_peak_rss_gib"])
        self.assertIsNone(baseline["compare_to"])
        self.assertEqual("E0", evidence["experiments"]["E1"]["compare_to"])
        self.assertTrue(
            all(
                evidence["experiments"]["E1"]["comparison_to_baseline"][key]
                for key in summarize_evidence.COMPARISON_KEYS
            )
        )

    def test_one_shared_allocation_is_recorded_as_uniform_and_not_a_confound(self):
        evidence = json.loads(self.summarize({"E0": COMPLETED, "E1": COMPLETED}).read_text())
        regime = evidence["allocation_regime"]
        self.assertTrue(regime["uniform_across_experiments"])
        self.assertEqual({}, regime["varying_keys"])
        self.assertEqual(32.0, regime["values"]["genotype_memory_gib"])
        self.assertEqual(50, regime["values"]["batch_size"])
        self.assertIn(
            "e0_oom_failure", regime["production_baseline_failure_is_recorded_separately"]
        )

    def test_a_differing_allocation_between_experiments_is_reported_as_varying(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name) / "run"
        root.mkdir()
        cohort = self.build(root, {"E0": COMPLETED, "E1": COMPLETED})
        record = root / "E1/experiment_result.private.json"
        payload = json.loads(record.read_text())
        payload["configuration"]["genotype_memory_bytes"] = 16 * 1024**3
        record.write_text(json.dumps(payload))
        output = Path(tmp.name) / "evidence.json"
        with unittest.mock.patch.object(
            sys,
            "argv",
            [
                "summarize_evidence.py",
                "--run-dir",
                str(root),
                "--validated-cohort",
                str(cohort),
                "--output",
                str(output),
            ],
        ):
            summarize_evidence.main()
        regime = json.loads(output.read_text())["allocation_regime"]
        self.assertFalse(regime["uniform_across_experiments"])
        self.assertEqual([16.0, 32.0], regime["varying_keys"]["genotype_memory_gib"])

    def test_a_failed_experiment_is_reported_without_inventing_results(self):
        evidence = json.loads(self.summarize({"E0": "FAILED"}).read_text())
        failed = evidence["experiments"]["E0"]
        self.assertEqual("FAILED", failed["status"])
        for key in ("output_variant_count", "genotype_summary", "tasks"):
            with self.subTest(key=key):
                self.assertNotIn(key, failed)


if __name__ == "__main__":
    unittest.main()
