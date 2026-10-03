import copy
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from build_sample_name_map import build_map  # noqa: E402
from hash_input_fastqs import build_row  # noqa: E402
from plan_joint_genotyping import (  # noqa: E402
    OBSERVED_FASTQ_BYTES,
    SYNTHETIC_REFERENCE_SHA256,
    bind_input_provenance,
    build_plan,
    make_intervals,
)
from resume_cache_inventory import snapshot, verify  # noqa: E402


def measured_sample(sample_id, read_group_id, size):
    return {
        "sample_id": sample_id,
        "read_group_id": read_group_id,
        "bytes": size,
        "fastq_1_checksum": "sha256:" + "1" * 64,
        "fastq_2_checksum": "sha256:" + "2" * 64,
    }


class JointPlanTests(unittest.TestCase):
    def setUp(self):
        self.args = dict(
            contigs=[("chr1", 41), ("short", 1)],
            samples=[
                measured_sample("zeta", "z1", 10),
                measured_sample("alpha", "a1", 20),
            ],
            reference={
                "fasta": {"sha256": SYNTHETIC_REFERENCE_SHA256},
                "fingerprint": "sha256:" + "b" * 64,
            },
            pipeline_sha="c" * 40,
            ploidy=2,
            window=20,
            task_memory_gib=16,
            concurrency=3,
            memory_budget_gib=110,
            free_bytes=1_000_000,
        )

    def test_interval_partition_and_reference_order(self):
        rows = make_intervals(self.args["contigs"], 20)
        self.assertEqual(
            [(r["contig"], r["start"], r["end"]) for r in rows],
            [("chr1", 1, 20), ("chr1", 21, 40), ("chr1", 41, 41), ("short", 1, 1)],
        )
        positions = [
            p for r in rows if r["contig"] == "chr1" for p in range(r["start"], r["end"] + 1)
        ]
        self.assertEqual(positions, list(range(1, 42)))

    def test_boundary_variant_is_owned_by_start(self):
        rows = make_intervals([("chr1", 41)], 20)
        # An indel beginning at 20 and extending into the next interval belongs only to first.
        self.assertEqual([r["rank"] for r in rows if r["start"] <= 20 <= r["end"]], [0])

    def test_research_plan_exposes_assumptions_and_uncompleted_real_replay(self):
        plan = build_plan(**self.args)
        self.assertEqual(plan["status"], "eligible_for_research_execution")
        self.assertEqual(plan["resources"]["retained_bytes_estimated"], 180)
        self.assertEqual(plan["commercial_validation"], "not_established")
        self.assertEqual(plan["evidence"]["new_production_profile_real_replay"], "not_completed")
        self.assertEqual(plan["canonical_sample_order"], ["alpha", "zeta"])

    def test_insufficient_space_blocks(self):
        self.args["free_bytes"] = 179
        self.assertIn(
            "insufficient_snapshot_free_space", build_plan(**self.args)["blocking_reasons"]
        )

    def test_excessive_concurrency_blocks(self):
        self.args["memory_budget_gib"] = 47
        self.assertIn(
            "genotype_concurrency_exceeds_launch_memory_budget",
            build_plan(**self.args)["blocking_reasons"],
        )

    def test_21_samples_need_input_bound_review(self):
        self.args["samples"] = [measured_sample(f"s{i}", f"r{i}", 10) for i in range(21)]
        blocked = build_plan(**self.args)
        self.assertEqual(blocked["status"], "blocked")
        review = {
            **blocked["review_template"],
            "reviewer": "test reviewer",
            "rationale": "synthetic validation",
        }
        self.assertEqual(
            build_plan(**self.args, review=review)["status"], "eligible_for_research_execution"
        )
        changed = copy.deepcopy(self.args)
        changed["samples"][0]["bytes"] += 1
        self.assertEqual(build_plan(**changed, review=review)["status"], "blocked")
        same_size_changed = copy.deepcopy(self.args)
        same_size_changed["samples"][0]["fastq_1_checksum"] = "sha256:" + "3" * 64
        self.assertEqual(build_plan(**same_size_changed, review=review)["status"], "blocked")

    def test_327_remains_no_go_even_with_review(self):
        self.args["samples"] = [measured_sample(f"s{i}", f"r{i}", 10) for i in range(327)]
        plan = build_plan(**self.args)
        review = {**plan["review_template"], "reviewer": "x", "rationale": "x"}
        self.assertEqual(build_plan(**self.args, review=review)["status"], "blocked")

    def test_unobserved_reference_requires_review_even_for_small_cohort(self):
        self.args["reference"]["fasta"]["sha256"] = "a" * 64
        plan = build_plan(**self.args)
        self.assertIn(
            "unobserved_reference_requires_input_bound_validation_review",
            plan["blocking_reasons"],
        )
        review = {**plan["review_template"], "reviewer": "x", "rationale": "new reference trial"}
        self.assertEqual(
            build_plan(**self.args, review=review)["status"], "eligible_for_research_execution"
        )

    def test_larger_input_requires_review_even_with_available_storage(self):
        self.args["samples"][0]["bytes"] = OBSERVED_FASTQ_BYTES + 1
        self.args["free_bytes"] = 10 * OBSERVED_FASTQ_BYTES
        self.assertIn(
            "input_larger_than_observed_requires_input_bound_validation_review",
            build_plan(**self.args)["blocking_reasons"],
        )

    def test_unknown_reference_size_blocks(self):
        self.args["contigs"] = [("big", 500_000_001)]
        self.args["window"] = 20_000_000
        self.assertIn(
            "reference_larger_than_research_envelope", build_plan(**self.args)["blocking_reasons"]
        )

    def test_duplicate_contigs_and_invalid_window_rejected(self):
        for contigs, window in [
            ([("x", 1), ("x", 1)], 20),
            ([("x", 1)], 0),
            ([("x", 1)], 20_000_001),
        ]:
            with self.subTest(contigs=contigs, window=window), self.assertRaises(ValueError):
                make_intervals(contigs, window)

    def test_tiny_window_is_refused_before_materializing_huge_interval_plan(self):
        with self.assertRaisesRegex(ValueError, "interval count 500000000"):
            make_intervals([("chr1", 500_000_000)], 1)

    def test_size_only_summary_cannot_authorize_inputs(self):
        del self.args["samples"][0]["fastq_1_checksum"]
        with self.assertRaisesRegex(ValueError, "measured FASTQ SHA-256"):
            build_plan(**self.args)

    def test_provenance_pairing_is_by_identity_not_completion_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for sample, group, digit in [("alpha", "a1", "4"), ("zeta", "z1", "5")]:
                path = Path(tmp) / f"{group}.tsv"
                path.write_text(
                    "\t".join(
                        [
                            "00000000",
                            sample,
                            group,
                            "library",
                            "ILLUMINA",
                            "",
                            "read1.gz",
                            "sha256:" + digit * 64,
                            "read2.gz",
                            "sha256:" + "6" * 64,
                        ]
                    )
                    + "\n"
                )
                paths.append(path)
            samples = [
                {k: r[k] for k in ("sample_id", "read_group_id", "bytes")}
                for r in self.args["samples"]
            ]
            bound = bind_input_provenance(samples, paths)
            self.assertEqual(bound[0]["sample_id"], "zeta")
            self.assertEqual(bound[0]["fastq_1_checksum"], "sha256:" + "5" * 64)
            self.assertEqual(bound[1]["fastq_1_checksum"], "sha256:" + "4" * 64)
            self.assertEqual(bind_input_provenance(samples, list(reversed(paths))), bound)
            for wrong in [paths[:1], paths + paths[:1]]:
                with self.subTest(wrong=wrong), self.assertRaises(ValueError):
                    bind_input_provenance(samples, wrong)
            paths[0].write_text(paths[0].read_text().replace("alpha\ta1", "other\ta1"))
            with self.assertRaisesRegex(ValueError, "match each planned read group"):
                bind_input_provenance(samples, paths)

    def test_cli_reuses_measured_provenance_and_refuses_same_size_input_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            read1, read2 = root / "r1.gz", root / "r2.gz"
            read1.write_bytes(b"AAA")
            read2.write_bytes(b"CCC")
            provenance = root / "reads.input_provenance.tsv"

            def measure():
                provenance.write_text(
                    build_row(
                        rank=0,
                        sample_id="sample",
                        read_group_id="group",
                        library_id="library",
                        platform="ILLUMINA",
                        platform_unit="",
                        fastq_1=read1,
                        fastq_2=read2,
                    )
                    + "\n"
                )

            measure()
            (root / "reference.fai").write_text("chr1\t41\n")
            (root / "samples.jsonl").write_text(
                json.dumps(
                    {
                        "sample_id": "sample",
                        "read_group_id": "group",
                        "bytes": 6,
                    }
                )
                + "\n"
            )
            (root / "reference.json").write_text(
                json.dumps(
                    {
                        "fasta": {"sha256": "a" * 64},
                        "fingerprint": "sha256:" + "b" * 64,
                    }
                )
            )
            output = root / "plan.json"
            command = [
                sys.executable,
                str(Path(__file__).resolve().parents[2] / "bin/plan_joint_genotyping.py"),
                "--fai",
                str(root / "reference.fai"),
                "--samples",
                str(root / "samples.jsonl"),
                "--input-provenance",
                str(provenance),
                "--reference-manifest",
                str(root / "reference.json"),
                "--pipeline-sha",
                "c" * 40,
                "--ploidy",
                "2",
                "--output-free-bytes",
                "1000000",
                "--output",
                str(output),
                "--intervals",
                str(root / "intervals.tsv"),
                "--require-eligible",
            ]
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 2)
            template = json.loads(output.read_text())["review_template"]
            review = root / "review.json"
            review.write_text(
                json.dumps({**template, "reviewer": "fixture", "rationale": "unit test"})
            )
            self.assertEqual(
                subprocess.run(command + ["--review", str(review)], capture_output=True).returncode,
                0,
            )
            read1.write_bytes(b"TTT")  # Same sample, name and byte count; different content.
            measure()
            self.assertEqual(
                subprocess.run(command + ["--review", str(review)], capture_output=True).returncode,
                2,
            )
            self.assertNotEqual(
                json.loads(output.read_text())["input_summary_sha256"],
                template["input_summary_sha256"],
            )


class SampleMapTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def entry(self, sample):
        vcf = sample + ".g.vcf.gz"
        with gzip.open(self.root / vcf, "wt") as handle:
            handle.write(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + sample
                + "\n"
            )
        (self.root / (vcf + ".tbi")).write_bytes(b"fake index for header contract only")
        return {"sample_id": sample, "gvcf": vcf, "index": vcf + ".tbi"}

    def test_nonlexical_input_keeps_identity_and_declares_canonical_order(self):
        z, a = self.entry("zeta"), self.entry("alpha")
        text, summary = build_map([z, a], self.root)
        self.assertEqual(
            text.splitlines(),
            ["alpha\talpha.g.vcf.gz\talpha.g.vcf.gz.tbi", "zeta\tzeta.g.vcf.gz\tzeta.g.vcf.gz.tbi"],
        )
        self.assertEqual(summary["sample_order"], ["alpha", "zeta"])

    def test_single_sample(self):
        row = self.entry("one")
        self.assertEqual(build_map([row], self.root)[1]["sample_count"], 1)

    def test_wrong_header_rejected(self):
        row = self.entry("one")
        row["sample_id"] = "other"
        with self.assertRaisesRegex(ValueError, "header"):
            build_map([row], self.root)

    def test_swapped_index_rejected(self):
        one, two = self.entry("one"), self.entry("two")
        one["index"] = two["index"]
        with self.assertRaisesRegex(ValueError, "tuple"):
            build_map([one, two], self.root)

    def test_duplicate_sample_rejected(self):
        one = self.entry("one")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            build_map([one, one], self.root)


class ResumeCacheTests(unittest.TestCase):
    def test_hidden_database_members_are_not_mistaken_for_nextflow_logs(self):
        for relative in [".metadata", ".hidden/fragment"]:
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                task = root / "aa" / "task"
                database = task / "interval.genomicsdb"
                database.mkdir(parents=True)
                (task / ".exitcode").write_text("0")
                (task / ".command.log").write_text("not an output")
                output = database / relative
                output.parent.mkdir(exist_ok=True)
                output.write_text("original")
                inventory = snapshot(root)
                self.assertEqual(len(inventory["files"]), 1)
                self.assertEqual(verify(root, inventory), [])
                (task / ".command.log").write_text("updated log")
                self.assertEqual(verify(root, inventory), [])
                output.write_text("modified")
                self.assertEqual(verify(root, inventory)[0]["reason"], "changed")
                output.unlink()
                self.assertEqual(verify(root, inventory)[0]["reason"], "missing_directory_member")
                output.write_text("original")
                (database / ".additional").write_text("new")
                self.assertEqual(
                    verify(root, inventory)[0]["reason"], "unexpected_directory_member"
                )

    def test_added_database_fragment_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / "aa" / "task"
            database = task / "interval.genomicsdb"
            database.mkdir(parents=True)
            (task / ".exitcode").write_text("0")
            (database / "fragment").write_text("data")
            inventory = snapshot(root)
            (database / "additional_fragment").write_text("new calls")
            self.assertEqual(
                verify(root, inventory),
                [
                    {
                        "path": "aa/task/interval.genomicsdb/additional_fragment",
                        "reason": "unexpected_directory_member",
                    }
                ],
            )

    def test_added_database_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / "aa" / "task"
            database = task / "interval.genomicsdb"
            database.mkdir(parents=True)
            (task / ".exitcode").write_text("0")
            (database / "fragment").write_text("data")
            inventory = snapshot(root)
            (database / "link").symlink_to(database / "fragment")
            self.assertEqual(verify(root, inventory)[0]["reason"], "unexpected_directory_member")

    def test_missing_directory_member_remains_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / "aa" / "task"
            database = task / "interval.genomicsdb"
            database.mkdir(parents=True)
            (task / ".exitcode").write_text("0")
            fragment = database / "fragment"
            fragment.write_text("data")
            inventory = snapshot(root)
            fragment.unlink()
            self.assertTrue(database.is_dir())
            self.assertEqual(verify(root, inventory)[0]["reason"], "missing_directory_member")

    def test_changed_and_missing_outputs_are_distinguished(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / "aa" / "task"
            task.mkdir(parents=True)
            (task / ".exitcode").write_text("0")
            output = task / "cohort.vcf"
            output.write_text("original")
            (task / "input.vcf").symlink_to(output)
            inventory = snapshot(root)
            self.assertEqual(len(inventory["files"]), 1)
            self.assertEqual(verify(root, inventory), [])
            output.write_text("modified")
            self.assertEqual(verify(root, inventory)[0]["reason"], "changed")
            output.unlink()
            self.assertEqual(verify(root, inventory)[0]["reason"], "missing")

    def test_failed_task_is_not_blessed_and_unsafe_paths_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / "aa" / "failed"
            task.mkdir(parents=True)
            (task / ".exitcode").write_text("1")
            (task / "partial.vcf").write_text("partial")
            self.assertEqual(snapshot(root)["files"], [])
            manifest = {
                "schema_version": 1,
                "contract": "stopped_local_cache_v1",
                "files": [{"path": "../outside", "bytes": 0, "sha256": ""}],
            }
            with self.assertRaises(ValueError):
                verify(root, manifest)


if __name__ == "__main__":
    unittest.main()
