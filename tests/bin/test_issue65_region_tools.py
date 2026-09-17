"""Issue #65: region asset, mappability, callable and descriptive stratification tools.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "benchmarks" / "issue65"))

import build_region_assets  # noqa: E402
import callable_regions  # noqa: E402
import mappability  # noqa: E402
import stratify_callset  # noqa: E402
from evidence_model import InvalidEvaluationError, validate_evaluation  # noqa: E402
from regions import MalformedRegionError, RegionSet, Stratum  # noqa: E402

# 0-based layout of chrA (40 bp): 0-9 ACGT..., 10-19 lower case, 20-24 N, 25-36 T x12, 37-39 GCG
CHR_A = "ACGTACGTAC" + "acgtacgtac" + "NNNNN" + "T" * 12 + "GCG"
CHR_B = "GGGGCCCCAT"


def _tmp(test: unittest.TestCase) -> Path:
    tmp = tempfile.TemporaryDirectory()
    test.addCleanup(tmp.cleanup)
    return Path(tmp.name)


def _bed(path: Path) -> list[tuple[str, int, int]]:
    return [
        (c, int(s), int(e))
        for c, s, e in (line.split("\t") for line in path.read_text().splitlines())
    ]


class RegionAssetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.dir = _tmp(self)
        self.fasta = self.dir / "ref.fa"
        self.fasta.write_text(f">chrA description\n{CHR_A[:25]}\n{CHR_A[25:]}\n>chrB\n{CHR_B}\n")

    def _gaps(self, rows: str) -> Path:
        path = self.dir / "gaps.txt.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("#accession.version\tstart\tstop\tgap_type\tlinkage_evidence\n" + rows)
        return path

    def test_sequence_derived_assets(self) -> None:
        out = self.dir / "out"
        manifest = build_region_assets.build(
            self.fasta, out, 10, 10, self._gaps("chrA\t21\t25\tscaffold\tpaired-ends\n")
        )
        self.assertEqual(_bed(out / "n_regions.bed"), [("chrA", 20, 25)])
        self.assertEqual(_bed(out / "repeat_windowmasker.bed"), [("chrA", 10, 20)])
        self.assertEqual(_bed(out / "homopolymer_ge10.bed"), [("chrA", 25, 37)])
        self.assertEqual(
            manifest["ncbi_gap_cross_check"],
            {"ncbi_gaps_listed": 1, "all_inside_fasta_n_runs": True},
        )
        # windows: chrA 0-10 GC .5, 10-20 GC .5, 20-30 5 N of 10 (defined half), 30-40 T*7+GCG, chrB GC .8
        self.assertEqual(_bed(out / "gc_45_100.bed"), [("chrA", 0, 20), ("chrB", 0, 10)])
        self.assertEqual(_bed(out / "gc_00_25.bed"), [("chrA", 20, 30)])
        self.assertEqual(_bed(out / "gc_30_35.bed"), [("chrA", 30, 40)])
        self.assertEqual(
            sum(a["bases"] for k, a in manifest["assets"].items() if k.startswith("gc_")), 50
        )
        again = build_region_assets.build(self.fasta, self.dir / "again", 10, 10, None)
        self.assertEqual(
            {k: v["sha256"] for k, v in again["assets"].items()},
            {k: v["sha256"] for k, v in manifest["assets"].items()},
        )

    def test_gap_that_is_not_an_n_run_fails_closed(self) -> None:
        with self.assertRaisesRegex(SystemExit, "1 of 1 NCBI gaps are not N runs"):
            build_region_assets.build(
                self.fasta, self.dir / "bad", 10, 10, self._gaps("chrA\t20\t25\tscaffold\tx\n")
            )

    def test_existing_output_directory_is_not_overwritten(self) -> None:
        (self.dir / "exists").mkdir()
        with self.assertRaises(FileExistsError):
            build_region_assets.build(self.fasta, self.dir / "exists", 10, 10, None)


class MappabilityTests(unittest.TestCase):
    def test_reads_skip_n_and_summary_partitions_reference(self) -> None:
        tmp = _tmp(self)
        fasta = tmp / "ref.fa"
        fasta.write_text(f">chrA\n{CHR_A}\n>chrB\n{CHR_B}\n")
        stream = io.StringIO()
        mappability.reads(fasta, 10, 5, stream)
        names = [line[1:] for line in stream.getvalue().splitlines()[::4]]
        self.assertEqual(
            names, ["chrA:0", "chrA:5", "chrA:10", "chrA:25", "chrA:30", "chrB:0"]
        )  # 15 and 20 contain N

        fai = tmp / "ref.fa.fai"
        fai.write_text("chrA\t40\t6\t40\t41\nchrB\t10\t53\t10\t11\n")
        sam = io.StringIO(
            "@HD\tVN:1.6\n"
            "chrA:0\t0\tchrA\t1\t60\t10M\t*\t0\t0\tX\tI\n"  # high
            "chrA:10\t0\tchrA\t11\t60\t10M\t*\t0\t0\tX\tI\n"  # high (lower case is not special)
            "chrA:5\t0\tchrA\t9\t60\t10M\t*\t0\t0\tX\tI\n"  # placed elsewhere: low
            "chrA:25\t0\tchrA\t26\t0\t10M\t*\t0\t0\tX\tI\n"  # MAPQ 0: low
            "chrA:25\t256\tchrA\t26\t60\t10M\t*\t0\t0\tX\tI\n"  # secondary: ignored
            "chrA:30\t4\t*\t0\t0\t*\t*\t0\t0\tX\tI\n"  # unmapped: low
            "chrB:0\t16\tchrB\t1\t30\t10M\t*\t0\t0\tX\tI\n"  # high (reverse strand is fine)
        )
        manifest = mappability.summarize(fai, 10, 5, 30, sam, tmp / "out")
        self.assertEqual(
            _bed(tmp / "out" / "high_mappability.bed"),
            [("chrA", 0, 5), ("chrA", 10, 15), ("chrB", 0, 5)],
        )
        self.assertEqual(
            _bed(tmp / "out" / "low_mappability.bed"), [("chrA", 5, 10), ("chrA", 25, 35)]
        )
        self.assertEqual(
            _bed(tmp / "out" / "mappability_undefined.bed"),
            [("chrA", 15, 25), ("chrA", 35, 40), ("chrB", 5, 10)],
        )
        self.assertEqual(manifest["reads_aligned"], 6)
        self.assertEqual(sum(a["bases"] for a in manifest["assets"].values()), 50)


class CallableTests(unittest.TestCase):
    def _depth(self, tmp: Path, rows: list[tuple[int, ...]], start: int = 1) -> Path:
        path = tmp / "depth.tsv.gz"
        with gzip.open(path, "wt") as handle:
            for offset, depths in enumerate(rows):
                handle.write("\t".join(["chrA", str(start + offset), *map(str, depths)]) + "\n")
        return path

    def _build(self, tmp: Path, depth: Path, region: str = "chrA:1-8"):
        return callable_regions.build(
            depth=depth,
            samples=["s1", "s2"],
            region=region,
            min_depth=3,
            max_depth_factor=2.0,
            min_sample_fraction=1.0,
            mapq=20,
            baseq=10,
            output=tmp / "out",
        )

    def test_sample_and_cohort_rules(self) -> None:
        tmp = _tmp(self)
        # medians: s1 = 5, s2 = 4 -> max depth 10 and 8
        rows = [(5, 4), (5, 4), (0, 0), (2, 4), (11, 4), (10, 8), (5, 9), (5, 4)]
        summary = self._build(tmp, self._depth(tmp, rows))
        self.assertEqual(summary["samples"]["s1"]["region_median_depth"], 5)
        self.assertEqual(summary["samples"]["s2"]["max_depth"], 8.0)
        self.assertEqual(
            _bed(tmp / "out" / "cohort_callable.bed"),
            [("chrA", 0, 2), ("chrA", 5, 6), ("chrA", 7, 8)],
        )
        self.assertEqual(
            _bed(tmp / "out" / "cohort_non_callable.bed"), [("chrA", 2, 5), ("chrA", 6, 7)]
        )
        self.assertEqual(summary["samples"]["s1"]["callable_bases"], 5)
        self.assertEqual(
            summary["bases"]["cohort_callable"] + summary["bases"]["cohort_non_callable"], 8
        )
        self.assertEqual(_bed(tmp / "out" / "median_depth_00_03.bed"), [("chrA", 2, 3)])
        self.assertIn("not genotype accuracy", summary["rule"]["not_a_claim"])

    def test_gapped_or_short_depth_stream_fails_closed(self) -> None:
        tmp = _tmp(self)
        with self.assertRaisesRegex(SystemExit, "not contiguous"):
            self._build(tmp, self._depth(tmp, [(5, 5), (5, 5)], start=2))
        with self.assertRaisesRegex(SystemExit, "ended at 2"):
            self._build(tmp, self._depth(tmp, [(5, 5), (5, 5)]))


class StratifyTests(unittest.TestCase):
    CONTIGS = (("chrA", 100),)
    SHA = "0" * 64

    def _strata(self, universe: RegionSet) -> list[Stratum]:
        low = RegionSet("gc_low", self.CONTIGS, {"chrA": [(0, 50)]})
        high = RegionSet("gc_high", self.CONTIGS, {"chrA": [(50, 100)]})
        repeat = RegionSet("repeat", self.CONTIGS, {"chrA": [(8, 29)]})
        return [
            Stratum("gc_low", "partition", "gc", low.intersect(universe, "gc_low")),
            Stratum("gc_high", "partition", "gc", high.intersect(universe, "gc_high")),
            Stratum("repeat", "tag", "repeat", repeat),
        ]

    CALLS = (
        "chrA\t10\tA\tG\tPASS\t0/1\t./.\t0/0\n"
        "chrA\t29\tAC\tA\tPASS\t1/1\t0/0\t0/0\n"  # indel straddling the repeat edge
        "chrA\t60\tC\tT\tPASS\t0|1\t0/0\t.|.\n"
        "chrA\t95\tC\tT\tPASS\t0/1\t0/1\t0/1\n"  # outside the universe
    )
    MASKED = (
        "chrA\t10\tA\tG\tPASS\t0/1\t./.\t./.\n"
        "chrA\t29\tAC\tA\tPASS\t1/1\t./.\t./.\n"
        "chrA\t60\tC\tT\tPASS\t0|1\t0/0\t.|.\n"
        "chrA\t95\tC\tT\tPASS\t./.\t0/1\t0/1\n"
    )

    def _run(self, calls: str, masked: str | None):
        universe = RegionSet("universe", self.CONTIGS, {"chrA": [(0, 90)]})
        return stratify_callset.stratify(
            calls=io.StringIO(calls),
            masked_calls=io.StringIO(masked) if masked is not None else None,
            universe=universe,
            strata=self._strata(universe),
            reference_sha256=self.SHA,
            evaluation_id="t",
            query_dataset={"name": "calls"},
            limitations=["test"],
        )

    def test_counts_masking_and_edges(self) -> None:
        record = self._run(self.CALLS, self.MASKED)
        m = record["metrics"]
        self.assertEqual(
            (record["eligible_denominator"], record["exclusion_reasons"]),
            (3, {"outside_universe": 1}),
        )
        self.assertEqual((m["variant_records"], m["snp_records"], m["indel_records"]), (3, 2, 1))
        self.assertEqual(
            (m["genotype_cells"], m["missing_genotype_cells"], m["masked_genotype_cells"]),
            (9, 2, 3),
        )
        self.assertEqual(
            (m["non_reference_calls"], m["heterozygous_calls"], m["homozygous_alt_calls"]),
            (3, 2, 1),
        )
        self.assertEqual(m["heterozygous_fraction_of_non_reference_calls"], round(2 / 3, 6))
        by = m["by_stratum"]
        self.assertEqual(by["gc_low"]["variant_records"] + by["gc_high"]["variant_records"], 3)
        self.assertEqual(by["repeat"]["variant_records"], 1)
        self.assertEqual(
            {
                s["name"]: s["records_straddling_edge_not_counted"]
                for s in record["strata"]
                if "records_straddling_edge_not_counted" in s
            }["repeat"],
            1,
        )
        self.assertEqual(by["gc_high"]["bases"], 40)
        self.assertEqual(by["gc_high"]["variant_records_per_mb"], round(1 / 40 * 1e6, 3))
        for word in ("precision", "recall", "accuracy"):
            self.assertNotIn(word, json.dumps(m))
        self.assertEqual(
            json.dumps(record, sort_keys=True),
            json.dumps(self._run(self.CALLS, self.MASKED), sort_keys=True),
        )

    def test_filter_column_is_split_into_pass_failed_and_unfiltered(self) -> None:
        calls = (
            "chrA\t10\tA\tG\tPASS\t0/1\t0/0\t0/0\n"
            "chrA\t20\tA\tG\tQD2\t0/1\t0/0\t0/0\n"
            "chrA\t30\tA\tG\t.\t0/1\t0/0\t0/0\n"  # raw callset: no filter applied
            "chrA\t40\tA\tG\tQD2;FS60\t0/1\t0/0\t0/0\n"
        )
        m = self._run(calls, None)["metrics"]
        self.assertEqual(m["pass_records"], 1)
        self.assertEqual(m["failed_filter_records"], 2)
        self.assertEqual(m["unfiltered_records"], 1)
        self.assertEqual(
            m["pass_records"] + m["failed_filter_records"] + m["unfiltered_records"],
            m["variant_records"],
        )
        self.assertNotIn("filtered_records", m)

    def test_each_stratum_is_split_by_variant_type(self) -> None:
        by = self._run(self.CALLS, self.MASKED)["metrics"]["by_stratum"]
        self.assertEqual((by["snp"]["variant_records"], by["indel"]["variant_records"]), (2, 1))
        # gc_low holds the SNP at 10 and the indel at 29; gc_high the SNP at 60.
        self.assertEqual(by["gc_low:snp"]["variant_records"], 1)
        self.assertEqual(by["gc_low:indel"]["variant_records"], 1)
        self.assertEqual(by["gc_high:indel"]["variant_records"], 0)
        for name in ("gc_low", "gc_high", "repeat"):
            self.assertEqual(
                by[f"{name}:snp"]["variant_records"] + by[f"{name}:indel"]["variant_records"],
                by[name]["variant_records"],
                name,
            )
            self.assertEqual(by[f"{name}:snp"]["bases"], by[name]["bases"])
        # A density quoted for SNPs counts SNP records only.
        self.assertEqual(
            by["gc_low:snp"]["variant_records_per_mb"], round(1 / by["gc_low"]["bases"] * 1e6, 3)
        )

    def test_unmasked_only_has_no_mask_metrics(self) -> None:
        self.assertNotIn("masked_genotype_cells", self._run(self.CALLS, None)["metrics"])

    def test_mismatched_or_inconsistent_masked_stream_fails_closed(self) -> None:
        with self.assertRaisesRegex(stratify_callset.MalformedCallsError, "different record"):
            self._run(self.CALLS, self.MASKED.replace("chrA\t60\tC\tT", "chrA\t61\tC\tT"))
        with self.assertRaisesRegex(stratify_callset.MalformedCallsError, "fewer missing"):
            self._run(self.CALLS, self.MASKED.replace("./.\t./.\n", "0/0\t0/0\n", 1))
        with self.assertRaisesRegex(stratify_callset.MalformedCallsError, "more records"):
            self._run(self.CALLS, self.MASKED + "chrA\t99\tC\tT\tPASS\t0/1\t0/1\t0/1\n")
        with self.assertRaisesRegex(stratify_callset.MalformedCallsError, "genotypes, expected 3"):
            self._run(self.CALLS + "chrA\t70\tC\tT\tPASS\t0/1\n", None)

    def test_a_partition_that_does_not_cover_the_universe_is_refused(self) -> None:
        universe = RegionSet("universe", self.CONTIGS, {"chrA": [(0, 90)]})
        strata = self._strata(universe)[:1] + self._strata(universe)[2:]
        with self.assertRaisesRegex(MalformedRegionError, "covers"):
            stratify_callset.stratify(
                calls=io.StringIO(self.CALLS),
                masked_calls=None,
                universe=universe,
                strata=strata,
                reference_sha256=self.SHA,
                evaluation_id="t",
                query_dataset={},
                limitations=["x"],
            )

    def test_a_comparison_class_must_name_both_datasets(self) -> None:
        record = self._run(self.CALLS, None)
        for evidence_class in ("caller_concordance", "cross_platform_concordance"):
            broken = dict(record, evidence_class=evidence_class, metrics={"shared_variants": 1})
            with self.assertRaisesRegex(InvalidEvaluationError, "must name its comparator_dataset"):
                validate_evaluation(broken)
        # A self-consistency record describes one callset against the reference sequence.
        self_consistency = self._run(self.CALLS, None)
        self_consistency["evidence_class"] = "reference_sample_self_consistency"
        validate_evaluation(self_consistency)
        self.assertIsNone(self_consistency["comparator_dataset"])

    def test_truth_words_cannot_be_added_to_a_descriptive_record(self) -> None:
        record = self._run(self.CALLS, None)
        record["metrics"]["by_stratum"]["repeat"]["precision"] = 0.99
        with self.assertRaisesRegex(InvalidEvaluationError, "reserved for independent truth"):
            validate_evaluation(record)


if __name__ == "__main__":
    unittest.main()


class DeliveryMatrixTests(unittest.TestCase):
    """The matrix derivation, on evaluation records built from the synthetic fixture."""

    @classmethod
    def setUpClass(cls) -> None:
        import assemble_evidence
        from compare import evaluate
        from regions import read_fai as fai
        from variants import Reference

        cls.assemble = assemble_evidence
        fixture = REPO_ROOT / "tests" / "bin" / "fixtures" / "issue65"
        contigs = fai(fixture / "reference.fa.fai")
        region = RegionSet.from_bed(fixture / "evaluation_region.bed", "evaluation_region", contigs)
        strata = [
            Stratum(
                "cohort_non_callable",
                "tag",
                "callable",
                RegionSet.from_bed(fixture / "low_mappability.bed", "cohort_non_callable", contigs),
            ),
            Stratum(
                "core", "tag", "core", RegionSet.from_bed(fixture / "gc_high.bed", "core", contigs)
            ),
        ]
        reference = Reference(fixture / "reference.fa")
        description = {"sha256": "0" * 64}
        try:

            def run(evaluation_id, cls_):
                return evaluate(
                    evaluation_id=evaluation_id,
                    evidence_class=cls_,
                    reference=reference,
                    reference_description=description,
                    region=region,
                    strata=strata,
                    query={
                        "path": fixture / "query.vcf",
                        "sample": "QUERY",
                        "dataset": {"name": "q"},
                    },
                    comparator={
                        "path": fixture / "truth.vcf",
                        "sample": "TRUTH",
                        "dataset": {"name": "c"},
                    },
                    limitations=["synthetic"],
                )

            cls.evaluations = [
                run("downsampling.X.f050_s65", "downsampling_stability"),
                run("caller.X.bcftools_vs_haplotypecaller", "caller_concordance"),
                *assemble_evidence.not_evaluated_records(description),
            ]
        finally:
            reference.close()
        cls.document = assemble_evidence.matrix(cls.evaluations)
        cls.rows = {row["scope"]: row for row in cls.document["rows"]}

    def test_no_scope_is_marked_unsupported_by_a_benchmark_parameter(self) -> None:
        statuses = {row["overall_status"] for row in self.document["rows"]}
        self.assertNotIn("unsupported", statuses)
        flagged = [row for row in self.document["rows"] if row["fails_benchmark_callable_rule"]]
        self.assertTrue(flagged)
        for row in flagged:
            self.assertTrue(row["caveats"])
            self.assertIn("not calibrated", row["caveats"][0])

    def test_no_scope_is_supported_without_truth(self) -> None:
        self.assertNotIn("supported", {row["overall_status"] for row in self.document["rows"]})
        truth = [c for c in self.document["cells"] if c["evidence_class"] == "independent_truth"]
        self.assertTrue(truth)
        self.assertEqual({c["status"] for c in truth}, {"not_evaluated"})

    def test_statuses_follow_the_rule(self) -> None:
        self.assertEqual(self.rows["window:snp"]["overall_status"], "supported_with_caveat")
        self.assertEqual(
            self.rows["window:snp"]["evidence_classes_evaluated"],
            ["caller_concordance", "downsampling_stability"],
        )
        self.assertEqual(
            self.rows["window:snp"]["key_metrics"]["call_retention.f050_s65"], round(9 / 14, 6)
        )  # SNP units only: 9 shared, 3 absent, 2 no-call
        # The benchmark callable rule is a stratification parameter, not a delivery gate:
        # a scope that fails it is flagged, not downgraded.
        non_callable = self.rows["cohort_non_callable:snp"]
        self.assertEqual(non_callable["overall_status"], "supported_with_caveat")
        self.assertTrue(non_callable["fails_benchmark_callable_rule"])
        self.assertFalse(self.rows["core:snp"]["fails_benchmark_callable_rule"])
        self.assertEqual(
            self.rows["low_mappability:snp"]["overall_status"], "not_evaluated"
        )  # no stratum by that name here
        self.assertEqual(
            self.rows["core:hom_alt"]["overall_status"], "supported_with_caveat"
        )  # 121 ins and 160 in gc_high
        self.assertEqual(
            self.rows["genome_outside_window:snp+indel"]["overall_status"], "not_evaluated"
        )

    def test_a_multi_class_row_gets_a_neutral_combined_claim(self) -> None:
        row = self.rows["window:snp"]
        self.assertGreater(len(row["evidence_classes_evaluated"]), 1)
        claim = row["claim_allowed"]
        self.assertIn("combined caveated evidence", claim)
        self.assertIn("no accuracy claim", claim)
        for evidence_class in row["evidence_classes_evaluated"]:
            self.assertIn(evidence_class, claim)
        # The per-class cells still carry each class's own claim.
        cells = [c for c in self.document["cells"] if c["status"] == "supported_with_caveat"]
        self.assertTrue(cells)
        for cell in cells:
            self.assertNotIn("combined caveated evidence", cell["claim_allowed"])
            self.assertIn("not accuracy", cell["claim_allowed"])

    def test_matrix_is_deterministic_and_tabulates(self) -> None:
        again = self.assemble.matrix(self.evaluations)
        self.assertEqual(
            json.dumps(again, sort_keys=True), json.dumps(self.document, sort_keys=True)
        )
        tsv = self.assemble.matrix_tsv(self.document).splitlines()
        self.assertEqual(len(tsv), len(self.document["rows"]) + 1)
        self.assertTrue(tsv[0].startswith("scope\tsplit\toverall_status"))
        self.assertNotIn("precision", "\n".join(tsv))

    def test_not_evaluated_records_are_valid_and_carry_reasons(self) -> None:
        records = [e for e in self.evaluations if e["not_evaluated_reason"]]
        self.assertEqual(
            {r["evidence_class"] for r in records},
            {
                "independent_truth",
                "technical_replicate_concordance",
                "cross_platform_concordance",
                "reference_sample_self_consistency",
                "descriptive_stratification",
                "downsampling_stability",
            },
        )
        for record in records:
            validate_evaluation(record)
            self.assertEqual(record["metrics"], {})


class CommittedEvidenceTests(unittest.TestCase):
    """The published Issue #65 records obey the contract and match their manifest."""

    EVIDENCE = REPO_ROOT / "docs" / "evidence" / "issue65"

    def test_records_validate_and_match_manifest_hashes(self) -> None:
        import hashlib

        import evidence_model as model

        manifest = json.loads((self.EVIDENCE / "evidence_manifest.json").read_text())
        for name, digest in manifest["outputs_sha256"].items():
            self.assertEqual(
                hashlib.sha256((self.EVIDENCE / name).read_bytes()).hexdigest(), digest, name
            )
        evaluations = json.loads((self.EVIDENCE / "evaluations.json").read_text())
        for record in evaluations:
            model.validate_evaluation(record)
        classes = {r["evidence_class"] for r in evaluations if r["not_evaluated_reason"] is None}
        self.assertNotIn("independent_truth", classes)
        self.assertNotIn("technical_replicate_concordance", classes)
        text = (self.EVIDENCE / "evaluations.json").read_text()
        for word in ('"precision"', '"recall"', '"true_positive"', '"f1"'):
            self.assertNotIn(word, text)

    def test_matrix_claims_nothing_stronger_than_its_evidence(self) -> None:
        import evidence_model as model

        document = json.loads((self.EVIDENCE / "delivery_support_matrix.json").read_text())
        self.assertNotIn("supported", {row["overall_status"] for row in document["rows"]})
        for cell in document["cells"]:
            again = model.delivery_row(
                scope=cell["scope"],
                evidence_class=cell["evidence_class"],
                status=cell["status"],
                evidence_refs=cell["evidence"],
                rationale=cell["rationale"],
            )
            self.assertEqual(again, cell)
