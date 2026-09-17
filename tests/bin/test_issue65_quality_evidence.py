"""Issue #65: evidence-class contract, region strata and truth/concordance accounting.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "benchmarks" / "issue65"))

import compare  # noqa: E402
import evidence_model as model  # noqa: E402
from regions import (  # noqa: E402
    MalformedRegionError,
    RegionSet,
    Stratum,
    read_fai,
    validate_partition,
)
from variants import MalformedVcfError, Reference, normalize_allele, read_vcf  # noqa: E402

FIXTURE = REPO_ROOT / "tests" / "bin" / "fixtures" / "issue65"
REFERENCE_DESCRIPTION = {"accession": "synthetic", "sha256": "0" * 64}


def _contigs():
    return read_fai(FIXTURE / "reference.fa.fai")


def _region(name: str) -> RegionSet:
    return RegionSet.from_bed(FIXTURE / f"{name}.bed", name, _contigs())


def _strata() -> list[Stratum]:
    return [
        Stratum("repeat", "tag", "repeat", _region("repeat")),
        Stratum("low_mappability", "tag", "mappability", _region("low_mappability")),
        Stratum("gc_low", "partition", "gc", _region("gc_low")),
        Stratum("gc_high", "partition", "gc", _region("gc_high")),
    ]


def _evaluate(evidence_class: str = "synthetic_truth_fixture", **overrides):
    reference = Reference(FIXTURE / "reference.fa")
    try:
        arguments = dict(
            evaluation_id="synthetic",
            evidence_class=evidence_class,
            reference=reference,
            reference_description=REFERENCE_DESCRIPTION,
            region=_region("evaluation_region"),
            strata=_strata(),
            query={"path": FIXTURE / "query.vcf", "sample": "QUERY", "dataset": {"name": "query"}},
            comparator={
                "path": FIXTURE / "truth.vcf",
                "sample": "TRUTH",
                "dataset": {"name": "truth"},
            },
            limitations=["synthetic fixture"],
        )
        arguments.update(overrides)
        return compare.evaluate(**arguments)
    finally:
        reference.close()


class SyntheticTruthAccountingTests(unittest.TestCase):
    """Expected numbers come from the comments in make_synthetic_fixture.py, counted by hand."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.record = _evaluate()
        cls.metrics = cls.record["metrics"]

    def test_totals(self) -> None:
        # TP: 10, 20, 99, 121, 140 C, 140 T, 160, chrU:10, 200, 210, 215.
        m = self.metrics
        self.assertEqual(m["true_positive"], 11)
        self.assertEqual(m["false_positive"], 2)  # 30, 70
        self.assertEqual(m["false_negative"], 3)  # 40, 80, 205 (query calls only `*` there)
        self.assertEqual(m["query_nocall_at_truth_variant"], 2)  # 50, 90
        self.assertEqual(m["truth_nocall_at_query_variant"], 1)  # 60
        self.assertEqual(m["genotype_concordant"], 10)
        self.assertEqual(m["genotype_discordant"], 1)  # 20
        self.assertEqual(m["matched_after_normalization"], 1)  # 99 vs 103

    def test_rates_use_the_stated_denominators(self) -> None:
        m = self.metrics
        precision, recall = 11 / 13, 11 / 16
        self.assertEqual(m["precision"], round(precision, 6))  # TP / (TP + FP)
        self.assertEqual(m["recall"], round(recall, 6))  # TP / (TP + FN + query no-call)
        self.assertEqual(
            m["f1"],
            round(
                2
                * round(precision, 6)
                * round(recall, 6)
                / (round(precision, 6) + round(recall, 6)),
                6,
            ),
        )
        self.assertEqual(m["genotype_concordance"], round(10 / 11, 6))

    def test_denominators_and_exclusions_add_up(self) -> None:
        r = self.record
        self.assertEqual(r["eligible_denominator"], 19)
        self.assertEqual(
            r["exclusion_reasons"],
            {
                "coordinate_mismatch": 1,  # chrZ:5, written by both sides
                "outside_region": 1,  # 250
                "ref_mismatch": 2,  # 180 (query only) and 220 (both sides)
                "region_boundary": 1,  # 238
                "symbolic_allele": 4,  # 190 <DEL>, 205 `*`, 210 `*`, 215 `*`
            },
        )
        self.assertEqual(r["excluded_denominator"], 9)
        m = self.metrics
        self.assertEqual(
            r["eligible_denominator"],
            m["true_positive"]
            + m["false_positive"]
            + m["false_negative"]
            + m["query_nocall_at_truth_variant"]
            + m["truth_nocall_at_query_variant"],
        )

    def test_strata_tags_overlap_and_partitions_add_up(self) -> None:
        by = self.metrics["by_stratum"]
        self.assertEqual(by["repeat"]["true_positive"], 1)  # the normalized homopolymer deletion
        self.assertEqual(
            (
                by["low_mappability"]["true_positive"],
                by["low_mappability"]["false_positive"],
                by["low_mappability"]["false_negative"],
            ),
            (2, 1, 1),
        )
        for metric in (
            "true_positive",
            "false_positive",
            "false_negative",
            "query_nocall_at_truth_variant",
        ):
            self.assertEqual(
                by["gc_low"][metric] + by["gc_high"][metric], self.metrics[metric], metric
            )
            self.assertEqual(by["snp"][metric] + by["indel"][metric], self.metrics[metric], metric)
        self.assertEqual(by["indel"]["true_positive"], 2)
        # Dosage comes from the truth call, or the query call where truth is not positive:
        # TP 1/1 at 20 (truth dosage 2), 121 and 160; FP at 30 and 70 are query 0/1.
        self.assertEqual((by["hom_alt"]["true_positive"], by["het"]["true_positive"]), (3, 8))
        self.assertEqual((by["hom_alt"]["false_positive"], by["het"]["false_positive"]), (0, 2))
        self.assertEqual(by["indel:hom_alt"]["true_positive"], 1)
        self.assertEqual(by["repeat:indel"]["true_positive"], 1)
        self.assertEqual(by["repeat:snp"]["true_positive"], 0)
        for metric in ("true_positive", "false_positive", "false_negative"):
            self.assertEqual(
                by["low_mappability:snp"][metric] + by["low_mappability:indel"][metric],
                by["low_mappability"][metric],
            )
        # Tags are not additive: repeat and low_mappability together cover fewer units than the total.
        self.assertLess(
            by["repeat"]["true_positive"] + by["low_mappability"]["true_positive"],
            self.metrics["true_positive"],
        )

    def test_true_negatives_are_not_counted_and_say_why(self) -> None:
        self.assertNotIn("true_negative", self.metrics)
        self.assertIn("not counted", self.record["definitions"]["true_negatives"])

    def test_output_is_deterministic(self) -> None:
        again = _evaluate()
        self.assertEqual(json.dumps(again, sort_keys=True), json.dumps(self.record, sort_keys=True))
        self.assertEqual(self.record["region_definition_hash"], again["region_definition_hash"])


class ConcordanceVocabularyTests(unittest.TestCase):
    def test_same_accounting_without_accuracy_words(self) -> None:
        record = _evaluate("caller_concordance")
        m = record["metrics"]
        self.assertEqual(m["shared_variants"], 11)
        self.assertEqual(m["only_in_query"], 2)
        self.assertEqual(m["only_in_comparator"], 3)
        for word in ("precision", "recall", "f1", "true_positive", "false_positive"):
            self.assertNotIn(word, json.dumps(m))
        self.assertEqual(m["site_agreement_rate"], round(11 / 19, 6))

    def test_downsampling_names_retention_not_recall(self) -> None:
        m = _evaluate("downsampling_stability")["metrics"]
        self.assertEqual(m["only_in_full_depth"], 3)
        self.assertEqual(m["call_retention"], round(11 / 16, 6))
        self.assertNotIn("recall", m)

    def test_descriptive_class_cannot_compare(self) -> None:
        with self.assertRaises(ValueError):
            _evaluate("descriptive_stratification")


class MultiAllelicSymbolicTests(unittest.TestCase):
    """A symbolic ALT excludes itself only; an excluded unit is counted once."""

    HEADER = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    )

    def _vcf(self, sample: str, rows: list[tuple]) -> Path:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / f"{sample}.vcf"
        body = "".join(
            "\t".join([contig, str(pos), ".", ref, alt, "50", "PASS", ".", "GT", gt]) + "\n"
            for contig, pos, ref, alt, gt in rows
        )
        path.write_text(self.HEADER.format(sample=sample) + body)
        return path

    def _run(self, query_rows: list[tuple], comparator_rows: list[tuple]):
        return _evaluate(
            query={
                "path": self._vcf("QUERY", query_rows),
                "sample": "QUERY",
                "dataset": {"n": "q"},
            },
            comparator={
                "path": self._vcf("TRUTH", comparator_rows),
                "sample": "TRUTH",
                "dataset": {"n": "c"},
            },
        )

    # chrT:10 is C in the fixture reference and chrT:20 is C as well.
    def test_sequence_allele_survives_a_symbolic_partner(self) -> None:
        for name, query, comparator in (
            ("query only", [("chrT", 10, "C", "G,*", "0/1")], [("chrT", 10, "C", "G", "0/1")]),
            ("comparator only", [("chrT", 10, "C", "G", "0/1")], [("chrT", 10, "C", "G,*", "0/1")]),
            ("both", [("chrT", 10, "C", "G,*", "0/1")], [("chrT", 10, "C", "G,*", "0/1")]),
        ):
            with self.subTest(name):
                record = self._run(query, comparator)
                self.assertEqual(record["metrics"]["true_positive"], 1)
                self.assertEqual(record["eligible_denominator"], 1)
                # `*` is not called by GT 0/1, so it is not a comparison unit at all.
                self.assertEqual(record["exclusion_reasons"], {})

    def test_a_called_symbolic_allele_is_the_only_excluded_unit(self) -> None:
        record = self._run(
            [("chrT", 10, "C", "G,*", "0/2")],  # query calls only the symbolic allele
            [("chrT", 10, "C", "G", "0/1")],
        )
        m = record["metrics"]
        self.assertEqual((m["true_positive"], m["false_negative"]), (0, 1))
        self.assertEqual(record["exclusion_reasons"], {"symbolic_allele": 1})
        self.assertEqual(
            record["definitions"]["excluded_units_by_side"]["symbolic_allele"],
            {"query_only": 1, "comparator_only": 0, "both": 0},
        )

    def test_gt_one_two_counts_the_sequence_allele_and_excludes_the_symbolic_one(self) -> None:
        record = self._run(
            [("chrT", 10, "C", "G,*", "1/2")],
            [("chrT", 10, "C", "G,*", "1/2")],
        )
        m = record["metrics"]
        self.assertEqual((m["true_positive"], m["genotype_concordant"]), (1, 1))
        # Both sides call the same `*`: one excluded unit, not two.
        self.assertEqual(record["exclusion_reasons"], {"symbolic_allele": 1})
        self.assertEqual(
            record["definitions"]["excluded_units_by_side"]["symbolic_allele"]["both"], 1
        )

    def test_the_same_invalid_unit_on_both_sides_is_counted_once(self) -> None:
        rows = [
            ("chrZ", 5, "A", "C", "0/1"),
            ("chrT", 20, "A", "G", "0/1"),
        ]  # chrT:20 REF is C, so "A" is a mismatch
        record = self._run(rows, rows)
        self.assertEqual(record["exclusion_reasons"], {"coordinate_mismatch": 1, "ref_mismatch": 1})
        self.assertEqual(record["excluded_denominator"], 2)
        for reason in ("coordinate_mismatch", "ref_mismatch"):
            self.assertEqual(
                record["definitions"]["excluded_units_by_side"][reason]["both"], 1, reason
            )

    def test_reason_precedence_is_declared_and_applied(self) -> None:
        self.assertEqual(
            compare.EXCLUSION_PRECEDENCE,
            (
                "coordinate_mismatch",
                "ref_mismatch",
                "symbolic_allele",
                "outside_region",
                "region_boundary",
            ),
        )
        # A symbolic ALT on a record whose site is already unusable takes the site's reason.
        record = self._run(
            [("chrZ", 5, "A", "C,*", "1/2"), ("chrT", 20, "A", "G,*", "1/2")],
            [],
        )
        self.assertEqual(record["exclusion_reasons"], {"coordinate_mismatch": 2, "ref_mismatch": 2})
        self.assertEqual(record["definitions"]["exclusion_precedence"][0], "coordinate_mismatch")

    def test_exclusion_reasons_always_sum_to_the_denominator(self) -> None:
        record = self._run(
            [
                ("chrT", 10, "C", "G,*", "1/2"),
                ("chrT", 250, "A", "C", "0/1"),
                ("chrZ", 5, "A", "C", "0/1"),
            ],
            [("chrT", 10, "C", "G", "0/1"), ("chrT", 250, "A", "C", "0/1")],
        )
        self.assertEqual(sum(record["exclusion_reasons"].values()), record["excluded_denominator"])
        self.assertEqual(
            record["exclusion_reasons"],
            {"coordinate_mismatch": 1, "outside_region": 1, "symbolic_allele": 1},
        )


class NormalizationTests(unittest.TestCase):
    def test_python_normalization_matches_pinned_bcftools_norm(self) -> None:
        reference = Reference(FIXTURE / "reference.fa")
        try:
            ours = sorted(
                normalize_allele(reference, r.contig, r.pos, r.ref, alt, i).key
                for _h, r in read_vcf(FIXTURE / "normalization_input.vcf")
                for i, alt in enumerate(r.alts, start=1)
            )
        finally:
            reference.close()
        # Written by `bcftools norm --no-version -f reference.fa -m -any` (bcftools 1.24, pinned container).
        golden = sorted(
            (r.contig, r.pos, r.ref, alt)
            for _h, r in read_vcf(FIXTURE / "normalization_bcftools_norm.vcf")
            for alt in r.alts
        )
        self.assertEqual(ours, golden)
        self.assertIn(("chrT", 99, "CA", "C"), ours)
        self.assertIn(("chrT", 18, "A", "AC"), ours)


class RegionTests(unittest.TestCase):
    def test_boundary_inside_outside(self) -> None:
        region = _region("evaluation_region")
        self.assertEqual(region.contains_span("chrT", 240, 1), "inside")  # 0-based 239
        self.assertEqual(region.contains_span("chrT", 241, 1), "outside")
        self.assertEqual(region.contains_span("chrT", 238, 4), "boundary")
        self.assertEqual(region.contains_span("chrU", 1, 1), "inside")
        self.assertEqual(region.contains_span("chrZ", 1, 1), "outside")

    def test_set_algebra_and_complement(self) -> None:
        contigs = _contigs()
        a = RegionSet("a", contigs, {"chrT": [(0, 10), (5, 20), (30, 40)]})
        b = RegionSet("b", contigs, {"chrT": [(15, 35)]})
        self.assertEqual(a.intervals["chrT"], [(0, 20), (30, 40)])
        self.assertEqual(a.intersect(b, "i").intervals["chrT"], [(15, 20), (30, 35)])
        self.assertEqual(a.subtract(b, "s").intervals["chrT"], [(0, 15), (35, 40)])
        self.assertEqual(a.complement("c").bases(), 300 + 100 - 30)

    def test_region_hash_tracks_membership(self) -> None:
        contigs = _contigs()
        first = RegionSet("r", contigs, {"chrT": [(0, 10)]})
        same = RegionSet("r", contigs, {"chrT": [(0, 5), (5, 10)]})
        other = RegionSet("r", contigs, {"chrT": [(0, 11)]})
        sha = "0" * 64
        self.assertEqual(first.definition_hash(sha), same.definition_hash(sha))
        self.assertNotEqual(first.definition_hash(sha), other.definition_hash(sha))
        self.assertNotEqual(first.definition_hash(sha), first.definition_hash("1" * 64))

    def test_partition_must_be_disjoint_and_complete(self) -> None:
        universe = _region("evaluation_region")
        validate_partition([_region("gc_low"), _region("gc_high")], universe)
        with self.assertRaisesRegex(MalformedRegionError, "overlap"):
            validate_partition([_region("gc_low"), _region("low_mappability")], universe)
        with self.assertRaisesRegex(MalformedRegionError, "covers"):
            validate_partition([_region("gc_low")], universe)

    def test_malformed_bed_fails_closed(self) -> None:
        contigs = _contigs()
        for text, expected in (
            ("chrT\t10\n", "not a BED interval"),
            ("chrT\t-1\t5\n", "not a BED interval"),
            ("chrT\t20\t10\n", "outside"),
            ("chrT\t0\t301\n", "outside"),
            ("chrZ\t0\t5\n", "not in the reference"),
        ):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "x.bed"
                path.write_text(text)
                with self.assertRaisesRegex(MalformedRegionError, expected):
                    RegionSet.from_bed(path, "x", contigs)


class MalformedVcfTests(unittest.TestCase):
    def _write(self, body: str) -> Path:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "x.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            + body
        )
        return path

    def test_malformed_rows_fail_closed(self) -> None:
        for body, expected in (
            ("chrT\t10\t.\tC\tG\t50\tPASS\t.\tGT\n", "columns"),
            ("chrT\tten\t.\tC\tG\t50\tPASS\t.\tGT\t0/1\n", "not an integer"),
            ("chrT\t10\t.\tC\tG\t50\tPASS\t.\tDP\t5\n", "no GT"),
            ("chrT\t10\t.\tC\tG\t50\tPASS\t.\tGT\tx/1\n", "not a genotype"),
        ):
            with self.subTest(body=body), self.assertRaisesRegex(MalformedVcfError, expected):
                list(read_vcf(self._write(body)))

    def test_a_missing_sample_fails_closed(self) -> None:
        with self.assertRaisesRegex(MalformedVcfError, "not in the VCF"):
            list(read_vcf(self._write(""), ["OTHER"]))


class EvidenceContractTests(unittest.TestCase):
    def _record(self, **overrides):
        record = _evaluate("caller_concordance")
        record.update(overrides)
        return record

    def test_accuracy_words_are_refused_outside_truth(self) -> None:
        record = self._record()
        record["metrics"]["precision"] = 0.9
        with self.assertRaisesRegex(model.InvalidEvaluationError, "reserved for independent truth"):
            model.validate_evaluation(record)
        stratum = self._record()
        stratum["metrics"]["by_stratum"]["snp"]["recall"] = 0.5
        with self.assertRaisesRegex(model.InvalidEvaluationError, "stratum 'snp'"):
            model.validate_evaluation(stratum)

    def test_denominator_accounting_is_enforced(self) -> None:
        record = self._record(excluded_denominator=4)
        with self.assertRaisesRegex(model.InvalidEvaluationError, "sum to excluded_denominator"):
            model.validate_evaluation(record)

    def test_not_evaluated_is_first_class(self) -> None:
        record = compare.not_evaluated(
            evaluation_id="truth",
            evidence_class="independent_truth",
            reference_description=REFERENCE_DESCRIPTION,
            region_set="all",
            reason="no independent truth set exists for this reference",
            limitations=["none found in the audit"],
        )
        self.assertEqual(record["metrics"], {})
        with self.assertRaisesRegex(model.InvalidEvaluationError, "cannot carry metrics"):
            model.validate_evaluation(record | {"metrics": {"precision": 1.0}})
        with self.assertRaises(model.InvalidEvaluationError):
            model.validate_evaluation(record | {"not_evaluated_reason": "  "})

    def test_unknown_class_and_missing_fields_fail(self) -> None:
        with self.assertRaisesRegex(model.InvalidEvaluationError, "unknown evidence_class"):
            model.validate_evaluation(self._record(evidence_class="vibes"))
        record = self._record()
        del record["limitations"]
        with self.assertRaisesRegex(model.InvalidEvaluationError, "lacks"):
            model.validate_evaluation(record)

    def test_delivery_status_is_bounded_by_evidence(self) -> None:
        row = model.delivery_row(
            scope="callable SNP",
            evidence_class="caller_concordance",
            status="supported_with_caveat",
            evidence_refs=["e1"],
            rationale="r",
        )
        self.assertIn("not accuracy", row["claim_allowed"])
        with self.assertRaisesRegex(model.InvalidEvaluationError, "at most supported_with_caveat"):
            model.delivery_row(
                scope="s",
                evidence_class="downsampling_stability",
                status="supported",
                evidence_refs=["e"],
                rationale="r",
            )
        with self.assertRaisesRegex(model.InvalidEvaluationError, "only validate the harness"):
            model.delivery_row(
                scope="s",
                evidence_class="synthetic_truth_fixture",
                status="supported",
                evidence_refs=["e"],
                rationale="r",
            )
        combined = model.delivery_row(
            scope="s",
            evidence_class=["caller_concordance", "downsampling_stability"],
            status="supported_with_caveat",
            evidence_refs=["e1", "e2"],
            rationale="r",
        )
        self.assertIn("combined caveated evidence", combined["claim_allowed"])
        self.assertIn("no accuracy claim", combined["claim_allowed"])
        with self.assertRaisesRegex(model.InvalidEvaluationError, "at most supported_with_caveat"):
            model.delivery_row(
                scope="s",
                evidence_class=["independent_truth", "caller_concordance"],
                status="supported",
                evidence_refs=["e"],
                rationale="r",
            )
        with self.assertRaisesRegex(model.InvalidEvaluationError, "cites no evidence"):
            model.delivery_row(
                scope="s",
                evidence_class=None,
                status="not_evaluated",
                evidence_refs=["e"],
                rationale="r",
            )
        with self.assertRaisesRegex(model.InvalidEvaluationError, "must cite evidence"):
            model.delivery_row(
                scope="s",
                evidence_class="independent_truth",
                status="supported",
                evidence_refs=[],
                rationale="r",
            )
        self.assertEqual(
            model.delivery_row(
                scope="s",
                evidence_class="independent_truth",
                status="supported",
                evidence_refs=["t"],
                rationale="r",
            )["status"],
            "supported",
        )


if __name__ == "__main__":
    unittest.main()
