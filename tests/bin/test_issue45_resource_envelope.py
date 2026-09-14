"""The 327 projection must keep observation, trend, assumption and projection apart."""

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks/issue45"))
import project_resource_envelope as envelope  # noqa: E402

EVIDENCE = ROOT / "docs/evidence/issue45"


class FitTests(unittest.TestCase):
    def test_linear_and_power_fits_recover_exact_relationships(self):
        linear = envelope.linear_fit([1, 2, 4], [3, 5, 9])
        self.assertAlmostEqual(1.0, linear["intercept"])
        self.assertAlmostEqual(2.0, linear["slope"])
        self.assertAlmostEqual(1.0, linear["r2"])
        power = envelope.power_fit([1, 2, 4], [3, 12, 48])
        self.assertAlmostEqual(3.0, power["coefficient"])
        self.assertAlmostEqual(2.0, power["exponent"])

    def test_harmonic_number_is_wattersons_a_n(self):
        self.assertAlmostEqual(1 + 1 / 2 + 1 / 3, envelope.harmonic(4))

    def test_ratio_range_spans_all_three_trend_forms(self):
        result = envelope.ratio_range([13, 26, 51], [10, 20, 60], target=102)
        ratios = result["ratio_to_last_observed"].values()
        self.assertEqual(min(ratios), result["ratio_min"])
        self.assertEqual(max(ratios), result["ratio_max"])
        # The upper pair (26 -> 51) is steeper than the three-point line here.
        self.assertEqual(result["ratio_to_last_observed"]["upper_pair_slope"], result["ratio_max"])

    def test_variant_scenarios_are_ordered_and_start_from_the_observed_floor(self):
        scenarios = envelope.variant_scenarios(1000, 130, 500)
        counts = [scenario["variants"] for scenario in scenarios.values()]
        self.assertEqual(1000, counts[0])
        self.assertEqual(sorted(counts), counts)
        for scenario in scenarios.values():
            self.assertTrue(scenario["assumption"])


class CommittedEnvelopeTests(unittest.TestCase):
    """Re-derive the committed envelope, so the evidence cannot drift from its inputs."""

    @classmethod
    def setUpClass(cls):
        import csv

        with (EVIDENCE / "candidate_cohort.tsv").open() as handle:
            rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row["library_strategy"] == "WGS" and row["library_layout"] == "PAIRED"
            ]
        cls.envelope = envelope.build_envelope(
            json.loads((EVIDENCE / "benchmark_results.json").read_text()),
            json.loads((EVIDENCE / "sample_scaling_20260914.json").read_text()),
            json.loads((EVIDENCE / "host_and_generation_observed.json").read_text()),
            sum(int(size) for row in rows for size in row["fastq_bytes"].split(";")),
        )

    def test_committed_evidence_matches_a_fresh_derivation(self):
        committed = json.loads((EVIDENCE / "resource_envelope_327.json").read_text())
        self.assertEqual(json.loads(json.dumps(self.envelope)), committed)

    def test_the_five_parts_are_separate_sections(self):
        for section in ("observed", "empirical_trends", "assumptions", "projected", "uncertainty"):
            self.assertIn(section, self.envelope)

    def test_memory_projection_is_anchored_at_the_observed_51_sample_window(self):
        observed = self.envelope["observed"]["scaling_levels"][-1]
        floor = self.envelope["projected"]["per_variant_scenario"]["S0_no_new_sites"]
        low, high = self.envelope["empirical_trends"]["genotype_rss_vs_cells"][
            "slope_range_gib_per_million_cells"
        ]
        added = floor["window_genotype_cells_millions"] - observed["genotype_cells_millions"]
        self.assertAlmostEqual(
            observed["genotype_peak_rss_gib"] + low * added,
            floor["genotype_peak_rss_gib_20mb_window"][0],
            delta=0.1,
        )
        self.assertAlmostEqual(
            observed["genotype_peak_rss_gib"] + high * added,
            floor["genotype_peak_rss_gib_20mb_window"][1],
            delta=0.1,
        )

    def test_sample_count_alone_is_not_the_memory_predictor(self):
        trends = self.envelope["empirical_trends"]
        self.assertIn("genotype_rss_vs_n_alone_not_used", trends)
        self.assertIn("composition", trends["genotype_rss_vs_n_alone_not_used"]["reason"])


if __name__ == "__main__":
    unittest.main()
