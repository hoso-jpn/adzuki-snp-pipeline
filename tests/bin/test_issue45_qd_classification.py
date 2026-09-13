"""The QD difference gate must not excuse a difference for being numerically close."""

import sys
import unittest
from pathlib import Path

HELPERS = Path(__file__).resolve().parents[2] / "benchmarks/issue45"
sys.path.insert(0, str(HELPERS))
import classify_qd_differences as qd  # noqa: E402


def record(qual="525.07", qd_value="30.00", filt=".", gt="0/1", ad="10,10", dp="20"):
    return [
        "NC_1",
        "1000",
        ".",
        "T",
        "C",
        qual,
        filt,
        f"AC=1;AN=2;QD={qd_value}",
        "GT:AD:DP",
        f"{gt}:{ad}:{dp}",
    ]


class ClassificationTests(unittest.TestCase):
    def test_above_the_helper_margin_is_directly_explained(self):
        self.assertEqual(
            ("directly_explained_high_qd_jitter", "directly_explained"),
            qd.classify(record(qd_value="28.0"), record(qd_value="31.0"), 40.0),
        )

    def test_inside_the_margin_but_above_gatk_threshold_is_reconstruction_consistent(self):
        self.assertEqual(
            ("reconstruction_margin_above_gatk_threshold", "reconstruction_limit_consistent"),
            qd.classify(record(qd_value="28.0"), record(qd_value="31.0"), 35.004),
        )

    def test_below_gatk_threshold_stays_unresolved_however_close(self):
        for raw in (34.999, 34.5, 30.0, 2.0):
            with self.subTest(raw=raw):
                category, rollup = qd.classify(
                    record(qd_value="28.0"), record(qd_value="31.0"), raw
                )
                self.assertEqual("unresolved", rollup)
                self.assertEqual("near_threshold_below_gatk_threshold", category)

    def test_exactly_at_the_threshold_errs_toward_unresolved(self):
        _, rollup = qd.classify(record(qd_value="28.0"), record(qd_value="31.0"), 35.0)
        self.assertEqual("unresolved", rollup)

    def test_an_output_matching_the_deterministic_value_is_not_excused(self):
        """If one side equals raw QD, jitter did not fire on both sides."""
        _, rollup = qd.classify(record(qd_value="35.00"), record(qd_value="31.0"), 35.004)
        self.assertEqual("unresolved", rollup)

    def test_a_qd_implausible_as_a_jitter_draw_is_not_excused(self):
        for implausible in ("1.0", "80.0"):
            with self.subTest(implausible=implausible):
                _, rollup = qd.classify(
                    record(qd_value=implausible), record(qd_value="31.0"), 35.004
                )
                self.assertEqual("unresolved", rollup)

    def test_unreconstructable_or_non_finite_is_unresolved(self):
        self.assertEqual(
            ("cannot_reconstruct_raw_qd", "unresolved"),
            qd.classify(record(qd_value="28.0"), record(qd_value="31.0"), None),
        )
        self.assertEqual(
            ("non_finite_or_missing_annotation", "unresolved"),
            qd.classify(record(qd_value="nan"), record(qd_value="31.0"), 40.0),
        )
        missing = record()
        missing[7] = "AC=1;AN=2"
        self.assertEqual(
            ("non_finite_or_missing_annotation", "unresolved"),
            qd.classify(missing, record(qd_value="31.0"), 40.0),
        )

    def test_gatk_constants_match_the_pinned_implementation(self):
        self.assertEqual(35.0, qd.GATK_MAX_QD_BEFORE_FIXING)
        self.assertEqual(30.0, qd.GATK_IDEAL_HIGH_QD)
        self.assertEqual(3.0, qd.GATK_HIGH_QD_SD)
        self.assertGreater(qd.HELPER_MARGIN, qd.GATK_MAX_QD_BEFORE_FIXING)


if __name__ == "__main__":
    unittest.main()
