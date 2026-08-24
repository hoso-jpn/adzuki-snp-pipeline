"""Regression guard for Issue #30's dedicated resource labels.

Not a `bin/` CLI script test -- `nextflow.config`'s `withLabel` blocks
are plain text (Groovy DSL) with no CLI surface to exercise, and the
Python scripts these three processes run take no `--memory`-style flag
reflecting `task.memory`, so there is no rendered command line (unlike
GATK's `-Xmx`) to assert against. What can regress silently, without
any nf-test or unit test noticing, is a module file's `label` directive
being reverted back to `process_low` (or `nextflow.config`'s dedicated
block being deleted) -- neither would produce a syntax error, and the
synthetic fixture's tiny inputs would keep passing under any of these
labels' memory value regardless. These tests read the real source
files directly and would fail if either regression happened.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# (module file, expected label)
MODULE_LABELS = (
    ("modules/local/classify_normalized_variants.nf", "process_variant_classification"),
    ("modules/local/summarize_filter_qc.nf", "process_variant_qc_summary"),
    ("modules/local/build_gs_panel.nf", "process_gs_panel"),
    ("modules/local/gatk_haplotypecaller.nf", "process_haplotypecaller"),
)

# (label, expected cpus, expected memory closure text, expected time)
LABEL_CONTRACTS = (
    ("process_variant_classification", "2", "{ 12.GB * task.attempt }", "'2h'"),
    ("process_variant_qc_summary", "2", "{ 12.GB * task.attempt }", "'2h'"),
    ("process_gs_panel", "2", "{ 8.GB * task.attempt }", "'2h'"),
    # process_high itself is unchanged (still cpus=8) -- Issue #30 only
    # benchmarked GATK_HAPLOTYPECALLER, not GATK_GENOTYPEGVCFS (this
    # label's other, unbenchmarked consumer), so the cpus value actually
    # measured lives on process_haplotypecaller instead of here.
    ("process_high", "8", "{ 16.GB * task.attempt }", "'24h'"),
    # Issue #30's real-data GATK_HAPLOTYPECALLER 4-vs-8-cpu benchmark
    # (same real BAM, isolated): ~3.5% wall-time difference, identical
    # gVCF variant records, but double the theoretical concurrency at
    # 4 cpus (floor(32/4)=8 vs floor(32/8)=4 on this host) -- see
    # nextflow.config's own comment for the full numbers.
    ("process_haplotypecaller", "4", "{ 16.GB * task.attempt }", "'24h'"),
)

TEST_PROFILE_LABELS = (
    "process_variant_classification",
    "process_variant_qc_summary",
    "process_gs_panel",
    "process_haplotypecaller",
)


def _extract_label_block(config_text: str, label: str) -> str:
    """Return the `withLabel: <label> { ... }` block's raw text, braces included."""
    marker = f"withLabel: {label} {{"
    start = config_text.index(marker)
    depth = 0
    for index in range(start, len(config_text)):
        if config_text[index] == "{":
            depth += 1
        elif config_text[index] == "}":
            depth -= 1
            if depth == 0:
                return config_text[start : index + 1]
    raise AssertionError(f"unbalanced braces in withLabel: {label} block")


class ModuleLabelTests(unittest.TestCase):
    """Each process must still declare the dedicated label Issue #30 assigned it."""

    def test_each_module_declares_its_dedicated_label(self) -> None:
        for relative_path, expected_label in MODULE_LABELS:
            with self.subTest(module=relative_path):
                text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
                self.assertIn(
                    f"label '{expected_label}'",
                    text,
                    f"{relative_path} no longer declares label '{expected_label}'",
                )
                self.assertNotIn(
                    "label 'process_low'",
                    text,
                    f"{relative_path} was reverted back to label 'process_low'",
                )

    def test_haplotypecaller_does_not_share_process_high(self) -> None:
        # The exact mistake this guard exists for: GATK_HAPLOTYPECALLER
        # silently sharing process_high with GATK_GENOTYPEGVCFS again
        # would change GATK_GENOTYPEGVCFS's own scheduler concurrency as
        # an unreviewed side effect of a future HaplotypeCaller-only
        # change -- the exact coupling this Issue's label split exists
        # to prevent (PR #31 review feedback).
        text = (REPO_ROOT / "modules/local/gatk_haplotypecaller.nf").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("label 'process_high'", text)


class NextflowConfigResourceContractTests(unittest.TestCase):
    """nextflow.config's production values for the three new labels."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.config_text = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")

    def test_each_label_has_expected_cpus_memory_time(self) -> None:
        for label, expected_cpus, expected_memory, expected_time in LABEL_CONTRACTS:
            with self.subTest(label=label):
                block = _extract_label_block(self.config_text, label)
                self.assertRegex(block, rf"cpus\s*=\s*{re.escape(expected_cpus)}\b")
                self.assertIn(expected_memory, block)
                self.assertRegex(block, rf"time\s*=\s*{re.escape(expected_time)}")

    def test_memory_scales_with_task_attempt_not_a_fixed_value(self) -> None:
        # Issue #8/#11's own OOM-retry contract (errorStrategy retries on
        # 137/140/143, maxRetries = 1) only helps if a retried attempt
        # actually gets more memory than the one that was killed -- a
        # fixed (non-attempt-scaled) value would silently defeat that.
        for label, *_rest in LABEL_CONTRACTS:
            with self.subTest(label=label):
                block = _extract_label_block(self.config_text, label)
                self.assertIn("task.attempt", block)


class TestProfileHasNoUnboundedResourceRequest(unittest.TestCase):
    """conf/test.config must override all three new labels to tiny CI-sized values.

    Without this, CI would request nextflow.config's real 8-12 GiB
    production values for every push/PR, rather than a value sized for
    the synthetic fixture.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.config_text = (REPO_ROOT / "conf" / "test.config").read_text(encoding="utf-8")

    def test_each_label_is_overridden_in_test_profile(self) -> None:
        for label in TEST_PROFILE_LABELS:
            with self.subTest(label=label):
                self.assertIn(f"withLabel: {label}", self.config_text)


if __name__ == "__main__":
    unittest.main()
