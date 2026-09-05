"""The paths tooling writes into the workspace must stay git-ignored.

Run with: python3 -m unittest discover -s tests/bin -v

Issue #42 made an unclean working tree a hard failure: a run manifest
records the HEAD commit as the code that produced the run, and that is
only true if the tree matched HEAD. The escape hatch for everything that
legitimately appears in the workspace without being part of the pipeline
-- a run's own outputs, the test framework's scratch, CI installers that
unpack themselves in place -- is `.gitignore`.

That makes these ignore rules load-bearing rather than cosmetic: if one
of them regresses, the provenance guard fires on the pipeline's own
tooling and every full-pipeline test fails with a message about
uncommitted changes. This pins them, cheaply, so that failure shows up
here instead of as a confusing red CI run.
"""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Each entry is a path something other than a human writes into the
# working tree, with the reason it appears.
IGNORED_TOOLING_PATHS: tuple[tuple[str, str], ...] = (
    ("work/ab/cd/.command.sh", "a Nextflow task work directory"),
    ("results/variants/raw/cohort.raw.vcf.gz", "a run's own published output"),
    (".nextflow.log", "Nextflow's own log"),
    (".nextflow/cache/db", "Nextflow's own cache"),
    (".nf-test/tests/abc/output.json", "nf-test's working directory"),
    (
        ".nf-test-0123456789abcdef0123456789abcdef.nf",
        "the mock script nf-test generates per process test, in the project root",
    ),
    (
        "setup-nextflow/action.yml",
        "nf-core/setup-nextflow checking its own action repository out into the "
        "workspace before installing the engine",
    ),
    ("bin/__pycache__/manifest_utils.cpython-312.pyc", "compiled Python bytecode"),
)


def _is_ignored(relative_path: str) -> bool:
    """Ask git itself whether this repository would ignore the path."""
    completed = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "check-ignore", "-q", "--no-index", relative_path],
        capture_output=True,
        check=False,
    )
    # 0 = ignored, 1 = not ignored, anything else = git could not answer.
    assert completed.returncode in (0, 1), completed.stderr.decode("utf-8", "replace")
    return completed.returncode == 0


class ProvenanceGuardIgnoreRuleTests(unittest.TestCase):
    def test_tooling_paths_are_ignored(self) -> None:
        for relative_path, reason in IGNORED_TOOLING_PATHS:
            with self.subTest(path=relative_path):
                self.assertTrue(
                    _is_ignored(relative_path),
                    f"{relative_path} is not git-ignored, so the Issue #42 "
                    f"working-tree guard would treat {reason} as a change to "
                    "the pipeline and fail every full-pipeline run",
                )

    def test_pipeline_source_is_not_ignored(self) -> None:
        # The rules above must not have grown so broad that a real source
        # change stops counting as one -- that would silently defeat the
        # guard instead of the tooling false positives it exists for.
        for relative_path in (
            "main.nf",
            "nextflow.config",
            "bin/build_run_manifest.py",
            "modules/local/build_run_manifest.nf",
            "workflows/adzuki_snp_pipeline.nf",
        ):
            with self.subTest(path=relative_path):
                self.assertFalse(_is_ignored(relative_path))


if __name__ == "__main__":
    unittest.main()
