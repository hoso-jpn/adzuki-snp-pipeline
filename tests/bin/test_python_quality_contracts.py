import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PYTHON_PATHS = ("bin", "tests/bin", "tests/scripts")


class PythonQualityContractsTest(unittest.TestCase):
    def test_ruff_is_exactly_pinned(self) -> None:
        requirements = (ROOT / "requirements-dev.txt").read_text(encoding="utf-8")
        self.assertRegex(requirements, r"(?m)^ruff==\d+\.\d+\.\d+$")

    def test_ruff_targets_python_312(self) -> None:
        config = (ROOT / "ruff.toml").read_text(encoding="utf-8")
        self.assertIn('target-version = "py312"', config)

    def test_ci_runs_lint_and_format_for_every_managed_path(self) -> None:
        workflow = (ROOT / ".github/workflows/lint.yml").read_text(encoding="utf-8")
        self.assertIn('python-version: "3.12"', workflow)
        self.assertIn("pip install --requirement requirements-dev.txt", workflow)

        paths = " ".join(PYTHON_PATHS)
        self.assertIn(f"ruff check {paths}", workflow)
        self.assertIn(f"ruff format --check {paths}", workflow)

    def test_pull_request_checks_are_not_restricted_to_a_base_branch(self) -> None:
        # Issue #53: a `pull_request` `branches:` filter matches the base
        # branch, so restricting it to `main` silently gave stacked PRs no
        # check run at all. Re-adding such a filter would remove CI from
        # exactly the PRs that most need it, without failing anything
        # visibly, so both workflows are pinned here.
        for workflow_name in (".github/workflows/lint.yml", ".github/workflows/test.yml"):
            with self.subTest(workflow=workflow_name):
                workflow = (ROOT / workflow_name).read_text(encoding="utf-8")
                trigger_start = workflow.index("  pull_request:")
                trigger_block = workflow[trigger_start : workflow.index("jobs:", trigger_start)]
                self.assertNotIn(
                    "branches:",
                    trigger_block,
                    f"{workflow_name} restricts pull_request to a base branch again",
                )

    def test_readme_uses_one_nf_test_invocation_form(self) -> None:
        # Issue #53: the README mixed `./nf-test` and `nf-test`, which is
        # only reproducible if the reader happens to install it the same
        # way. CI puts nf-test on PATH, so the PATH form is the contract.
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertNotIn("./nf-test", readme)
        self.assertIn("nf-test test", readme)

    def test_local_secrets_and_caches_are_ignored_but_examples_are_not(self) -> None:
        ignored = (
            ".env",
            ".env.local",
            ".venv/bin/python",
            "venv/bin/python",
            ".ruff_cache/state.json",
            ".mypy_cache/state.json",
            ".pytest_cache/state.json",
            "bin/__pycache__/tool.cpython-312.pyc",
        )
        for path in ignored:
            with self.subTest(path=path):
                result = subprocess.run(
                    ["git", "check-ignore", "--no-index", "--quiet", path],
                    cwd=ROOT,
                    check=False,
                )
                self.assertEqual(0, result.returncode)

        retained = (
            ".env.example",
            ".env.test.example",
            "tests/data/reads/sample_R1.fastq.gz",
            "tests/bin/fixtures/sample.vcf.gz",
        )
        for path in retained:
            with self.subTest(path=path):
                result = subprocess.run(
                    ["git", "check-ignore", "--no-index", "--quiet", path],
                    cwd=ROOT,
                    check=False,
                )
                self.assertEqual(1, result.returncode)


if __name__ == "__main__":
    unittest.main()
