"""Unit tests for bin/build_gs_panel_manifest.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import re
import subprocess
import sys
import tempfile
import types
import unittest
from datetime import UTC, datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "build_gs_panel_manifest.py"

RUN_ID_PATTERN = re.compile(r"^\d{8}T\d{6}Z-[0-9a-f]{8}$")

PLOIDY_CLI_ARGS = ["--sample-ploidy", "2"]

SNP_FILTER_CLI_ARGS = [
    "--snp-filter-qd-min", "2.0",
    "--snp-filter-qual-min", "30.0",
    "--snp-filter-sor-max", "3.0",
    "--snp-filter-fs-max", "60.0",
    "--snp-filter-mq-min", "40.0",
    "--snp-filter-mq-rank-sum-min", "-12.5",
    "--snp-filter-read-pos-rank-sum-min", "-8.0",
]

# Issue #52: one --container-<process-name> flag per GS-lineage process
# (schema v2), replacing the old --bcftools-container/--gatk-container/
# --python-container triple. Deliberately not all sharing one image
# string, so a test bug that accidentally reads the wrong process's value
# would show up as a mismatched assertion rather than passing by accident.
# `build_gs_panel_manifest` is the manifest-writing process itself, which
# records its own container too (BUILD_GS_PANEL_MANIFEST passes its own
# task.container from inside its own script).
CONTAINER_CLI_ARGS = [
    "--container-gs-normalize-variants", "bcftools:1.24-a",
    "--container-classify-normalized-variants", "python:3.12-b",
    "--container-gs-index-classified-variants", "bcftools:1.24-c",
    "--container-gatk-variantfiltration-gs", "gatk:4.6.2.0-d",
    "--container-gatk-selectpassvariants-gs", "gatk:4.6.2.0-e",
    "--container-build-gs-panel", "python:3.12-f",
    "--container-reconcile-gs-panel-accounting", "python:3.12-g",
    "--container-build-gs-panel-manifest", "python:3.12-h",
]

CONTAINERS_KWARG = {
    "gs_normalize_variants": "bcftools:1.24-a",
    "classify_normalized_variants": "python:3.12-b",
    "gs_index_classified_variants": "bcftools:1.24-c",
    "gatk_variantfiltration_gs": "gatk:4.6.2.0-d",
    "gatk_selectpassvariants_gs": "gatk:4.6.2.0-e",
    "build_gs_panel": "python:3.12-f",
    "reconcile_gs_panel_accounting": "python:3.12-g",
    "build_gs_panel_manifest": "python:3.12-h",
}


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


manifest_module = _load_module("build_gs_panel_manifest", SCRIPT_PATH)


def _write_record_accounting(path: Path, panel_status: str) -> None:
    rows = [
        "cohort_id\tmetric\tvalue",
        f"cohort\tpanel_status\t{panel_status}",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


class NewRunIdTests(unittest.TestCase):
    """Tests for new_run_id: format and injectability."""

    def test_format_with_injected_time_and_suffix(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        run_id = manifest_module.new_run_id(now=moment, suffix="a1b2c3d4")
        self.assertEqual(run_id, "20260814T123456Z-a1b2c3d4")

    def test_real_call_matches_expected_pattern(self) -> None:
        run_id = manifest_module.new_run_id()
        self.assertRegex(run_id, RUN_ID_PATTERN)


class UtcNowIsoTests(unittest.TestCase):
    """Tests for utc_now_iso: ISO-8601 formatting with injected time."""

    def test_format_with_injected_time(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        self.assertEqual(manifest_module.utc_now_iso(now=moment), "2026-08-14T12:34:56Z")


class ChecksumFilesTests(unittest.TestCase):
    """Tests for checksum_files: filename-only keys and duplicate rejection."""

    def test_checksums_keyed_by_filename_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matrix.tsv.gz"
            path.write_bytes(b"hello")

            checksums = manifest_module.checksum_files([path])

            self.assertIn("matrix.tsv.gz", checksums)
            self.assertTrue(checksums["matrix.tsv.gz"].startswith("sha256:"))
            self.assertNotIn(str(path), checksums["matrix.tsv.gz"])

    def test_duplicate_filename_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_a, tempfile.TemporaryDirectory() as tmp_b:
            path_a = Path(tmp_a) / "same_name.tsv"
            path_b = Path(tmp_b) / "same_name.tsv"
            path_a.write_bytes(b"a")
            path_b.write_bytes(b"b")

            with self.assertRaises(ValueError):
                manifest_module.checksum_files([path_a, path_b])


class ReadPanelStatusTests(unittest.TestCase):
    """Tests for read_panel_status."""

    def test_reads_panel_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            _write_record_accounting(path, "populated")

            self.assertEqual(manifest_module.read_panel_status(path), "populated")

    def test_missing_metric_raises_with_cause(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accounting.tsv"
            path.write_text("cohort_id\tmetric\tvalue\n", encoding="utf-8")

            with self.assertRaises(manifest_module.MalformedAccountingError) as raised:
                manifest_module.read_panel_status(path)
            self.assertIn("panel_status", str(raised.exception))


class BuildManifestTests(unittest.TestCase):
    """Tests for build_manifest: schema, git_commit null handling, and hash."""

    def _build(self, **overrides: object) -> dict[str, object]:
        defaults = dict(
            cohort_id="cohort",
            pipeline_version="0.2.0-dev",
            git_commit="",
            containers=dict(CONTAINERS_KWARG),
            sample_ploidy=2,
            snp_filter_params={"snp_filter_qd_min": 2.0},
            panel_status="populated",
            checksums={},
            run_id="20260814T000000Z-deadbeef",
            generated_at="2026-08-14T00:00:00Z",
        )
        defaults.update(overrides)
        return manifest_module.build_manifest(**defaults)

    def test_schema_version_is_two(self) -> None:
        self.assertEqual(self._build()["schema_version"], 2)

    def test_recorded_processes_are_the_whole_gs_lineage_including_this_one(
        self,
    ) -> None:
        # Issue #52 review (P2): the manifest's `containers` contract is
        # "every process the GS lineage runs", which includes the process
        # that writes the manifest. A regression that quietly dropped the
        # self-entry would otherwise still satisfy every other test here.
        self.assertEqual(
            set(manifest_module.CONTAINER_PROCESS_NAMES),
            set(CONTAINERS_KWARG),
        )
        self.assertIn(
            "build_gs_panel_manifest", manifest_module.CONTAINER_PROCESS_NAMES
        )
        self.assertEqual(len(manifest_module.CONTAINER_PROCESS_NAMES), 8)
        self.assertIn("build_gs_panel_manifest", self._build()["containers"])

    def test_containers_are_recorded_per_process_not_merged_by_tool(self) -> None:
        manifest = self._build()
        self.assertEqual(manifest["containers"], CONTAINERS_KWARG)
        # gs_normalize_variants and gs_index_classified_variants both use
        # bcftools by default, but a divergent override (as configured
        # here) must be visible as two distinct entries, never merged.
        self.assertNotEqual(
            manifest["containers"]["gs_normalize_variants"],
            manifest["containers"]["gs_index_classified_variants"],
        )

    def test_empty_git_commit_becomes_null(self) -> None:
        self.assertIsNone(self._build(git_commit="")["git_commit"])

    def test_nonempty_git_commit_is_preserved(self) -> None:
        self.assertEqual(self._build(git_commit="abc123")["git_commit"], "abc123")

    def test_manifest_hash_is_deterministic_for_identical_content(self) -> None:
        first = self._build()
        second = self._build()
        self.assertEqual(first["manifest_hash"], second["manifest_hash"])

    def test_manifest_hash_changes_when_content_changes(self) -> None:
        first = self._build(panel_status="populated")
        second = self._build(panel_status="empty")
        self.assertNotEqual(first["manifest_hash"], second["manifest_hash"])

    def test_manifest_hash_is_not_included_in_its_own_hash_input(self) -> None:
        manifest = self._build()
        recomputed = manifest_module.canonical_json_hash(
            {key: value for key, value in manifest.items() if key != "manifest_hash"}
        )
        self.assertEqual(manifest["manifest_hash"], recomputed)

    def test_sample_ploidy_is_recorded_in_parameters(self) -> None:
        manifest = self._build(sample_ploidy=2)
        self.assertEqual(manifest["parameters"]["sample_ploidy"], 2)

    def test_genotype_encoding_schema_is_recorded_verbatim(self) -> None:
        manifest = self._build()
        self.assertEqual(manifest["genotype_encoding"], manifest_module.GENOTYPE_ENCODING_SCHEMA)
        self.assertEqual(manifest["genotype_encoding"]["missing_token"], "nan")
        self.assertEqual(manifest["genotype_encoding"]["ploidy"], "diploid_only")


class WriteJsonAtomicTests(unittest.TestCase):
    """Tests for write_json_atomic: no temp file left behind, valid JSON written."""

    def test_writes_valid_json_and_leaves_no_temp_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "manifest.json"
            manifest_module.write_json_atomic(output_path, {"a": 1})

            self.assertEqual(json.loads(output_path.read_text(encoding="utf-8")), {"a": 1})
            leftover_temp_files = list(Path(tmp).glob(".*.tmp"))
            self.assertEqual(leftover_temp_files, [])


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_valid_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "populated")
            matrix_path = tmp_path / "matrix.tsv.gz"
            matrix_path.write_bytes(b"fake matrix content")
            output_path = tmp_path / "manifest.json"

            exit_code = manifest_module.main(
                [
                    "--cohort-id", "cohort",
                    "--pipeline-version", "0.2.0-dev",
                    "--git-commit", "",
                    *CONTAINER_CLI_ARGS,
                    *PLOIDY_CLI_ARGS,
                    *SNP_FILTER_CLI_ARGS,
                    "--record-accounting", str(accounting_path),
                    "--checksum-file", str(matrix_path),
                    "--output", str(output_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            manifest = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema_version"], 2)
            self.assertIsNone(manifest["git_commit"])
            self.assertEqual(manifest["panel_status"], "populated")
            self.assertIn("matrix.tsv.gz", manifest["checksums"])
            self.assertRegex(manifest["run_id"], RUN_ID_PATTERN)
            self.assertEqual(manifest["parameters"]["snp_filter_qd_min"], 2.0)
            self.assertEqual(manifest["parameters"]["sample_ploidy"], 2)
            self.assertEqual(manifest["genotype_encoding"], manifest_module.GENOTYPE_ENCODING_SCHEMA)
            self.assertEqual(manifest["containers"], CONTAINERS_KWARG)

    def test_main_fails_clearly_for_a_missing_accounting_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_accounting = tmp_path / "does-not-exist.tsv"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(
                    [
                        "--cohort-id", "cohort",
                        "--pipeline-version", "0.2.0-dev",
                        *CONTAINER_CLI_ARGS,
                        *PLOIDY_CLI_ARGS,
                        *SNP_FILTER_CLI_ARGS,
                        "--record-accounting", str(missing_accounting),
                        "--output", str(tmp_path / "manifest.json"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_accounting), stderr.getvalue())

    def _run_main_with_ploidy(self, ploidy: str) -> tuple[int, Path, str]:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "populated")
            output_path = tmp_path / "manifest.json"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(
                    [
                        "--cohort-id", "cohort",
                        "--pipeline-version", "0.2.0-dev",
                        *CONTAINER_CLI_ARGS,
                        "--sample-ploidy", ploidy,
                        *SNP_FILTER_CLI_ARGS,
                        "--record-accounting", str(accounting_path),
                        "--output", str(output_path),
                    ]
                )

            return exit_code, output_path, stderr.getvalue()

    def test_main_fails_fast_for_haploid_ploidy(self) -> None:
        exit_code, output_path, stderr = self._run_main_with_ploidy("1")

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("diploid", stderr)

    def test_main_fails_fast_for_polyploid_ploidy(self) -> None:
        exit_code, output_path, stderr = self._run_main_with_ploidy("3")

        self.assertEqual(exit_code, 1)
        self.assertFalse(output_path.exists())
        self.assertIn("diploid", stderr)

    def test_ploidy_error_names_the_current_schema_version(self) -> None:
        # Issue #52 review (P3): this message hardcoded "(v1)" while
        # SCHEMA_VERSION had already moved to 2. Asserting against the
        # constant (rather than a literal "v2") keeps the message and the
        # schema version from ever being two independently edited facts.
        _exit_code, _output_path, stderr = self._run_main_with_ploidy("3")

        self.assertIn(f"v{manifest_module.SCHEMA_VERSION}", stderr)
        self.assertNotIn("(v1)", stderr)

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "empty")
            output_path = tmp_path / "manifest.json"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--cohort-id", "cohort",
                    "--pipeline-version", "0.2.0-dev",
                    *CONTAINER_CLI_ARGS,
                    *PLOIDY_CLI_ARGS,
                    *SNP_FILTER_CLI_ARGS,
                    "--record-accounting", str(accounting_path),
                    "--output", str(output_path),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["panel_status"], "empty")

    def test_main_fails_fast_when_a_container_identity_is_a_host_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "empty")
            output_path = tmp_path / "manifest.json"
            stderr = io.StringIO()

            args = [
                "--cohort-id", "cohort",
                "--pipeline-version", "0.2.0-dev",
                *CONTAINER_CLI_ARGS,
                *PLOIDY_CLI_ARGS,
                *SNP_FILTER_CLI_ARGS,
                "--record-accounting", str(accounting_path),
                "--output", str(output_path),
            ]
            # Overwrite one container identity with a Singularity-style
            # absolute host path -- exactly the kind of value that must
            # never reach a manifest published as a shareable artifact.
            index = args.index("--container-build-gs-panel") + 1
            args[index] = "/opt/images/python.sif"

            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(args)

            self.assertEqual(exit_code, 1)
            self.assertFalse(output_path.exists())
            self.assertIn("build_gs_panel", stderr.getvalue())
            self.assertIn("host filesystem path", stderr.getvalue())

    def test_main_fails_fast_when_a_container_identity_is_a_file_uri(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "empty")
            output_path = tmp_path / "manifest.json"
            stderr = io.StringIO()

            args = [
                "--cohort-id", "cohort",
                "--pipeline-version", "0.2.0-dev",
                *CONTAINER_CLI_ARGS,
                *PLOIDY_CLI_ARGS,
                *SNP_FILTER_CLI_ARGS,
                "--record-accounting", str(accounting_path),
                "--output", str(output_path),
            ]
            # Issue #52 review (P1): a Singularity/Apptainer image path
            # wrapped in a file:// URI discloses host filesystem layout
            # exactly as fully as the bare path above, but starts with
            # none of `/`, `./`, `../`, `~` -- so the prefix check alone
            # let it straight through into a shareable manifest.
            index = args.index("--container-gatk-variantfiltration-gs") + 1
            args[index] = "file:///opt/images/gatk.sif"

            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(args)

            self.assertEqual(exit_code, 1)
            self.assertFalse(output_path.exists())
            self.assertIn("gatk_variantfiltration_gs", stderr.getvalue())
            self.assertIn("file://", stderr.getvalue())
            self.assertIn("host filesystem path", stderr.getvalue())

    def test_main_fails_fast_when_a_container_identity_embeds_credentials(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            accounting_path = tmp_path / "accounting.tsv"
            _write_record_accounting(accounting_path, "empty")
            output_path = tmp_path / "manifest.json"
            stderr = io.StringIO()

            args = [
                "--cohort-id", "cohort",
                "--pipeline-version", "0.2.0-dev",
                *CONTAINER_CLI_ARGS,
                *PLOIDY_CLI_ARGS,
                *SNP_FILTER_CLI_ARGS,
                "--record-accounting", str(accounting_path),
                "--output", str(output_path),
            ]
            index = args.index("--container-reconcile-gs-panel-accounting") + 1
            args[index] = "https://user:hunter2@example.com/python:3.12"

            with contextlib.redirect_stderr(stderr):
                exit_code = manifest_module.main(args)

            self.assertEqual(exit_code, 1)
            self.assertFalse(output_path.exists())
            self.assertIn("credentials", stderr.getvalue())


class ValidateContainerIdentityTests(unittest.TestCase):
    """Tests for validate_container_identity: the manifest's own leakage guard."""

    def test_accepts_a_plain_digest_pinned_image_reference(self) -> None:
        value = "quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:" + "a" * 64
        self.assertEqual(
            manifest_module.validate_container_identity("gs_normalize_variants", value),
            value,
        )

    def test_rejects_empty_value(self) -> None:
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity("build_gs_panel", "")

    def test_rejects_whitespace_only_value(self) -> None:
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity("build_gs_panel", "   ")

    def test_rejects_absolute_host_path(self) -> None:
        with self.assertRaises(ValueError) as raised:
            manifest_module.validate_container_identity(
                "gatk_variantfiltration_gs", "/data/images/gatk.sif"
            )
        self.assertIn("host filesystem path", str(raised.exception))

    def test_rejects_relative_host_path(self) -> None:
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity(
                "gatk_variantfiltration_gs", "./gatk.sif"
            )

    def test_rejects_home_relative_host_path(self) -> None:
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity(
                "gatk_variantfiltration_gs", "~/images/gatk.sif"
            )

    def test_rejects_file_uri_host_path(self) -> None:
        with self.assertRaises(ValueError) as raised:
            manifest_module.validate_container_identity(
                "gatk_variantfiltration_gs", "file:///opt/images/gatk.sif"
            )
        message = str(raised.exception)
        self.assertIn("file://", message)
        self.assertIn("host filesystem path", message)
        self.assertIn("/opt/images/gatk.sif", message)

    def test_rejects_file_uri_home_directory_path(self) -> None:
        with self.assertRaises(ValueError) as raised:
            manifest_module.validate_container_identity(
                "classify_normalized_variants",
                "file:///home/user/containers/python.sif",
            )
        self.assertIn("host filesystem path", str(raised.exception))

    def test_rejects_file_uri_regardless_of_scheme_case(self) -> None:
        # A `file://` URI is scheme-case-insensitive per RFC 3986, so a
        # `FILE://`-cased value leaks exactly as much host layout.
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity(
                "build_gs_panel_manifest", "FILE:///opt/images/python.sif"
            )

    def test_rejects_embedded_registry_credentials(self) -> None:
        with self.assertRaises(ValueError) as raised:
            manifest_module.validate_container_identity(
                "classify_normalized_variants",
                "https://svc-account:s3cr3t@registry.internal/python:3.12",
            )
        self.assertIn("credentials", str(raised.exception))

    def test_rejects_value_containing_whitespace(self) -> None:
        with self.assertRaises(ValueError):
            manifest_module.validate_container_identity(
                "build_gs_panel", "python:3.12 (local build)"
            )


if __name__ == "__main__":
    unittest.main()
