"""Unit tests for bin/manifest_utils.py, the shared manifest mechanics.

Run with: python3 -m unittest discover -s tests/bin -v

Issue #42 extracted seven functions that `bin/build_run_manifest.py`
and `bin/build_gs_panel_manifest.py` had each defined identically, plus
the GS manifest's container-identity leakage guard, into one module.
Two things have to stay true after that move, and this file pins both:

  * both manifests must go on producing byte-identical documents --
    hence the golden `manifest_hash` regressions below, computed by
    running the *pre-refactor* implementations (main@402f437) against
    fixed inputs with a fixed run_id/timestamp;
  * both scripts must be importing the same objects, not their own
    surviving copies -- hence the identity assertions.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
import types
import unittest
from datetime import UTC, datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# See tests/bin/test_build_run_manifest.py for why bin/ goes on sys.path
# here but needs no PYTHONPATH in production.
sys.path.insert(0, str(REPO_ROOT / "bin"))

import manifest_utils  # noqa: E402


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run_manifest_module = _load_module(
    "build_run_manifest_for_utils_tests", REPO_ROOT / "bin" / "build_run_manifest.py"
)
gs_manifest_module = _load_module(
    "build_gs_panel_manifest_for_utils_tests",
    REPO_ROOT / "bin" / "build_gs_panel_manifest.py",
)

SHARED_FUNCTION_NAMES = (
    "new_run_id",
    "utc_now_iso",
    "sha256_file",
    "checksum_files",
    "_json_default",
    "canonical_json_hash",
    "write_json_atomic",
)

# Fixed inputs + fixed run_id/generated_at, so the resulting document --
# and therefore its hash -- is fully determined by the manifest-building
# code alone. The expected hashes were produced by running the
# pre-refactor scripts (main@402f437, before manifest_utils.py existed)
# against exactly these arguments.
FIXED_RUN_ID = "20260814T000000Z-deadbeef"
FIXED_GENERATED_AT = "2026-08-14T00:00:00Z"

PRE_REFACTOR_RUN_MANIFEST_HASH = (
    "sha256:fa9991dd8ed5878765f5b1db63f561698bd0de57bc5be8e5440757868c1fb765"
)
PRE_REFACTOR_GS_MANIFEST_HASH = (
    "sha256:81bb294dd4d19a93e090ebac9e7eaaaed867ce45f4fda20b8698149328f54545"
)


def _build_run_manifest_fixture() -> dict[str, object]:
    return run_manifest_module.build_manifest(
        cohort_id="cohort",
        pipeline_version="0.2.0",
        git_commit="",
        nextflow_version="26.04.6",
        containers={
            "bwa_mem2": "bwa:1",
            "samtools": "sam:1",
            "gatk": "gatk:1",
            "python": "py:1",
        },
        reference={
            "reference_id": "ref",
            "fasta": {"filename": "r.fa", "checksum": "sha256:aa"},
        },
        parameters={"sample_ploidy": 2, "enable_gs_panel": True},
        samples=[{"sample_id": "sample_a", "read_group_id": "rg1"}],
        cohort_accounting={"number_of_samples": "2"},
        variant_type_accounting={"raw_all_records": "2"},
        gs_panel=None,
        checksums={"a.vcf.gz": "sha256:bb"},
        run_id=FIXED_RUN_ID,
        generated_at=FIXED_GENERATED_AT,
    )


def _build_gs_manifest_fixture() -> dict[str, object]:
    return gs_manifest_module.build_manifest(
        cohort_id="cohort",
        pipeline_version="0.2.0",
        git_commit="",
        containers={"gs_normalize_variants": "bcftools:1"},
        sample_ploidy=2,
        snp_filter_params={"snp_filter_qd_min": 2.0},
        panel_status="empty",
        checksums={"m.tsv.gz": "sha256:cc"},
        run_id=FIXED_RUN_ID,
        generated_at=FIXED_GENERATED_AT,
    )


class SharedImplementationTests(unittest.TestCase):
    """Both scripts must use the shared functions, not private copies."""

    def test_run_manifest_uses_the_shared_functions(self) -> None:
        for name in SHARED_FUNCTION_NAMES:
            with self.subTest(function=name):
                self.assertIs(
                    getattr(run_manifest_module, name), getattr(manifest_utils, name)
                )

    def test_gs_manifest_uses_the_shared_functions(self) -> None:
        for name in SHARED_FUNCTION_NAMES:
            with self.subTest(function=name):
                self.assertIs(
                    getattr(gs_manifest_module, name), getattr(manifest_utils, name)
                )

    def test_gs_manifest_uses_the_shared_container_identity_guard(self) -> None:
        self.assertIs(
            gs_manifest_module.validate_container_identity,
            manifest_utils.validate_container_identity,
        )


class RefactorRegressionTests(unittest.TestCase):
    """The extracted mechanics must not have changed either document."""

    def test_run_manifest_payload_is_unchanged(self) -> None:
        manifest = _build_run_manifest_fixture()

        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["run_id"], FIXED_RUN_ID)
        self.assertEqual(manifest["generated_at"], FIXED_GENERATED_AT)
        self.assertIsNone(manifest["git_commit"])
        self.assertEqual(
            manifest["containers"],
            {"bwa_mem2": "bwa:1", "samtools": "sam:1", "gatk": "gatk:1", "python": "py:1"},
        )

    def test_run_manifest_hash_matches_pre_refactor_implementation(self) -> None:
        self.assertEqual(
            _build_run_manifest_fixture()["manifest_hash"],
            PRE_REFACTOR_RUN_MANIFEST_HASH,
        )

    def test_gs_manifest_hash_matches_pre_refactor_implementation(self) -> None:
        self.assertEqual(
            _build_gs_manifest_fixture()["manifest_hash"],
            PRE_REFACTOR_GS_MANIFEST_HASH,
        )

    def test_manifest_hash_is_over_the_document_without_its_own_hash(self) -> None:
        # The self-referential hash both manifests carry: recomputing it
        # over the document minus the hash field must reproduce it.
        manifest = _build_run_manifest_fixture()
        without_hash = {
            key: value for key, value in manifest.items() if key != "manifest_hash"
        }
        self.assertEqual(
            manifest_utils.canonical_json_hash(without_hash), manifest["manifest_hash"]
        )


class ChecksumFilesTests(unittest.TestCase):
    def test_checksums_are_keyed_by_filename_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            nested = Path(tmp) / "deep" / "nested"
            nested.mkdir(parents=True)
            target = nested / "artifact.tsv"
            target.write_text("payload\n", encoding="utf-8")

            checksums = manifest_utils.checksum_files([target])

            self.assertEqual(list(checksums), ["artifact.tsv"])
            self.assertNotIn(tmp, json.dumps(checksums))
            self.assertTrue(checksums["artifact.tsv"].startswith("sha256:"))

    def test_duplicate_filename_from_different_directories_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "a"
            second = Path(tmp) / "b"
            first.mkdir()
            second.mkdir()
            (first / "same.tsv").write_text("one\n", encoding="utf-8")
            (second / "same.tsv").write_text("two\n", encoding="utf-8")

            with self.assertRaises(ValueError) as raised:
                manifest_utils.checksum_files(
                    [first / "same.tsv", second / "same.tsv"]
                )

            self.assertIn("duplicate checksum file name", str(raised.exception))


class WriteJsonAtomicTests(unittest.TestCase):
    """The atomic-write contract: no partial file, no temp file, no clobber."""

    def test_writes_readable_json_with_trailing_newline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "manifest.json"
            manifest_utils.write_json_atomic(output, {"b": 2, "a": 1})

            text = output.read_text(encoding="utf-8")
            self.assertTrue(text.endswith("\n"))
            self.assertEqual(json.loads(text), {"a": 1, "b": 2})

    def test_leaves_no_temp_file_on_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "manifest.json"
            manifest_utils.write_json_atomic(output, {"a": 1})

            self.assertEqual([p.name for p in Path(tmp).iterdir()], ["manifest.json"])

    def test_unserializable_payload_writes_nothing_at_all(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "manifest.json"

            with self.assertRaises(TypeError):
                manifest_utils.write_json_atomic(output, {"bad": object()})

            self.assertFalse(output.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_failed_write_leaves_no_temp_file_and_preserves_existing_final(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "manifest.json"
            existing = {"schema_version": 1, "run_id": "previous"}
            manifest_utils.write_json_atomic(output, existing)

            original_replace = os.replace

            def failing_replace(source: object, destination: object) -> None:
                raise OSError("simulated rename failure")

            manifest_utils.os.replace = failing_replace
            try:
                with self.assertRaises(OSError):
                    manifest_utils.write_json_atomic(output, {"run_id": "replacement"})
            finally:
                manifest_utils.os.replace = original_replace

            # The previous manifest is intact...
            self.assertEqual(json.loads(output.read_text(encoding="utf-8")), existing)
            # ...and no partial temp file was left beside it.
            self.assertEqual(
                sorted(p.name for p in Path(tmp).iterdir()), ["manifest.json"]
            )


class TimeHelperTests(unittest.TestCase):
    def test_run_id_is_sortable_and_uses_the_injected_moment(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        self.assertEqual(
            manifest_utils.new_run_id(moment, "abcd1234"), "20260814T123456Z-abcd1234"
        )

    def test_utc_now_iso_uses_the_injected_moment(self) -> None:
        moment = datetime(2026, 8, 14, 12, 34, 56, tzinfo=UTC)
        self.assertEqual(manifest_utils.utc_now_iso(moment), "2026-08-14T12:34:56Z")


if __name__ == "__main__":
    unittest.main()
