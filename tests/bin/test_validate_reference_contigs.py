"""Unit tests for bin/validate_reference_contigs.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "validate_reference_contigs.py"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


validate_module = _load_module("validate_reference_contigs", SCRIPT_PATH)


def _write_fai(path: Path, rows: list[tuple[str, int]]) -> None:
    lines = [f"{name}\t{length}\t0\t70\t71" for name, length in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_dict(path: Path, rows: list[tuple[str, int]]) -> None:
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    lines += [f"@SQ\tSN:{name}\tLN:{length}\tM5:deadbeef\tUR:file:./ref.fa" for name, length in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


_STANDARD_CONTIGS = [("Chr01", 1000), ("Chr02", 2000), ("Chr03", 3000)]


class ParseFaiTests(unittest.TestCase):
    def test_parses_name_and_length_in_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            _write_fai(fai_path, _STANDARD_CONTIGS)
            records = validate_module.parse_fai(fai_path)
        self.assertEqual(
            [(record.name, record.length) for record in records],
            _STANDARD_CONTIGS,
        )

    def test_rejects_non_integer_length(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            fai_path.write_text("Chr01\tnot-a-number\t0\t70\t71\n", encoding="utf-8")
            with self.assertRaises(validate_module.MalformedFaiError):
                validate_module.parse_fai(fai_path)

    def test_rejects_too_few_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            fai_path.write_text("Chr01\n", encoding="utf-8")
            with self.assertRaises(validate_module.MalformedFaiError):
                validate_module.parse_fai(fai_path)

    def test_rejects_duplicate_contig_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            _write_fai(fai_path, [("Chr01", 1000), ("Chr01", 2000)])
            with self.assertRaises(validate_module.DuplicateContigError):
                validate_module.parse_fai(fai_path)

    def test_skips_blank_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            fai_path.write_text("Chr01\t1000\t0\t70\t71\n\n\nChr02\t2000\t0\t70\t71\n", encoding="utf-8")
            records = validate_module.parse_fai(fai_path)
        self.assertEqual(len(records), 2)


class ParseDictTests(unittest.TestCase):
    def test_parses_sn_and_ln_in_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            _write_dict(dict_path, _STANDARD_CONTIGS)
            records = validate_module.parse_dict(dict_path)
        self.assertEqual(
            [(record.name, record.length) for record in records],
            _STANDARD_CONTIGS,
        )

    def test_ignores_non_sq_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            dict_path.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:Chr01\tLN:1000\n@CO\tsome comment\n",
                encoding="utf-8",
            )
            records = validate_module.parse_dict(dict_path)
        self.assertEqual(records, [validate_module.ContigRecord(name="Chr01", length=1000)])

    def test_rejects_missing_sn(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            dict_path.write_text("@SQ\tLN:1000\n", encoding="utf-8")
            with self.assertRaises(validate_module.MalformedDictError):
                validate_module.parse_dict(dict_path)

    def test_rejects_missing_ln(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            dict_path.write_text("@SQ\tSN:Chr01\n", encoding="utf-8")
            with self.assertRaises(validate_module.MalformedDictError):
                validate_module.parse_dict(dict_path)

    def test_rejects_non_integer_ln(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            dict_path.write_text("@SQ\tSN:Chr01\tLN:not-a-number\n", encoding="utf-8")
            with self.assertRaises(validate_module.MalformedDictError):
                validate_module.parse_dict(dict_path)

    def test_rejects_duplicate_contig_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dict_path = Path(tmp) / "ref.dict"
            _write_dict(dict_path, [("Chr01", 1000), ("Chr01", 2000)])
            with self.assertRaises(validate_module.DuplicateContigError):
                validate_module.parse_dict(dict_path)


class FindFirstMismatchTests(unittest.TestCase):
    def test_identical_lists_have_no_mismatch(self) -> None:
        fai = [validate_module.ContigRecord(name=n, length=l) for n, l in _STANDARD_CONTIGS]
        dict_records = [validate_module.ContigRecord(name=n, length=l) for n, l in _STANDARD_CONTIGS]
        self.assertIsNone(validate_module.find_first_mismatch(fai, dict_records))

    def test_reports_count_mismatch(self) -> None:
        fai = [validate_module.ContigRecord(name="Chr01", length=1000)]
        dict_records = [
            validate_module.ContigRecord(name="Chr01", length=1000),
            validate_module.ContigRecord(name="Chr02", length=2000),
        ]
        message = validate_module.find_first_mismatch(fai, dict_records)
        assert message is not None
        self.assertIn("contig count differs", message)
        self.assertIn("1 contig", message)
        self.assertIn("2 contig", message)

    def test_reports_first_mismatch_index_on_name_difference(self) -> None:
        fai = [
            validate_module.ContigRecord(name="Chr01", length=1000),
            validate_module.ContigRecord(name="Chr02", length=2000),
        ]
        dict_records = [
            validate_module.ContigRecord(name="Chr01", length=1000),
            validate_module.ContigRecord(name="Chr03", length=2000),
        ]
        message = validate_module.find_first_mismatch(fai, dict_records)
        assert message is not None
        self.assertIn("index 1", message)
        self.assertIn("Chr02", message)
        self.assertIn("Chr03", message)

    def test_reports_length_mismatch_even_when_names_match(self) -> None:
        fai = [validate_module.ContigRecord(name="Chr01", length=1000)]
        dict_records = [validate_module.ContigRecord(name="Chr01", length=999)]
        message = validate_module.find_first_mismatch(fai, dict_records)
        assert message is not None
        self.assertIn("index 0", message)
        self.assertIn("length=1000", message)
        self.assertIn("length=999", message)

    def test_same_contig_set_but_different_order_is_a_mismatch(self) -> None:
        # This is the critical case a naive set-equality comparison would miss.
        fai = [validate_module.ContigRecord(name=n, length=l) for n, l in _STANDARD_CONTIGS]
        reordered = [_STANDARD_CONTIGS[1], _STANDARD_CONTIGS[0], _STANDARD_CONTIGS[2]]
        dict_records = [validate_module.ContigRecord(name=n, length=l) for n, l in reordered]
        message = validate_module.find_first_mismatch(fai, dict_records)
        self.assertIsNotNone(message)


class MainCliTests(unittest.TestCase):
    def _run(self, fai_rows: list[tuple[str, int]], dict_rows: list[tuple[str, int]]) -> tuple[int, str]:
        with tempfile.TemporaryDirectory() as tmp:
            fai_path = Path(tmp) / "ref.fa.fai"
            dict_path = Path(tmp) / "ref.dict"
            _write_fai(fai_path, fai_rows)
            _write_dict(dict_path, dict_rows)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = validate_module.main(
                    ["--fai", str(fai_path), "--dict", str(dict_path)]
                )
            return exit_code, stderr.getvalue()

    def test_matching_generated_style_pair_succeeds(self) -> None:
        exit_code, stderr = self._run(_STANDARD_CONTIGS, _STANDARD_CONTIGS)
        self.assertEqual(exit_code, 0)
        self.assertEqual(stderr, "")

    def test_wrong_order_fails_with_actionable_message(self) -> None:
        reordered = [_STANDARD_CONTIGS[2], _STANDARD_CONTIGS[1], _STANDARD_CONTIGS[0]]
        exit_code, stderr = self._run(_STANDARD_CONTIGS, reordered)
        self.assertEqual(exit_code, 1)
        self.assertIn("inconsistent", stderr)
        self.assertIn("first mismatch at index", stderr)

    def test_length_mismatch_fails(self) -> None:
        wrong_length = [("Chr01", 1000), ("Chr02", 9999), ("Chr03", 3000)]
        exit_code, stderr = self._run(_STANDARD_CONTIGS, wrong_length)
        self.assertEqual(exit_code, 1)
        self.assertIn("inconsistent", stderr)

    def test_missing_contig_fails(self) -> None:
        missing_one = _STANDARD_CONTIGS[:2]
        exit_code, stderr = self._run(_STANDARD_CONTIGS, missing_one)
        self.assertEqual(exit_code, 1)
        self.assertIn("contig count differs", stderr)

    def test_extra_contig_fails(self) -> None:
        extra_one = [*_STANDARD_CONTIGS, ("Chr04", 4000)]
        exit_code, stderr = self._run(_STANDARD_CONTIGS, extra_one)
        self.assertEqual(exit_code, 1)
        self.assertIn("contig count differs", stderr)

    def test_empty_fai_fails(self) -> None:
        exit_code, stderr = self._run([], _STANDARD_CONTIGS)
        self.assertEqual(exit_code, 1)
        self.assertIn("contains no contigs", stderr)

    def test_empty_dict_fails(self) -> None:
        exit_code, stderr = self._run(_STANDARD_CONTIGS, [])
        self.assertEqual(exit_code, 1)
        self.assertIn("contains no @SQ records", stderr)


if __name__ == "__main__":
    unittest.main()
