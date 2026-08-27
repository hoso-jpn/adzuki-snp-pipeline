"""Unit and structural regression tests for the locus-stream classifier."""
from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import stat
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "classify_normalized_variants.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


classifier = _load_module("classify_normalized_variants", SCRIPT_PATH)
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FILTER=<ID=Old,Description=\"stale\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b\n"
)


def _row(chrom: str, pos: int, ref: str = "A", alt: str = "G",
         filter_: str = "Old") -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filter_}\tDP=9\tGT\t0/1\t1/1\n"


class Harness:
    def run_text(self, text: str):
        self.directory = tempfile.TemporaryDirectory()
        root = Path(self.directory.name)
        source = root / "input.vcf.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(text)
        paths = root / "output.vcf", root / "accounting.tsv", root / "summary.txt"
        code = classifier.main([
            "--normalized-vcf", str(source), "--cohort-id", "cohort",
            "--output", str(paths[0]), "--accounting-output", str(paths[1]),
            "--summary-output", str(paths[2]),
        ])
        return code, paths

    @staticmethod
    def data_lines(path: Path) -> list[str]:
        return [line for line in path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#")]

    @staticmethod
    def metrics(path: Path) -> dict[str, str]:
        return {fields[1]: fields[2] for fields in
                (line.split("\t") for line in path.read_text().splitlines()[1:])}


class ClassifyVariantTests(unittest.TestCase):
    def test_existing_classifications(self) -> None:
        cases = (("A", "G", "snp"), ("AT", "GC", "mnp"),
                 ("A", "ATT", "indel"), ("ATT", "A", "indel"),
                 ("A", "*", "symbolic_or_star"),
                 ("A", "<DEL>", "symbolic_or_star"),
                 ("A", "A[chr2:100[", "symbolic_or_star"),
                 ("A", ".", "no_alt"))
        for ref, alt, expected in cases:
            with self.subTest(alt=alt):
                self.assertEqual(classifier.classify_variant(ref, alt), expected)

    def test_multiallelic_raises(self) -> None:
        with self.assertRaises(classifier.MalformedVcfError):
            classifier.classify_variant("A", "G,C")


class StreamingContractTests(unittest.TestCase, Harness):
    def tearDown(self) -> None:
        if hasattr(self, "directory"):
            self.directory.cleanup()

    def test_mixed_fixture_preserves_classification_header_and_filter_contracts(self) -> None:
        root = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(root))
        paths = root / "out.vcf", root / "account.tsv", root / "summary.txt"
        stats = classifier.stream_classify_vcf(
            FIXTURES_DIR / "classify_normal_mixed.vcf.gz", "cohort", *paths)
        self.assertEqual(stats.total_input_records, 5)
        self.assertEqual(stats.class_counts,
                         {"snp": 2, "mnp": 1, "indel": 1,
                          "symbolic_or_star": 1, "no_alt": 0})
        text = paths[0].read_text()
        self.assertNotIn("##FILTER=", text)
        self.assertIn("sample_a\tsample_b", text)
        self.assertTrue(all(line.split("\t")[6] == "." for line in self.data_lines(paths[0])))
        self.assertNotEqual(paths[0].read_bytes()[:2], b"\x1f\x8b")
        self.assertEqual(stat.S_IMODE(paths[0].stat().st_mode), 0o644)

    def test_same_locus_non_adjacent_duplicate_excludes_every_occurrence(self) -> None:
        code, paths = self.run_text(HEADER + _row("chr1", 100, "A", "G")
                                    + _row("chr1", 100, "A", "C")
                                    + _row("chr1", 100, "A", "G"))
        self.assertEqual(code, 0)
        rows = self.data_lines(paths[0])
        self.assertEqual([(r.split("\t")[3], r.split("\t")[4]) for r in rows], [("A", "C")])
        metrics = self.metrics(paths[1])
        self.assertEqual(metrics["duplicate_key_records"], "2")
        self.assertEqual(metrics["distinct_duplicate_keys"], "1")
        self.assertEqual(metrics["output_records"], "1")
        self.assertIn("chr1:100 A>G", paths[2].read_text())

    def test_locus_boundaries_preserve_output_order(self) -> None:
        text = HEADER + _row("chr1", 100, "A", "G") + _row("chr1", 100, "A", "C")
        text += _row("chr1", 101, "C", "T") + _row("chr1", 102, "G", "A")
        code, paths = self.run_text(text)
        self.assertEqual(code, 0)
        self.assertEqual([int(r.split("\t")[1]) for r in self.data_lines(paths[0])],
                         [100, 100, 101, 102])

    def test_empty_vcf_preserves_header_and_zero_accounting(self) -> None:
        code, paths = self.run_text(HEADER)
        self.assertEqual(code, 0)
        self.assertEqual(self.data_lines(paths[0]), [])
        self.assertIn("sample_a\tsample_b", paths[0].read_text())
        self.assertTrue(all(value == "0" for value in self.metrics(paths[1]).values()))

    def test_no_alt_is_excluded(self) -> None:
        code, paths = self.run_text(HEADER + _row("chr1", 1, "A", ".") + _row("chr1", 2))
        self.assertEqual(code, 0)
        self.assertEqual(len(self.data_lines(paths[0])), 1)
        self.assertEqual(self.metrics(paths[1])["no_alt_records"], "1")

    def test_unsorted_position_has_no_final_outputs(self) -> None:
        self._assert_atomic_failure(HEADER + _row("chr1", 100) + _row("chr1", 200)
                                    + _row("chr1", 100), "POS decreased")

    def test_contig_reentry_has_no_final_outputs(self) -> None:
        self._assert_atomic_failure(HEADER + _row("chr1", 100) + _row("chr2", 100)
                                    + _row("chr1", 200), "re-entered")

    def test_late_malformed_row_has_no_final_or_temporary_outputs(self) -> None:
        text = HEADER + "".join(_row("chr1", i) for i in range(1, 2001)) + "chr1\t2001\n"
        code, paths = self.run_text(text)
        self.assertNotEqual(code, 0)
        self.assertTrue(all(not path.exists() for path in paths))
        self.assertEqual(list(paths[0].parent.glob(".*.tmp")), [])

    def test_multiallelic_and_malformed_rows_fail(self) -> None:
        self._assert_atomic_failure(HEADER + _row("chr1", 1, alt="G,C"), "multiple ALT")
        self._assert_atomic_failure(HEADER + "chr1\t1\n", "expected at least 10")

    def _assert_atomic_failure(self, text: str, message: str) -> None:
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            code, paths = self.run_text(text)
        self.assertNotEqual(code, 0)
        self.assertIn(message, stderr.getvalue())
        self.assertTrue(all(not path.exists() for path in paths))


class StructuralMemoryRegressionTests(unittest.TestCase, Harness):
    def tearDown(self) -> None:
        if hasattr(self, "directory"):
            self.directory.cleanup()

    def test_large_synthetic_vcf_uses_streaming_production_api(self) -> None:
        # A meaningful record count exercises many locus flushes without a flaky RSS gate.
        code, paths = self.run_text(HEADER + "".join(_row("chr1", i) for i in range(1, 50001)))
        self.assertEqual(code, 0)
        self.assertEqual(self.metrics(paths[1])["output_records"], "50000")
        source = SCRIPT_PATH.read_text(encoding="utf-8")
        for legacy in ("candidates =", "output_records = tuple", '"\\n".join(lines)',
                       "sample_fields = tuple"):
            self.assertNotIn(legacy, source)
        self.assertIn('line.split("\\t", 9)', source)


if __name__ == "__main__":
    unittest.main()
