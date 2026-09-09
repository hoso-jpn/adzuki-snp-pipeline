"""Unit tests for bin/build_gs_panel.py.

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import os
import random
import subprocess
import sys
import tempfile
import time
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "build_gs_panel.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def _load_module(name: str, path: Path) -> types.ModuleType:
    """Load a bin/ script by path, without needing it to be a package."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build_module = _load_module("build_gs_panel", SCRIPT_PATH)


class ClassifyGenotypeTests(unittest.TestCase):
    """Tests for classify_genotype: the dosage table and every non-standard case."""

    def test_hom_ref_is_dosage_negative_one(self) -> None:
        cell = build_module.classify_genotype("0/0")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "-1")

    def test_het_is_dosage_zero(self) -> None:
        cell = build_module.classify_genotype("0/1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "0")

    def test_het_reversed_order_is_also_dosage_zero(self) -> None:
        cell = build_module.classify_genotype("1/0")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "0")

    def test_hom_alt_is_dosage_one(self) -> None:
        cell = build_module.classify_genotype("1/1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "1")

    def test_fully_missing_dot(self) -> None:
        cell = build_module.classify_genotype(".")
        self.assertEqual(cell.category, "missing")
        self.assertEqual(cell.dosage, "nan")

    def test_fully_missing_slash(self) -> None:
        cell = build_module.classify_genotype("./.")
        self.assertEqual(cell.category, "missing")

    def test_partially_missing_is_missing_not_standard(self) -> None:
        cell = build_module.classify_genotype("0/.")
        self.assertEqual(cell.category, "missing")
        self.assertEqual(cell.dosage, "nan")

    def test_phased_hom_ref_has_same_dosage_as_unphased(self) -> None:
        cell = build_module.classify_genotype("0|0")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "-1")
        self.assertTrue(cell.is_phased)

    def test_phased_het_has_same_dosage_as_unphased(self) -> None:
        # Phase records haplotype origin, not allele count: per the VCF
        # spec, 0|1 and 0/1 have identical allele content and must
        # resolve to the identical dosage.
        cell = build_module.classify_genotype("0|1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "0")
        self.assertTrue(cell.is_phased)

    def test_phased_hom_alt_has_same_dosage_as_unphased(self) -> None:
        cell = build_module.classify_genotype("1|1")
        self.assertEqual(cell.category, "standard")
        self.assertEqual(cell.dosage, "1")
        self.assertTrue(cell.is_phased)

    def test_unphased_calls_have_is_phased_false(self) -> None:
        self.assertFalse(build_module.classify_genotype("0/1").is_phased)

    def test_phased_missing_is_still_missing(self) -> None:
        cell = build_module.classify_genotype("0|.")
        self.assertEqual(cell.category, "missing")
        self.assertTrue(cell.is_phased)

    def test_haploid_is_non_diploid(self) -> None:
        cell = build_module.classify_genotype("0")
        self.assertEqual(cell.category, "non_diploid")

    def test_triploid_is_non_diploid(self) -> None:
        cell = build_module.classify_genotype("0/0/1")
        self.assertEqual(cell.category, "non_diploid")

    def test_phased_triploid_is_still_non_diploid(self) -> None:
        cell = build_module.classify_genotype("0|0|1")
        self.assertEqual(cell.category, "non_diploid")
        self.assertTrue(cell.is_phased)

    def test_non_biallelic_index_is_its_own_category(self) -> None:
        cell = build_module.classify_genotype("0/2")
        self.assertEqual(cell.category, "non_biallelic_index")


class ParseGsPassVcfTests(unittest.TestCase):
    """Tests for parse_gs_pass_vcf: sample names, records, and GT extraction."""

    def test_standard_fixture_reads_samples_and_records_in_order(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")

        self.assertEqual(vcf.sample_names, ("sample_a", "sample_b"))
        self.assertEqual(len(vcf.records), 3)
        self.assertEqual(vcf.records[0].variant_key, "chrPanel:100:A:G")
        self.assertEqual(vcf.records[0].sample_genotypes, ("0/0", "0/1"))

    def test_empty_vcf_has_zero_records_but_keeps_sample_names(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")

        self.assertEqual(vcf.sample_names, ("sample_a", "sample_b"))
        self.assertEqual(vcf.records, ())

    def test_missing_gt_format_raises_with_cause(self) -> None:
        with self.assertRaises(build_module.MalformedVcfError) as raised:
            build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_no_gt_format.vcf.gz")

        self.assertIn("FORMAT field has no GT subfield", str(raised.exception))

    def test_zero_sample_columns_raises_with_cause(self) -> None:
        with self.assertRaises(build_module.MalformedVcfError) as raised:
            build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_no_samples.vcf.gz")

        self.assertIn("expected at least 10", str(raised.exception))


class BuildMatrixTests(unittest.TestCase):
    """Tests for build_matrix_rows and write_matrix: on-disk shape and header."""

    def test_standard_fixture_matrix_rows(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_matrix_rows(vcf)

        self.assertEqual(
            rows,
            [
                ["chrPanel:100:A:G", "-1", "0"],
                ["chrPanel:200:C:T", "1", "-1"],
                ["chrPanel:300:G:A", "0", "0"],
            ],
        )

    def test_missing_fixture_matrix_rows(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_matrix_rows(vcf)

        self.assertEqual(
            rows,
            [
                ["chrPanel:100:A:G", "nan", "0"],
                ["chrPanel:200:C:T", "nan", "1"],
            ],
        )

    def test_matrix_header_is_written_even_with_zero_variants(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")

        with tempfile.TemporaryDirectory() as tmp:
            matrix_path = Path(tmp) / "matrix.tsv.gz"
            build_module.write_matrix(matrix_path, vcf)

            with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
                lines = handle.read().splitlines()

        self.assertEqual(lines, ["variant_key\tsample_a\tsample_b"])

    def test_identical_content_produces_byte_identical_compressed_output(self) -> None:
        # gzip embeds the current wall-clock time by default, which would
        # otherwise make two runs over the same input produce different
        # bytes (and therefore different checksums) despite identical
        # logical content -- verify the fix (mtime=0) actually holds.
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")

        with tempfile.TemporaryDirectory() as tmp:
            first_path = Path(tmp) / "first.tsv.gz"
            second_path = Path(tmp) / "second.tsv.gz"
            build_module.write_matrix(first_path, vcf)
            build_module.write_matrix(second_path, vcf)

            self.assertEqual(first_path.read_bytes(), second_path.read_bytes())

    def test_dosage_tokens_round_trip_through_float(self) -> None:
        for token in ("-1", "0", "1", "nan"):
            # every cell token must be float()-parseable without special-casing
            value = float(token)
            if token == "nan":
                self.assertNotEqual(value, value)  # NaN != NaN
            else:
                self.assertEqual(str(int(value)), token)


class SampleMetadataTests(unittest.TestCase):
    """Tests for build_sample_metadata_rows: per-sample missing/non-standard counts."""

    def test_standard_fixture_has_zero_missing(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "0", "0.000000", "0"])
        self.assertEqual(rows[1], ["cohort", "1", "sample_b", "0", "0.000000", "0"])

    def test_missing_fixture_counts_missing_per_sample(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        # sample_a: both variants missing (./. and 0/.) -> 2/2
        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "2", "1.000000", "0"])
        # sample_b: 0/2 variants missing -> 0/2
        self.assertEqual(rows[1], ["cohort", "1", "sample_b", "0", "0.000000", "0"])

    def test_non_standard_fixture_counts_non_standard_as_missing_too(self) -> None:
        # sample_a across the 4 variants: 0|1 (phased, but now a real
        # dosage -- not missing), 0 (haploid), 0/0/1 (triploid),
        # 0/2 (non-biallelic-index) -> 3 missing, all 3 non-standard.
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0], ["cohort", "0", "sample_a", "3", "0.750000", "3"])

    def test_empty_vcf_reports_na_rate_not_zero(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_empty.vcf.gz")
        rows = build_module.build_sample_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0][4], "NA")


class VariantMetadataTests(unittest.TestCase):
    """Tests for build_variant_metadata_rows: per-variant missing counts and order."""

    def test_standard_fixture_variant_order_and_keys(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        rows = build_module.build_variant_metadata_rows("cohort", vcf)

        self.assertEqual(
            [row[2] for row in rows],
            ["chrPanel:100:A:G", "chrPanel:200:C:T", "chrPanel:300:G:A"],
        )
        self.assertEqual([row[1] for row in rows], ["0", "1", "2"])

    def test_missing_fixture_missing_counts_per_variant(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_missing.vcf.gz")
        rows = build_module.build_variant_metadata_rows("cohort", vcf)

        self.assertEqual(rows[0][8], "1")  # chrPanel:100 -- sample_a missing
        self.assertEqual(rows[0][9], "0.500000")


class GenotypeAccountingTests(unittest.TestCase):
    """Tests for build_genotype_accounting_rows: reason-specific totals."""

    def test_standard_fixture_accounting(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_standard.vcf.gz")
        metrics = {
            row[1]: row[2] for row in build_module.build_genotype_accounting_rows("cohort", vcf)
        }

        self.assertEqual(metrics["total_genotype_cells"], "6")
        self.assertEqual(metrics["standard_hom_ref_calls"], "2")
        self.assertEqual(metrics["standard_het_calls"], "3")
        self.assertEqual(metrics["standard_hom_alt_calls"], "1")
        self.assertEqual(metrics["total_treated_as_missing"], "0")

    def test_non_standard_fixture_accounting_breaks_down_every_reason(self) -> None:
        # var1: sample_a=0|1 (phased, but a real dosage: standard het),
        #       sample_b=0/0 (standard hom-ref)
        # var2: sample_a=0, sample_b=1 (both haploid -> non_diploid)
        # var3: sample_a=0/0/1 (triploid -> non_diploid), sample_b=0/1 (standard het)
        # var4: sample_a=0/2 (non_biallelic_index), sample_b=0/0 (standard hom-ref)
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        metrics = {
            row[1]: row[2] for row in build_module.build_genotype_accounting_rows("cohort", vcf)
        }

        self.assertNotIn("phased_calls_treated_as_missing", metrics)
        self.assertEqual(metrics["non_diploid_calls_treated_as_missing"], "3")
        self.assertEqual(metrics["non_biallelic_index_calls_treated_as_missing"], "1")
        self.assertEqual(metrics["missing_calls"], "0")
        self.assertEqual(metrics["total_treated_as_missing"], "4")
        self.assertEqual(metrics["standard_hom_ref_calls"], "2")
        self.assertEqual(metrics["standard_het_calls"], "2")
        self.assertEqual(metrics["standard_hom_alt_calls"], "0")
        # exactly one call (var1 sample_a, 0|1) was phased
        self.assertEqual(metrics["phased_genotype_count"], "1")

    def test_phased_standard_call_is_not_double_counted_as_missing(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        metrics = {
            row[1]: row[2] for row in build_module.build_genotype_accounting_rows("cohort", vcf)
        }
        total_standard = (
            int(metrics["standard_hom_ref_calls"])
            + int(metrics["standard_het_calls"])
            + int(metrics["standard_hom_alt_calls"])
        )
        self.assertEqual(
            total_standard + int(metrics["total_treated_as_missing"]),
            int(metrics["total_genotype_cells"]),
        )

    def test_summary_text_mentions_every_reason(self) -> None:
        vcf = build_module.parse_gs_pass_vcf(FIXTURES_DIR / "build_panel_non_standard.vcf.gz")
        summary = build_module.build_genotype_accounting_summary_text("cohort", vcf)

        self.assertIn("non-diploid, treated as missing", summary)
        self.assertIn("non-biallelic-index, treated as missing", summary)
        self.assertIn("Phased genotype calls: 1", summary)
        self.assertNotIn("phased, treated as missing", summary)


class OutputContractTests(unittest.TestCase):
    """Tests for output headers."""

    def test_headers(self) -> None:
        self.assertEqual(
            build_module.GENOTYPE_ACCOUNTING_HEADER, ("cohort_id", "metric", "value")
        )
        self.assertEqual(
            build_module.SAMPLE_METADATA_HEADER,
            (
                "cohort_id",
                "sample_index",
                "sample_id",
                "missing_genotype_count",
                "missing_genotype_rate",
                "non_standard_genotype_count",
            ),
        )
        self.assertEqual(
            build_module.VARIANT_METADATA_HEADER,
            (
                "cohort_id",
                "variant_index",
                "variant_key",
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "missing_genotype_count",
                "missing_genotype_rate",
            ),
        )


class CliTests(unittest.TestCase):
    """Tests for the main() CLI entry point: success and failure exit codes."""

    def test_main_succeeds_and_writes_all_five_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "matrix.tsv.gz"
            sample_metadata_path = tmp_path / "sample_metadata.tsv"
            variant_metadata_path = tmp_path / "variant_metadata.tsv"
            accounting_path = tmp_path / "accounting.tsv"
            summary_path = tmp_path / "summary.txt"

            exit_code = build_module.main(
                [
                    "--gs-pass-vcf",
                    str(FIXTURES_DIR / "build_panel_standard.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--sample-ploidy",
                    "2",
                    "--matrix-output",
                    str(matrix_path),
                    "--sample-metadata-output",
                    str(sample_metadata_path),
                    "--variant-metadata-output",
                    str(variant_metadata_path),
                    "--genotype-accounting-output",
                    str(accounting_path),
                    "--genotype-accounting-summary-output",
                    str(summary_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            for path in (
                matrix_path,
                sample_metadata_path,
                variant_metadata_path,
                accounting_path,
                summary_path,
            ):
                self.assertTrue(path.exists(), f"{path} was not written")

    def test_main_fails_clearly_for_a_missing_input_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_vcf = tmp_path / "does-not-exist.vcf.gz"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = build_module.main(
                    [
                        "--gs-pass-vcf",
                        str(missing_vcf),
                        "--cohort-id",
                        "cohort",
                        "--sample-ploidy",
                        "2",
                        "--matrix-output",
                        str(tmp_path / "m.tsv.gz"),
                        "--sample-metadata-output",
                        str(tmp_path / "sm.tsv"),
                        "--variant-metadata-output",
                        str(tmp_path / "vm.tsv"),
                        "--genotype-accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--genotype-accounting-summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn(str(missing_vcf), stderr.getvalue())

    def test_main_fails_fast_for_haploid_ploidy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = build_module.main(
                    [
                        "--gs-pass-vcf",
                        str(FIXTURES_DIR / "build_panel_standard.vcf.gz"),
                        "--cohort-id",
                        "cohort",
                        "--sample-ploidy",
                        "1",
                        "--matrix-output",
                        str(tmp_path / "m.tsv.gz"),
                        "--sample-metadata-output",
                        str(tmp_path / "sm.tsv"),
                        "--variant-metadata-output",
                        str(tmp_path / "vm.tsv"),
                        "--genotype-accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--genotype-accounting-summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("diploid-only", stderr.getvalue())
            self.assertFalse((tmp_path / "m.tsv.gz").exists())

    def test_main_fails_fast_for_polyploid_ploidy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = build_module.main(
                    [
                        "--gs-pass-vcf",
                        str(FIXTURES_DIR / "build_panel_standard.vcf.gz"),
                        "--cohort-id",
                        "cohort",
                        "--sample-ploidy",
                        "4",
                        "--matrix-output",
                        str(tmp_path / "m.tsv.gz"),
                        "--sample-metadata-output",
                        str(tmp_path / "sm.tsv"),
                        "--variant-metadata-output",
                        str(tmp_path / "vm.tsv"),
                        "--genotype-accounting-output",
                        str(tmp_path / "a.tsv"),
                        "--genotype-accounting-summary-output",
                        str(tmp_path / "s.txt"),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("diploid-only", stderr.getvalue())

    def test_cli_subprocess_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            matrix_path = tmp_path / "matrix.tsv.gz"

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--gs-pass-vcf",
                    str(FIXTURES_DIR / "build_panel_empty.vcf.gz"),
                    "--cohort-id",
                    "cohort",
                    "--sample-ploidy",
                    "2",
                    "--matrix-output",
                    str(matrix_path),
                    "--sample-metadata-output",
                    str(tmp_path / "sm.tsv"),
                    "--variant-metadata-output",
                    str(tmp_path / "vm.tsv"),
                    "--genotype-accounting-output",
                    str(tmp_path / "a.tsv"),
                    "--genotype-accounting-summary-output",
                    str(tmp_path / "s.txt"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "variant_key\tsample_a\tsample_b\n")


# ==========================================================================
# Issue #44: the builder streams, and must go on producing exactly what the
# materializing implementation produced.
# ==========================================================================

MATRIX_OUTPUT_NAME = "matrix.tsv.gz"
PANEL_OUTPUT_NAMES = (
    MATRIX_OUTPUT_NAME,
    "sample_metadata.tsv",
    "variant_metadata.tsv",
    "genotype_accounting.tsv",
    "genotype_accounting_summary.txt",
)
PANEL_OUTPUT_FLAGS = (
    "--matrix-output",
    "--sample-metadata-output",
    "--variant-metadata-output",
    "--genotype-accounting-output",
    "--genotype-accounting-summary-output",
)

VCF_META_LINES = (
    "##fileformat=VCFv4.2",
    "##contig=<ID=chrPanel,length=100000000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
)


def _chrom_header(sample_names: tuple[str, ...]) -> str:
    fixed = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    return "\t".join([*fixed, *sample_names])


def _data_row(pos: int, genotype_fields: list[str], format_field: str = "GT") -> str:
    fixed = ["chrPanel", str(pos), ".", "A", "G", "100", "PASS", ".", format_field]
    return "\t".join([*fixed, *genotype_fields])


def _write_vcf(path: Path, lines: list[str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _panel_argv(vcf: Path, directory: Path, cohort_id: str = "cohort") -> list[str]:
    argv = [
        "--gs-pass-vcf", str(vcf),
        "--cohort-id", cohort_id,
        "--sample-ploidy", "2",
    ]
    for flag, name in zip(PANEL_OUTPUT_FLAGS, PANEL_OUTPUT_NAMES):
        argv.extend([flag, str(directory / name)])
    return argv


def _run_panel(vcf: Path, directory: Path, cohort_id: str = "cohort") -> tuple[int, str]:
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = build_module.main(_panel_argv(vcf, directory, cohort_id))
    return exit_code, stderr.getvalue()


#: Every genotype shape this encoding distinguishes, so generated
#: cohorts exercise all of them rather than only the clean cases.
_GENOTYPE_SHAPES = (
    "0/0", "0/1", "1/0", "1/1",          # standard, unphased
    "0|0", "0|1", "1|0", "1|1",          # standard, phased -- same dosages
    "./.", ".", "./1", "1/.",            # missing
    "0", "0/0/1", "0|1|1",               # non-diploid
    "0/2", "2/2", "1/3",                 # non-biallelic index
)


def _generated_vcf_lines(variants: int, samples: int, seed: int) -> list[str]:
    """A VCF exercising every genotype shape this encoding distinguishes."""
    rng = random.Random(seed)
    sample_names = tuple(f"sample_{index}" for index in range(samples))

    lines = [*VCF_META_LINES, _chrom_header(sample_names)]
    for variant in range(variants):
        calls = [rng.choice(_GENOTYPE_SHAPES) for _ in range(samples)]
        # A FORMAT that puts GT at a non-zero index must work identically.
        if variant % 3 == 1:
            lines.append(
                _data_row(
                    variant * 7 + 1,
                    [f"{rng.randint(1, 60)}:{call}" for call in calls],
                    format_field="DP:GT",
                )
            )
        else:
            lines.append(_data_row(variant * 7 + 1, calls))
    return lines


class StreamingEquivalenceTests(unittest.TestCase):
    """The streaming production path must reproduce the reference outputs exactly.

    `parse_gs_pass_vcf` and the `build_*_rows` family materialize the
    whole document and are no longer on the production path (Issue #44).
    They are kept as a declarative statement of what each output
    contains, which makes them an independent oracle: these tests assert
    the streaming CLI reproduces them byte for byte, so sample order,
    variant order, variant keys, dosage tokens, both metadata files and
    the accounting cannot drift without failing here.
    """

    def _assert_matches_reference(self, vcf: Path, directory: Path) -> None:
        exit_code, stderr = _run_panel(vcf, directory)
        self.assertEqual(exit_code, 0, stderr)

        reference = build_module.parse_gs_pass_vcf(vcf)

        with gzip.open(directory / MATRIX_OUTPUT_NAME, "rt", encoding="utf-8") as handle:
            streamed_matrix = handle.read()
        expected_matrix_rows = [
            "\t".join(["variant_key", *reference.sample_names]),
            *(
                "\t".join(row)
                for row in build_module.build_matrix_rows(reference)
            ),
        ]
        self.assertEqual(streamed_matrix, "\n".join(expected_matrix_rows) + "\n")

        expected_by_name = {
            "sample_metadata.tsv": (
                build_module.SAMPLE_METADATA_HEADER,
                build_module.build_sample_metadata_rows("cohort", reference),
            ),
            "variant_metadata.tsv": (
                build_module.VARIANT_METADATA_HEADER,
                build_module.build_variant_metadata_rows("cohort", reference),
            ),
            "genotype_accounting.tsv": (
                build_module.GENOTYPE_ACCOUNTING_HEADER,
                build_module.build_genotype_accounting_rows("cohort", reference),
            ),
        }
        for name, (header, rows) in expected_by_name.items():
            expected = "\n".join(
                ["\t".join(header), *("\t".join(row) for row in rows)]
            ) + "\n"
            self.assertEqual(
                (directory / name).read_text(encoding="utf-8"), expected, name
            )

        self.assertEqual(
            (directory / "genotype_accounting_summary.txt").read_text(encoding="utf-8"),
            build_module.build_genotype_accounting_summary_text("cohort", reference),
        )

    def test_every_committed_fixture_matches_the_reference_implementation(self) -> None:
        for name in (
            "build_panel_standard",
            "build_panel_missing",
            "build_panel_non_standard",
            "build_panel_empty",
        ):
            with self.subTest(fixture=name), tempfile.TemporaryDirectory() as tmp:
                self._assert_matches_reference(
                    FIXTURES_DIR / f"{name}.vcf.gz", Path(tmp)
                )

    def test_generated_cohorts_match_the_reference_implementation(self) -> None:
        # Wide enough to cover every genotype shape, several sample
        # counts, and enough variants to cross the matrix writer's
        # internal flush threshold more than once.
        for variants, samples, seed in ((1, 1, 1), (37, 3, 2), (500, 8, 3), (4000, 20, 4)):
            with self.subTest(variants=variants, samples=samples):
                with tempfile.TemporaryDirectory() as tmp:
                    directory = Path(tmp)
                    vcf = _write_vcf(
                        directory / "in.vcf.gz",
                        _generated_vcf_lines(variants, samples, seed),
                    )
                    self._assert_matches_reference(vcf, directory)

    def test_zero_variant_vcf_still_publishes_all_five_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "in.vcf.gz",
                [*VCF_META_LINES, _chrom_header(("sample_a", "sample_b"))],
            )
            exit_code, stderr = _run_panel(vcf, directory)

            self.assertEqual(exit_code, 0, stderr)
            for name in PANEL_OUTPUT_NAMES:
                self.assertTrue((directory / name).exists(), name)

            with gzip.open(directory / MATRIX_OUTPUT_NAME, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "variant_key\tsample_a\tsample_b\n")

            # A zero-variant panel reports NA rather than dividing by zero.
            sample_rows = (
                (directory / "sample_metadata.tsv").read_text(encoding="utf-8").splitlines()
            )
            self.assertEqual(sample_rows[1].split("\t")[4], "NA")

    def test_sample_and_variant_order_follow_the_file_not_a_sort(self) -> None:
        # Deliberately unsorted sample names and descending positions:
        # both must be preserved exactly as the VCF presented them.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            samples = ("zeta", "alpha", "mike")
            vcf = _write_vcf(
                directory / "in.vcf.gz",
                [
                    *VCF_META_LINES,
                    _chrom_header(samples),
                    _data_row(900, ["0/0", "0/1", "1/1"]),
                    _data_row(100, ["1/1", "0/0", "0/1"]),
                    _data_row(500, ["0/1", "1/1", "0/0"]),
                ],
            )
            exit_code, stderr = _run_panel(vcf, directory)
            self.assertEqual(exit_code, 0, stderr)

            with gzip.open(directory / MATRIX_OUTPUT_NAME, "rt", encoding="utf-8") as handle:
                matrix = handle.read().splitlines()
            self.assertEqual(matrix[0], "variant_key\tzeta\talpha\tmike")
            self.assertEqual(
                [row.split("\t")[0] for row in matrix[1:]],
                ["chrPanel:900:A:G", "chrPanel:100:A:G", "chrPanel:500:A:G"],
            )
            self.assertEqual(matrix[1], "chrPanel:900:A:G\t-1\t0\t1")

            sample_rows = (
                (directory / "sample_metadata.tsv").read_text(encoding="utf-8").splitlines()[1:]
            )
            self.assertEqual([row.split("\t")[2] for row in sample_rows], list(samples))

            variant_rows = (
                (directory / "variant_metadata.tsv").read_text(encoding="utf-8").splitlines()[1:]
            )
            self.assertEqual([row.split("\t")[1] for row in variant_rows], ["0", "1", "2"])
            self.assertEqual(
                [row.split("\t")[2] for row in variant_rows],
                ["chrPanel:900:A:G", "chrPanel:100:A:G", "chrPanel:500:A:G"],
            )

    def test_phased_calls_keep_their_dosages_and_are_counted_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "in.vcf.gz",
                [
                    *VCF_META_LINES,
                    _chrom_header(("sample_a", "sample_b")),
                    _data_row(100, ["0|0", "0|1"]),
                    _data_row(200, ["1|1", "1|0"]),
                ],
            )
            exit_code, stderr = _run_panel(vcf, directory)
            self.assertEqual(exit_code, 0, stderr)

            with gzip.open(directory / MATRIX_OUTPUT_NAME, "rt", encoding="utf-8") as handle:
                matrix = handle.read().splitlines()
            self.assertEqual(matrix[1], "chrPanel:100:A:G\t-1\t0")
            self.assertEqual(matrix[2], "chrPanel:200:A:G\t1\t0")

            accounting = dict(
                line.split("\t")[1:]
                for line in (directory / "genotype_accounting.tsv")
                .read_text(encoding="utf-8")
                .splitlines()[1:]
            )
            self.assertEqual(accounting["phased_genotype_count"], "4")
            # Phasing is informational: none of these are missing.
            self.assertEqual(accounting["total_treated_as_missing"], "0")
            self.assertEqual(accounting["missing_calls"], "0")


class MalformedVcfRejectionTests(unittest.TestCase):
    """Every shape mismatch must fail loudly, and publish nothing.

    Issue #44's review demonstrated the gap: a `#CHROM` header declaring
    two samples and a data row carrying one genotype was parsed without
    complaint into a matrix row one cell short of its own header. The
    old check only asked for "at least 10" fields, which cannot see it.
    """

    def _assert_rejected(self, lines: list[str], expected_message: str) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(directory / "in.vcf.gz", lines)
            exit_code, stderr = _run_panel(vcf, directory)

            self.assertEqual(exit_code, 1, stderr)
            self.assertIn(expected_message, stderr)
            # Nothing is published, and no staging debris survives.
            for name in PANEL_OUTPUT_NAMES:
                self.assertFalse((directory / name).exists(), f"published {name}")
            leftovers = sorted(
                path.name for path in directory.iterdir() if path.name.startswith(".")
            )
            self.assertEqual(leftovers, [], f"staging files left behind: {leftovers}")

    def test_sample_column_shortage_is_rejected(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["0/0"]),
            ],
            "expected exactly 11",
        )

    def test_sample_column_excess_is_rejected(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["0/0", "0/1", "1/1"]),
            ],
            "expected exactly 11",
        )

    def test_a_mismatch_only_on_the_final_row_is_still_rejected(self) -> None:
        # The failure must not depend on being near the start of the file:
        # a streaming reader that had already written 500 good rows must
        # still refuse to publish any of them.
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                *(_data_row(position, ["0/0", "0/1"]) for position in range(1, 501)),
                _data_row(501, ["0/0"]),
            ],
            "line 505",
        )

    def test_format_without_gt_is_rejected(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["30", "31"], format_field="DP"),
            ],
            "FORMAT field has no GT subfield",
        )

    def test_a_sample_field_missing_gts_subfield_is_rejected(self) -> None:
        # FORMAT declares DP:GT, so GT is at index 1 -- but sample_b
        # supplies only one subfield. This used to raise a bare
        # IndexError and escape as an unhandled traceback.
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["30:0/0", "31"], format_field="DP:GT"),
            ],
            "FORMAT subfield(s), but FORMAT places GT at index 1",
        )

    def test_a_sample_field_missing_gts_subfield_names_the_sample(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["30:0/0", "31"], format_field="DP:GT"),
            ],
            "sample 'sample_b'",
        )

    def test_data_row_before_the_chrom_header_is_rejected(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _data_row(100, ["0/0", "0/1"]),
                _chrom_header(("sample_a", "sample_b")),
            ],
            "data row seen before #CHROM header",
        )

    def test_malformed_fixed_column_count_is_rejected(self) -> None:
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                "chrPanel\t100\t.\tA\tG",
            ],
            "expected exactly 11",
        )

    def test_a_second_chrom_header_is_rejected(self) -> None:
        # It would silently redefine the sample count every later row is
        # checked against.
        self._assert_rejected(
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["0/0", "0/1"]),
                _chrom_header(("sample_a", "sample_b", "sample_c")),
            ],
            "a second #CHROM header line",
        )

    def test_a_vcf_with_no_chrom_header_is_rejected(self) -> None:
        self._assert_rejected(list(VCF_META_LINES), "no #CHROM header line found")

    def test_a_header_with_no_sample_columns_is_rejected(self) -> None:
        self._assert_rejected(
            [*VCF_META_LINES, _chrom_header(())], "#CHROM header has 9 fields"
        )

    def test_no_rejection_surfaces_as_a_bare_python_traceback(self) -> None:
        # Every case above must arrive as this script's own diagnosable
        # error, never as an unhandled IndexError/ValueError.
        for label, lines in {
            "short row": [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["0/0"]),
            ],
            "missing GT subfield": [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["30:0/0", "31"], format_field="DP:GT"),
            ],
        }.items():
            with self.subTest(case=label), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                vcf = _write_vcf(directory / "in.vcf.gz", lines)
                result = subprocess.run(
                    [sys.executable, str(SCRIPT_PATH), *_panel_argv(vcf, directory)],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 1)
                self.assertNotIn("Traceback", result.stderr)
                self.assertIn("build_gs_panel.py: error:", result.stderr)


class OutputPublicationTests(unittest.TestCase):
    """When a run fails, the previous panel must survive it intact."""

    def _write_existing_panel(self, directory: Path) -> dict[str, bytes]:
        vcf = _write_vcf(
            directory / "good.vcf.gz",
            [
                *VCF_META_LINES,
                _chrom_header(("sample_a", "sample_b")),
                _data_row(100, ["0/0", "0/1"]),
            ],
        )
        exit_code, stderr = _run_panel(vcf, directory)
        self.assertEqual(exit_code, 0, stderr)
        return {
            name: (directory / name).read_bytes() for name in PANEL_OUTPUT_NAMES
        }

    def test_a_failing_run_does_not_disturb_an_existing_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self._write_existing_panel(directory)

            bad = _write_vcf(
                directory / "bad.vcf.gz",
                [
                    *VCF_META_LINES,
                    _chrom_header(("sample_a", "sample_b")),
                    *(_data_row(position, ["0/0", "0/1"]) for position in range(1, 51)),
                    _data_row(51, ["0/0"]),
                ],
            )
            exit_code, _stderr = _run_panel(bad, directory)

            self.assertEqual(exit_code, 1)
            for name, content in before.items():
                self.assertEqual((directory / name).read_bytes(), content, name)

    def test_a_failure_partway_through_publishing_restores_every_output(self) -> None:
        # `os.replace` is atomic per file; five of them are not a
        # transaction. The publisher moves each existing output aside
        # before replacing it so that a failure midway can put the
        # previous panel back -- this pins that rollback.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            before = self._write_existing_panel(directory)

            replacement = _write_vcf(
                directory / "next.vcf.gz",
                [
                    *VCF_META_LINES,
                    _chrom_header(("sample_a", "sample_b")),
                    _data_row(700, ["1/1", "1/1"]),
                    _data_row(800, ["0/0", "0/0"]),
                ],
            )

            real_replace = build_module.os.replace
            calls = {"count": 0}

            def failing_replace(source: object, destination: object) -> None:
                # Fail once several replaces in, i.e. with some outputs
                # already swapped and others not.
                calls["count"] += 1
                if calls["count"] == 5:
                    raise OSError("simulated failure midway through publishing")
                real_replace(source, destination)

            build_module.os.replace = failing_replace
            try:
                with self.assertRaises(OSError):
                    _run_panel(replacement, directory)
            finally:
                build_module.os.replace = real_replace

            for name, content in before.items():
                self.assertEqual(
                    (directory / name).read_bytes(),
                    content,
                    f"{name} was not rolled back",
                )

    def test_publishing_leaves_no_staging_or_rollback_files_behind(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._write_existing_panel(directory)
            # A second successful run over the same paths exercises the
            # move-aside path and must clean up after itself.
            again = _write_vcf(
                directory / "again.vcf.gz",
                [
                    *VCF_META_LINES,
                    _chrom_header(("sample_a", "sample_b")),
                    _data_row(100, ["1/1", "1/1"]),
                ],
            )
            exit_code, stderr = _run_panel(again, directory)

            self.assertEqual(exit_code, 0, stderr)
            hidden = sorted(
                path.name for path in directory.iterdir() if path.name.startswith(".")
            )
            self.assertEqual(hidden, [])


class GzipContractTests(unittest.TestCase):
    """The compressed matrix must stay byte-for-byte what it always was.

    Issue #44's review established that swapping the one-shot
    `gzip.compress(payload, mtime=0)` for a streaming writer is not
    automatically byte-preserving: `gzip.GzipFile` writes `OS=255` where
    zlib's own gzip wrapper writes `OS=3`, so the decompressed bytes
    match while the file does not -- which would have invalidated every
    previously published matrix checksum.

    The streaming writer goes through `zlib.compressobj(wbits=31)`,
    which is that same zlib wrapper. These tests pin the equivalence
    against a live `gzip.compress` oracle rather than a hardcoded
    digest, so a future Python or zlib that changed it would fail here
    instead of silently changing published checksums.
    """

    def _build_matrix(self, directory: Path, lines: list[str]) -> Path:
        vcf = _write_vcf(directory / "in.vcf.gz", lines)
        exit_code, stderr = _run_panel(vcf, directory)
        self.assertEqual(exit_code, 0, stderr)
        return directory / MATRIX_OUTPUT_NAME

    def test_streamed_matrix_is_byte_identical_to_one_shot_gzip_compress(self) -> None:
        for variants, samples, seed in ((0, 2, 11), (3, 2, 12), (2500, 6, 13)):
            with self.subTest(variants=variants), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                matrix_path = self._build_matrix(
                    directory, _generated_vcf_lines(variants, samples, seed)
                )
                decompressed = gzip.decompress(matrix_path.read_bytes())

                self.assertEqual(
                    matrix_path.read_bytes(),
                    gzip.compress(decompressed, mtime=0),
                    "streamed gzip bytes diverged from the historical one-shot form",
                )

    def test_repeated_runs_over_one_input_produce_identical_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            lines = _generated_vcf_lines(400, 5, 14)
            vcf = _write_vcf(directory / "in.vcf.gz", lines)

            digests = set()
            for run in range(3):
                run_directory = directory / f"run{run}"
                run_directory.mkdir()
                exit_code, stderr = _run_panel(vcf, run_directory)
                self.assertEqual(exit_code, 0, stderr)
                digests.add((run_directory / MATRIX_OUTPUT_NAME).read_bytes())

            self.assertEqual(len(digests), 1, "matrix bytes were not deterministic")

    def test_the_gzip_header_carries_no_filename_and_no_host_path(self) -> None:
        # A gzip member may embed the source filename (FLG bit 3). The
        # builder writes into a hidden staging file whose name would
        # otherwise be recorded in the published artifact, and the
        # staging path is an absolute host path.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            matrix_path = self._build_matrix(
                directory, _generated_vcf_lines(50, 3, 15)
            )
            raw = matrix_path.read_bytes()

            self.assertEqual(raw[:2], b"\x1f\x8b")
            flags = raw[3]
            self.assertEqual(flags & 0x08, 0, "gzip FNAME flag is set")
            self.assertEqual(flags & 0x10, 0, "gzip FCOMMENT flag is set")
            # mtime must stay zeroed for reproducibility.
            self.assertEqual(raw[4:8], b"\x00\x00\x00\x00")

            self.assertNotIn(b".partial", raw[:512])
            self.assertNotIn(str(directory).encode(), raw)
            self.assertNotIn(b"matrix.tsv.gz", raw[:512])


class ProductionPathIsBoundedTests(unittest.TestCase):
    """The CLI must not reach for the materializing reference implementation."""

    MATERIALIZING_HELPERS = (
        "parse_gs_pass_vcf",
        "build_matrix_rows",
        "build_sample_metadata_rows",
        "build_variant_metadata_rows",
        "build_genotype_accounting_rows",
        "build_genotype_accounting_summary_text",
        "write_matrix",
    )

    def test_main_never_calls_a_whole_document_helper(self) -> None:
        # Structural, not statistical: if `main` ever routes back through
        # one of these, memory becomes O(variants x samples) again and
        # this fails immediately rather than at real-cohort scale.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_vcf(
                directory / "in.vcf.gz", _generated_vcf_lines(120, 4, 16)
            )

            originals = {
                name: getattr(build_module, name) for name in self.MATERIALIZING_HELPERS
            }

            def forbidden(name: str):
                def _raise(*args: object, **kwargs: object) -> None:
                    raise AssertionError(
                        f"the production path called the materializing helper {name}()"
                    )

                return _raise

            for name in self.MATERIALIZING_HELPERS:
                setattr(build_module, name, forbidden(name))
            try:
                exit_code, stderr = _run_panel(vcf, directory)
            finally:
                for name, original in originals.items():
                    setattr(build_module, name, original)

            self.assertEqual(exit_code, 0, stderr)
            for name in PANEL_OUTPUT_NAMES:
                self.assertTrue((directory / name).exists(), name)

    @unittest.skipUnless(
        Path("/proc/self/statm").exists(), "needs /proc to sample a child's RSS"
    )
    def test_peak_memory_does_not_grow_with_the_variant_count(self) -> None:
        # "The fixture is small so the memory was small" proves nothing.
        # This runs the real CLI over two cohorts 16x apart in variant
        # count, with the sample count held fixed, and compares the peak
        # RSS each run actually reached. The previous implementation
        # retained every record, every matrix row, the joined text and
        # its encoded bytes, so it grew close to linearly across this
        # range; a bounded implementation should barely move.
        small, large = 2_000, 32_000
        samples = 10

        measurements = {}
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for variants in (small, large):
                vcf = _write_generated_vcf(
                    directory / f"in_{variants}.vcf.gz", variants, samples, 17
                )
                run_directory = directory / f"out_{variants}"
                run_directory.mkdir()
                measurements[variants] = _peak_rss_kib_for_panel_run(
                    vcf, run_directory
                )

        growth = measurements[large] - measurements[small]
        self.assertLess(
            growth,
            8 * 1024,
            "peak RSS grew by "
            f"{growth} KiB between {small} and {large} variants "
            f"({measurements[small]} -> {measurements[large]} KiB), "
            "which suggests something is being retained per variant",
        )


def _peak_rss_kib_for_panel_run(vcf: Path, directory: Path) -> int:
    """Run the builder in a fresh interpreter and report its peak RSS.

    Sampled from `/proc/<pid>/statm` rather than taken from the child's
    own `ru_maxrss`, which would be the wrong number here: on Linux
    `ru_maxrss` is carried across `fork`/`exec`, so a child spawned from
    a test runner that had already grown large reports the *runner's*
    high-water mark and every measurement collapses to the same value.
    Current RSS is not inherited, so polling it measures this run.

    Polling can in principle miss a spike between samples. It does not
    hide the regression this guards against: an implementation that
    retains per-variant state holds it until the end of the run, so its
    peak is a plateau, not a spike.
    """
    process = subprocess.Popen(
        [sys.executable, str(SCRIPT_PATH), *_panel_argv(vcf, directory)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    statm = Path(f"/proc/{process.pid}/statm")
    page_kib = os.sysconf("SC_PAGE_SIZE") // 1024
    peak_pages = 0
    while process.poll() is None:
        try:
            peak_pages = max(peak_pages, int(statm.read_text().split()[1]))
        except (OSError, IndexError, ValueError):
            break
        time.sleep(0.002)

    _stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise AssertionError(f"builder failed: {stderr.decode(errors='replace')}")

    return peak_pages * page_kib


def _write_generated_vcf(path: Path, variants: int, samples: int, seed: int) -> Path:
    """Write a generated VCF straight to disk, never holding it in memory.

    The test process must stay small: it is the parent of the measured
    run, and materializing a large VCF here would dominate the machine's
    memory picture even though it has nothing to do with the builder.
    """
    rng = random.Random(seed)
    sample_names = tuple(f"sample_{index}" for index in range(samples))
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for line in VCF_META_LINES:
            handle.write(line + "\n")
        handle.write(_chrom_header(sample_names) + "\n")
        for variant in range(variants):
            calls = [rng.choice(_GENOTYPE_SHAPES) for _ in range(samples)]
            handle.write(_data_row(variant * 7 + 1, calls) + "\n")
    return path


if __name__ == "__main__":
    unittest.main()
