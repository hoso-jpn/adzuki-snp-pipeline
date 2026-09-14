"""Issue #64: genotype-level DP/GQ quality masking for the GS panel.

Covers the shared policy evaluator (bin/gs_genotype_quality.py), the builder's
enabled path (bin/build_gs_panel.py), and the independent verifier
(bin/verify_gs_genotype_quality_mask.py).

Run with: python3 -m unittest discover -s tests/bin -v
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import json
import os
import random
import struct
import subprocess
import sys
import tempfile
import time
import types
import unittest
import zlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BIN_DIR = REPO_ROOT / "bin"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
QUALITY_FIXTURE = FIXTURES_DIR / "build_panel_genotype_quality.vcf.gz"
DISABLED_GOLDEN = FIXTURES_DIR / "build_panel_genotype_quality_disabled_golden"

# bin/ scripts import their siblings; see tests/bin/test_build_gs_panel_manifest.py.
sys.path.insert(0, str(BIN_DIR))


def _load_module(name: str, path: Path) -> types.ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


quality = _load_module("gs_genotype_quality", BIN_DIR / "gs_genotype_quality.py")
build_module = _load_module("build_gs_panel", BIN_DIR / "build_gs_panel.py")
verify_module = _load_module(
    "verify_gs_genotype_quality_mask", BIN_DIR / "verify_gs_genotype_quality_mask.py"
)

#: The policy the committed fixture's exact expectations are written for.
FIXTURE_POLICY = quality.GenotypeQualityPolicy(
    min_dp=10,
    min_gq=20,
    missing_format_field=quality.UNEVALUATED,
    missing_value=quality.MASK,
    malformed_value=quality.MASK,
)

OUTPUTS = {
    "matrix": "cohort.gs_panel.genotype_matrix.tsv.gz",
    "sample_metadata": "cohort.gs_panel.sample_metadata.tsv",
    "variant_metadata": "cohort.gs_panel.variant_metadata.tsv",
    "accounting": "cohort.gs_panel.genotype_encoding_accounting.tsv",
    "summary": "cohort.gs_panel.genotype_encoding_accounting.summary.txt",
    "masked_vcf": "cohort.gs_panel.quality_masked.vcf.gz",
    "policy": "cohort.gs_panel.genotype_quality_policy.json",
}


def _policy_argv(policy) -> list[str]:
    argv = ["--genotype-quality-mask"]
    if policy.min_dp is not None:
        argv += ["--genotype-min-dp", str(policy.min_dp)]
    if policy.min_gq is not None:
        argv += ["--genotype-min-gq", str(policy.min_gq)]
    return argv + [
        "--genotype-missing-format-field-policy",
        policy.missing_format_field,
        "--genotype-missing-value-policy",
        policy.missing_value,
        "--genotype-malformed-value-policy",
        policy.malformed_value,
    ]


def _build_argv(
    vcf: Path, directory: Path, policy=None, extra: list[str] | None = None
) -> list[str]:
    argv = [
        "--gs-pass-vcf",
        str(vcf),
        "--cohort-id",
        "cohort",
        "--sample-ploidy",
        "2",
        "--matrix-output",
        str(directory / OUTPUTS["matrix"]),
        "--sample-metadata-output",
        str(directory / OUTPUTS["sample_metadata"]),
        "--variant-metadata-output",
        str(directory / OUTPUTS["variant_metadata"]),
        "--genotype-accounting-output",
        str(directory / OUTPUTS["accounting"]),
        "--genotype-accounting-summary-output",
        str(directory / OUTPUTS["summary"]),
    ]
    if policy is not None:
        argv += _policy_argv(policy)
        argv += [
            "--quality-masked-vcf-output",
            str(directory / OUTPUTS["masked_vcf"]),
            "--genotype-quality-policy-output",
            str(directory / OUTPUTS["policy"]),
        ]
    return argv + (extra or [])


def _build(
    vcf: Path, directory: Path, policy=None, extra: list[str] | None = None
) -> tuple[int, str]:
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        code = build_module.main(_build_argv(vcf, directory, policy, extra))
    return code, stderr.getvalue()


def _verify(vcf: Path, directory: Path) -> tuple[int, str]:
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        code = verify_module.main(
            [
                "--cohort-id",
                "cohort",
                "--genotype-quality-policy",
                str(directory / OUTPUTS["policy"]),
                "--gs-pass-vcf",
                str(vcf),
                "--quality-masked-vcf",
                str(directory / OUTPUTS["masked_vcf"]),
                "--matrix",
                str(directory / OUTPUTS["matrix"]),
                "--variant-metadata",
                str(directory / OUTPUTS["variant_metadata"]),
                "--sample-metadata",
                str(directory / OUTPUTS["sample_metadata"]),
                "--genotype-accounting",
                str(directory / OUTPUTS["accounting"]),
                "--output",
                str(directory / "verification.tsv"),
                "--summary-output",
                str(directory / "verification.summary.txt"),
            ]
        )
    return code, stderr.getvalue()


def _accounting(directory: Path) -> dict[str, str]:
    lines = (directory / OUTPUTS["accounting"]).read_text(encoding="utf-8").splitlines()[1:]
    return {line.split("\t")[1]: line.split("\t")[2] for line in lines}


def _matrix_rows(directory: Path) -> list[list[str]]:
    with gzip.open(directory / OUTPUTS["matrix"], "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t") for line in handle]


def _vcf_records(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]


QUALITY_HEADER = (
    "##fileformat=VCFv4.2",
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "##contig=<ID=chrQ,length=100000000>",
)


def _small_vcf(path: Path, rows: list[tuple[str, list[str]]], samples: int = 2) -> Path:
    """Write a VCF from (FORMAT, sample fields) rows with AC/AN/AF defined."""
    names = [f"s{index}" for index in range(1, samples + 1)]
    lines = [
        *QUALITY_HEADER,
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *names]
        ),
    ]
    for index, (format_field, calls) in enumerate(rows):
        lines.append(
            "\t".join(
                [
                    "chrQ",
                    str(100 * (index + 1)),
                    ".",
                    "A",
                    "G",
                    "50",
                    "PASS",
                    "AC=1;AF=0.5;AN=2",
                    format_field,
                    *calls,
                ]
            )
        )
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))
    return path


# ==========================================================================
# Policy and evaluator
# ==========================================================================


class PolicyTests(unittest.TestCase):
    def test_an_enabled_policy_needs_a_threshold(self) -> None:
        with self.assertRaises(quality.InvalidPolicyError):
            quality.GenotypeQualityPolicy(min_dp=None, min_gq=None)

    def test_thresholds_must_be_non_negative_integers(self) -> None:
        for bad in (-1, 1.5, True, "10"):
            with self.subTest(value=bad), self.assertRaises(quality.InvalidPolicyError):
                quality.GenotypeQualityPolicy(min_dp=bad, min_gq=None)

    def test_actions_are_an_explicit_closed_set(self) -> None:
        with self.assertRaises(quality.InvalidPolicyError):
            quality.GenotypeQualityPolicy(min_dp=5, min_gq=None, missing_value="pass")

    def test_missing_field_handling_defaults_to_reject_not_pass(self) -> None:
        policy = quality.GenotypeQualityPolicy(min_dp=5, min_gq=None)
        self.assertEqual(
            (policy.missing_format_field, policy.missing_value, policy.malformed_value),
            ("reject", "reject", "reject"),
        )

    def test_document_round_trips_and_hash_is_canonical(self) -> None:
        document = FIXTURE_POLICY.document()
        self.assertEqual(
            quality.policy_from_document(json.loads(json.dumps(document))), FIXTURE_POLICY
        )
        reordered = dict(reversed(list(document.items())))
        self.assertEqual(quality.canonical_json_hash(reordered), FIXTURE_POLICY.policy_hash())
        self.assertNotEqual(
            FIXTURE_POLICY.policy_hash(),
            quality.GenotypeQualityPolicy(min_dp=11, min_gq=20).policy_hash(),
        )

    def test_a_tampered_or_disabled_document_is_not_an_enabled_policy(self) -> None:
        self.assertIsNone(quality.policy_from_document(quality.DISABLED_POLICY_DOCUMENT))
        tampered = FIXTURE_POLICY.document() | {"comparison": "value > threshold passes"}
        with self.assertRaises(quality.InvalidPolicyError):
            quality.policy_from_document(tampered)
        with self.assertRaises(quality.InvalidPolicyError):
            quality.policy_from_document({"schema": "something_else"})


class EvaluatorTests(unittest.TestCase):
    FORMAT = "GT:AD:DP:GQ"

    def verdict(self, sample: str, policy=FIXTURE_POLICY, format_field: str = FORMAT):
        return quality.RowQualityEvaluator(policy, format_field).evaluate(sample.split(":"))

    def test_threshold_itself_passes_and_one_below_masks(self) -> None:
        self.assertEqual(self.verdict("0/1:5,5:10:20").statuses, ("pass", "pass"))
        self.assertFalse(self.verdict("0/1:5,5:10:20").masked)
        self.assertFalse(self.verdict("0/1:5,5:11:21").masked)
        self.assertEqual(self.verdict("0/1:5,5:9:20").reasons, ("low", "none"))
        self.assertEqual(self.verdict("0/1:5,5:10:19").reasons, ("none", "low"))
        self.assertEqual(self.verdict("0/1:5,5:9:19").reasons, ("low", "low"))

    def test_phasing_does_not_change_the_verdict(self) -> None:
        for sample in ("0|1:5,5:9:40", "0|1:5,5:30:40"):
            with self.subTest(sample=sample):
                self.assertEqual(self.verdict(sample), self.verdict(sample.replace("|", "/")))

    def test_missing_values_are_never_zero_or_pass(self) -> None:
        missing = self.verdict("0/1:5,5:.:30")
        self.assertEqual(missing.statuses[0], "value_missing")
        truncated = self.verdict("0/1:5,5")
        self.assertEqual(truncated.statuses, ("value_truncated", "value_truncated"))
        # With missing_value=unevaluated the call is kept, but never as a pass.
        policy = quality.GenotypeQualityPolicy(min_dp=10, min_gq=None, missing_value="unevaluated")
        kept = self.verdict("0/1:5,5:.:30", policy)
        self.assertFalse(kept.masked)
        self.assertTrue(kept.unevaluated)
        self.assertEqual(kept.statuses, ("value_missing",))

    def test_absent_format_field_is_its_own_status(self) -> None:
        verdict = self.verdict("0/1:5,5:30", format_field="GT:AD:GQ")
        self.assertEqual(verdict.statuses, ("format_field_absent", "pass"))
        self.assertTrue(verdict.unevaluated)

    def test_malformed_values(self) -> None:
        for value in ("abc", "1.5", "-1", "", "1e3", "٣"):
            with self.subTest(value=value):
                verdict = self.verdict(f"0/1:5,5:{value}:30")
                self.assertEqual(verdict.statuses[0], "value_malformed")
                self.assertEqual(verdict.reasons[0], "malformed")

    def test_reject_names_the_field_status_and_option(self) -> None:
        policy = quality.GenotypeQualityPolicy(min_dp=10, min_gq=None)
        evaluator = quality.RowQualityEvaluator(policy, self.FORMAT)
        for sample, expected in (
            ("0/1:5,5:.:30", "DP is value_missing .*missing_value=reject"),
            ("0/1:5,5:abc:30", "DP is value_malformed .*malformed_value=reject"),
        ):
            with (
                self.subTest(sample=sample),
                self.assertRaisesRegex(quality.GenotypeQualityPolicyError, expected),
            ):
                evaluator.evaluate(sample.split(":"))
        with self.assertRaisesRegex(
            quality.GenotypeQualityPolicyError, "missing_format_field=reject"
        ):
            quality.RowQualityEvaluator(policy, "GT:GQ").evaluate(["0/1", "30"])

    def test_reason_metrics_partition_is_complete_and_excludes_none_none(self) -> None:
        metrics = quality.reason_metrics(FIXTURE_POLICY)
        self.assertEqual(len(metrics), 15)
        self.assertNotIn("quality_masked_reason.dp_none.gq_none", metrics)
        dp_only = quality.GenotypeQualityPolicy(min_dp=5, min_gq=None)
        self.assertEqual(
            quality.reason_metrics(dp_only),
            (
                "quality_masked_reason.dp_low",
                "quality_masked_reason.dp_unavailable",
                "quality_masked_reason.dp_malformed",
            ),
        )

    def test_masked_genotype_keeps_its_separator(self) -> None:
        self.assertEqual(quality.masked_genotype("0/1"), "./.")
        self.assertEqual(quality.masked_genotype("1|0"), ".|.")

    def test_allele_counts_and_info_rewrite(self) -> None:
        self.assertEqual(quality.allele_counts(["./.", "0|1", "0/.", "1/1"]), (3, 5))
        self.assertEqual(
            quality.rewrite_allele_info("AC=4;AF=0.500;AN=8;DP=48", 3, 5),
            "AC=3;AF=0.600000;AN=5;DP=48",
        )
        self.assertEqual(quality.rewrite_allele_info("AN=8;AC=4;AF=0.5", 0, 0), "AN=0;AC=0;AF=.")
        with self.assertRaisesRegex(ValueError, "lacks AF"):
            quality.rewrite_allele_info("AC=1;AN=2", 1, 2)


# ==========================================================================
# Builder: disabled compatibility
# ==========================================================================


class DisabledPolicyCompatibilityTests(unittest.TestCase):
    """With the mask off, outputs are byte-identical to main before Issue #64.

    The golden files were produced by `bin/build_gs_panel.py` at
    main@7ef04cb1dfe38fca6d7f100123f4fc7b8ddc0b6c over the DP/GQ fixture,
    before any Issue #64 change, so this pins compatibility against the
    implementation that was actually published rather than against itself.
    """

    GOLDEN = {
        "matrix": "genotype_matrix.tsv.gz",
        "sample_metadata": "sample_metadata.tsv",
        "variant_metadata": "variant_metadata.tsv",
        "accounting": "genotype_encoding_accounting.tsv",
        "summary": "genotype_encoding_accounting.summary.txt",
    }

    def test_disabled_outputs_match_the_pre_issue_golden_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            code, stderr = _build(QUALITY_FIXTURE, directory)
            self.assertEqual(code, 0, stderr)
            for key, golden in self.GOLDEN.items():
                with self.subTest(output=key):
                    self.assertEqual(
                        (directory / OUTPUTS[key]).read_bytes(),
                        (DISABLED_GOLDEN / golden).read_bytes(),
                    )
            self.assertFalse((directory / OUTPUTS["masked_vcf"]).exists())
            self.assertFalse((directory / OUTPUTS["policy"]).exists())

    def test_quality_options_without_the_mask_are_refused_not_ignored(self) -> None:
        for extra in (
            ["--genotype-min-dp", "10"],
            ["--genotype-min-gq", "20"],
            ["--genotype-missing-value-policy", "mask"],
            ["--quality-masked-vcf-output", "x.vcf.gz"],
        ):
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as tmp:
                code, stderr = _build(QUALITY_FIXTURE, Path(tmp), extra=extra)
                self.assertEqual(code, 1)
                self.assertIn("without --genotype-quality-mask", stderr)
                self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_the_mask_without_a_threshold_or_outputs_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            code, stderr = _build(
                QUALITY_FIXTURE,
                Path(tmp),
                extra=[
                    "--genotype-quality-mask",
                    "--quality-masked-vcf-output",
                    str(Path(tmp) / "q.vcf.gz"),
                    "--genotype-quality-policy-output",
                    str(Path(tmp) / "p.json"),
                ],
            )
            self.assertEqual(code, 1)
            self.assertIn("needs min_dp, min_gq, or both", stderr)
        with tempfile.TemporaryDirectory() as tmp:
            code, stderr = _build(
                QUALITY_FIXTURE,
                Path(tmp),
                extra=["--genotype-quality-mask", "--genotype-min-dp", "5"],
            )
            self.assertEqual(code, 1)
            self.assertIn("needs --quality-masked-vcf-output", stderr)


# ==========================================================================
# Builder: enabled
# ==========================================================================


class EnabledFixtureTests(unittest.TestCase):
    """Exact expectations for the committed fixture under FIXTURE_POLICY.

    Each fixture row targets one case the Issue lists; the comments name it.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls._tmp.name)
        code, stderr = _build(QUALITY_FIXTURE, cls.directory, FIXTURE_POLICY)
        assert code == 0, stderr
        cls.accounting = _accounting(cls.directory)
        cls.matrix = _matrix_rows(cls.directory)
        cls.masked = _vcf_records(cls.directory / OUTPUTS["masked_vcf"])
        cls.original = _vcf_records(QUALITY_FIXTURE)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_matrix(self) -> None:
        self.assertEqual(
            self.matrix,
            [
                ["variant_key", "s1", "s2", "s3", "s4"],
                # thresholds exactly (s1), threshold + 1 phased (s4): all encoded
                ["chrQ:100:A:G", "-1", "0", "1", "0"],
                # low DP only, low GQ only, both low, phased low DP: all masked
                ["chrQ:200:C:T", "nan", "nan", "nan", "nan"],
                # DP '.', GQ '.', ./. and 0/. originally missing
                ["chrQ:300:G:A", "nan", "nan", "nan", "nan"],
                # malformed DP, malformed GQ, haploid, non-biallelic index
                ["chrQ:400:T:C", "nan", "nan", "nan", "nan"],
                # DP+GQ truncated, GQ truncated, two passing calls
                ["chrQ:500:A:T", "nan", "nan", "-1", "0"],
                # FORMAT without DP: unevaluated for DP, s4 still fails GQ
                ["chrQ:600:C:G", "-1", "0", "1", "nan"],
                # FORMAT without GQ: s2/s3 fail DP, phased s4 kept
                ["chrQ:700:G:T", "-1", "nan", "nan", "0"],
                # GT at FORMAT index 2; phased hom-ref failing GQ, phased hom-alt passing
                ["chrQ:800:T:A", "nan", "1", "0", "-1"],
            ],
        )

    def test_accounting_partitions(self) -> None:
        a = {key: int(value) for key, value in self.accounting.items() if value.isdigit()}
        encoded = (
            a["standard_hom_ref_calls"] + a["standard_het_calls"] + a["standard_hom_alt_calls"]
        )
        self.assertEqual(a["total_genotype_cells"], 32)
        self.assertEqual(a["total_genotype_cells"], encoded + a["total_treated_as_missing"])
        self.assertEqual(
            a["total_treated_as_missing"],
            a["missing_calls"]
            + a["non_diploid_calls_treated_as_missing"]
            + a["non_biallelic_index_calls_treated_as_missing"]
            + a["quality_masked_calls"],
        )
        self.assertEqual(a["quality_evaluated_calls"], encoded + a["quality_masked_calls"])
        self.assertEqual(
            a["quality_evaluated_calls"],
            a["quality_passed_calls"]
            + a["quality_kept_unevaluated_calls"]
            + a["quality_masked_calls"],
        )
        self.assertEqual(
            a["quality_masked_calls"],
            sum(a[metric] for metric in quality.reason_metrics(FIXTURE_POLICY)),
        )
        for field in ("dp", "gq"):
            self.assertEqual(
                a["quality_evaluated_calls"],
                sum(a[f"{field}_status.{status}"] for status in quality.FIELD_STATUSES),
            )
        self.assertEqual(
            a["quality_masked_calls"],
            a["quality_masked_hom_ref_calls"]
            + a["quality_masked_het_calls"]
            + a["quality_masked_hom_alt_calls"],
        )

    def test_accounting_values(self) -> None:
        expected = {
            "standard_hom_ref_calls": "5",
            "standard_het_calls": "6",
            "standard_hom_alt_calls": "3",
            "missing_calls": "2",
            "non_diploid_calls_treated_as_missing": "1",
            "non_biallelic_index_calls_treated_as_missing": "1",
            "total_treated_as_missing": "18",
            "phased_genotype_count": "5",
            "genotype_quality_policy_hash": FIXTURE_POLICY.policy_hash(),
            "quality_evaluated_calls": "28",
            "quality_passed_calls": "9",
            "quality_kept_unevaluated_calls": "5",
            "quality_masked_calls": "14",
            "quality_masked_hom_ref_calls": "4",
            "quality_masked_het_calls": "7",
            "quality_masked_hom_alt_calls": "3",
            "quality_masked_reason.dp_low.gq_none": "4",
            "quality_masked_reason.dp_none.gq_low": "3",
            "quality_masked_reason.dp_low.gq_low": "1",
            "quality_masked_reason.dp_unavailable.gq_none": "1",
            "quality_masked_reason.dp_none.gq_unavailable": "2",
            "quality_masked_reason.dp_unavailable.gq_unavailable": "1",
            "quality_masked_reason.dp_malformed.gq_none": "1",
            "quality_masked_reason.dp_none.gq_malformed": "1",
            "dp_status.format_field_absent": "4",
            "dp_status.value_missing": "1",
            "dp_status.value_truncated": "1",
            "dp_status.value_malformed": "1",
            "gq_status.format_field_absent": "4",
            "gq_status.value_truncated": "2",
        }
        for metric, value in expected.items():
            with self.subTest(metric=metric):
                self.assertEqual(self.accounting[metric], value)
        # The historical rows keep their order at the top of the file.
        rows = list(self.accounting)
        self.assertEqual(rows[:9], [row[1] for row in _accounting_rows_of_golden()])

    def test_masked_vcf_changes_only_masked_gts_and_recomputed_info(self) -> None:
        self.assertEqual(len(self.masked), len(self.original))
        # Row 100 has nothing masked: byte-identical, INFO included.
        self.assertEqual(self.masked[0], self.original[0])
        row200 = self.masked[1]
        self.assertEqual(row200[7], "AC=0;AF=.;AN=0;DP=25;QD=18.00")
        self.assertEqual([call.split(":")[0] for call in row200[9:]], ["./.", "./.", "./.", ".|."])
        # FORMAT subfields other than GT survive masking, so the evidence stays auditable.
        self.assertEqual(row200[9], "./.:3,2:5:40:40,0,50")
        # Originally missing and non-standard calls are not touched.
        self.assertEqual(self.masked[2][11:], self.original[2][11:])
        self.assertEqual(self.masked[2][7], "AC=0;AF=0.000000;AN=1;DP=22;QD=15.00")
        self.assertEqual(self.masked[3][11:], self.original[3][11:])
        # Truncated sample fields stay truncated.
        self.assertEqual(self.masked[4][9:11], ["./.:5,5", "./.:0,20:20"])
        # GT at a non-zero FORMAT index is masked in place.
        self.assertEqual(self.masked[7][9], "20,0:20:.|.:10:0,10,200")
        for original, masked in zip(self.original, self.masked):
            self.assertEqual(original[:7], masked[:7])
            self.assertEqual(original[8], masked[8])

    def test_masked_vcf_is_bgzf_with_an_eof_block_and_a_provenance_line(self) -> None:
        raw = (self.directory / OUTPUTS["masked_vcf"]).read_bytes()
        self.assertTrue(
            raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
        )
        offset = 0
        while offset < len(raw):
            self.assertEqual(raw[offset : offset + 4], b"\x1f\x8b\x08\x04")
            self.assertEqual(raw[offset + 12 : offset + 14], b"BC")
            block_size = struct.unpack("<H", raw[offset + 16 : offset + 18])[0] + 1
            offset += block_size
        self.assertEqual(offset, len(raw))
        with gzip.open(self.directory / OUTPUTS["masked_vcf"], "rt") as handle:
            header = [line for line in handle if line.startswith("##gs_genotype_quality_mask=")]
        self.assertEqual(header, [quality.mask_header_line(FIXTURE_POLICY) + "\n"])

    def test_metadata_carry_the_masked_counts(self) -> None:
        sample_lines = (self.directory / OUTPUTS["sample_metadata"]).read_text().splitlines()
        self.assertEqual(sample_lines[0].split("\t")[-1], "quality_masked_genotype_count")
        self.assertEqual(
            [line.split("\t")[3::3] for line in sample_lines[1:]],
            [["5", "5"], ["5", "5"], ["4", "2"], ["4", "2"]],
        )
        variant_lines = (self.directory / OUTPUTS["variant_metadata"]).read_text().splitlines()
        self.assertEqual(
            [line.split("\t")[-1] for line in variant_lines[1:]],
            ["0", "4", "2", "2", "2", "1", "2", "1"],
        )

    def test_policy_document_is_published_and_deterministic(self) -> None:
        document = json.loads((self.directory / OUTPUTS["policy"]).read_text())
        self.assertEqual(quality.policy_from_document(document), FIXTURE_POLICY)
        with tempfile.TemporaryDirectory() as tmp:
            code, stderr = _build(QUALITY_FIXTURE, Path(tmp), FIXTURE_POLICY)
            self.assertEqual(code, 0, stderr)
            for name in OUTPUTS.values():
                self.assertEqual(
                    (Path(tmp) / name).read_bytes(), (self.directory / name).read_bytes(), name
                )

    def test_summary_names_it_a_sensitivity_not_accuracy(self) -> None:
        summary = (self.directory / OUTPUTS["summary"]).read_text()
        self.assertIn("masked to nan: 14", summary)
        self.assertIn("not a measure of genotype accuracy", summary)

    def test_the_verifier_accepts_the_build(self) -> None:
        code, stderr = _verify(QUALITY_FIXTURE, self.directory)
        self.assertEqual(code, 0, stderr)


def _accounting_rows_of_golden() -> list[list[str]]:
    lines = (DISABLED_GOLDEN / "genotype_encoding_accounting.tsv").read_text().splitlines()[1:]
    return [line.split("\t") for line in lines]


class RejectPolicyTests(unittest.TestCase):
    """With reject, a call that cannot be judged stops the build and publishes nothing."""

    def test_each_unjudgeable_case_is_rejected_with_its_location(self) -> None:
        policy = quality.GenotypeQualityPolicy(min_dp=10, min_gq=20)
        for label, row, expected in (
            (
                "dp value missing",
                ("GT:DP:GQ", ["0/1:.:30", "0/0:12:30"]),
                "sample 's1': DP is value_missing",
            ),
            (
                "gq malformed",
                ("GT:DP:GQ", ["0/1:12:30", "0/0:12:x"]),
                "sample 's2': GQ is value_malformed",
            ),
            ("dp not in format", ("GT:GQ", ["0/1:30", "0/0:30"]), "DP is format_field_absent"),
            ("gq truncated", ("GT:DP:GQ", ["0/1:12", "0/0:12:30"]), "GQ is value_truncated"),
        ):
            with self.subTest(case=label), tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                vcf = _small_vcf(
                    directory / "in.vcf.gz", [("GT:DP:GQ", ["0/0:30:60", "0/1:30:60"]), row]
                )
                out = directory / "out"
                out.mkdir()
                code, stderr = _build(vcf, out, policy)
                self.assertEqual(code, 1)
                self.assertIn("line 11: chrQ:200", stderr)
                self.assertIn(expected, stderr)
                self.assertEqual(list(out.iterdir()), [])

    def test_calls_already_missing_are_not_evaluated_even_under_reject(self) -> None:
        policy = quality.GenotypeQualityPolicy(min_dp=10, min_gq=20)
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _small_vcf(directory / "in.vcf.gz", [("GT:DP:GQ", ["./.", "0/.:.:x"])])
            code, stderr = _build(vcf, directory, policy)
            self.assertEqual(code, 0, stderr)
            self.assertEqual(_accounting(directory)["quality_evaluated_calls"], "0")
            self.assertEqual(_accounting(directory)["missing_calls"], "2")

    def test_the_header_must_define_the_recomputed_info_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            lines = [
                line
                for line in gzip.open(QUALITY_FIXTURE, "rt").read().splitlines()
                if not line.startswith("##INFO=<ID=AF,")
            ]
            vcf = directory / "in.vcf.gz"
            vcf.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
            out = directory / "out"
            out.mkdir()
            code, stderr = _build(vcf, out, FIXTURE_POLICY)
            self.assertEqual(code, 1)
            self.assertIn("does not define AF as INFO", stderr)


# ==========================================================================
# Streaming vs reference, and bounded memory, with the policy on
# ==========================================================================

_QUALITY_SHAPES = ("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1", "./.", "0/.", "0", "0/2")
_VALUES = ("0", "5", "9", "10", "11", "19", "20", "21", "60", ".", "x", "1.5")


def _generated_quality_lines(variants: int, samples: int, seed: int):
    rng = random.Random(seed)
    names = tuple(f"sample_{index}" for index in range(samples))
    yield from QUALITY_HEADER
    yield "\t".join(
        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *names]
    )
    formats = ("GT:AD:DP:GQ", "DP:GT:GQ", "GT:GQ", "GT:DP", "GQ:DP:GT")
    for variant in range(variants):
        format_field = formats[variant % len(formats)]
        keys = format_field.split(":")
        calls = []
        for _ in names:
            values = {
                "GT": rng.choice(_QUALITY_SHAPES),
                "AD": "3,4",
                "DP": rng.choice(_VALUES),
                "GQ": rng.choice(_VALUES),
            }
            fields = [values[key] for key in keys]
            if rng.random() < 0.1:
                fields = fields[: keys.index("GT") + 1]
            calls.append(":".join(fields))
        yield "\t".join(
            [
                "chrQ",
                str(variant * 11 + 1),
                ".",
                "A",
                "G",
                "50",
                "PASS",
                "AC=1;AF=0.5;AN=2",
                format_field,
                *calls,
            ]
        )


def _write_generated(path: Path, variants: int, samples: int, seed: int) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for line in _generated_quality_lines(variants, samples, seed):
            handle.write(line + "\n")
    return path


GENERATED_POLICY = quality.GenotypeQualityPolicy(
    min_dp=10,
    min_gq=20,
    missing_format_field="mask",
    missing_value="unevaluated",
    malformed_value="mask",
)


class EnabledStreamingEquivalenceTests(unittest.TestCase):
    def test_generated_cohorts_match_the_reference_and_verify(self) -> None:
        policies = (
            GENERATED_POLICY,
            quality.GenotypeQualityPolicy(
                min_dp=10,
                min_gq=None,
                missing_value="mask",
                malformed_value="unevaluated",
                missing_format_field="unevaluated",
            ),
            quality.GenotypeQualityPolicy(
                min_dp=None,
                min_gq=20,
                missing_value="mask",
                malformed_value="mask",
                missing_format_field="mask",
            ),
        )
        for (variants, samples, seed), policy in zip(
            ((40, 3, 1), (700, 9, 2), (3000, 17, 3)), policies
        ):
            with (
                self.subTest(variants=variants, samples=samples),
                tempfile.TemporaryDirectory() as tmp,
            ):
                directory = Path(tmp)
                vcf = _write_generated(directory / "in.vcf.gz", variants, samples, seed)
                code, stderr = _build(vcf, directory, policy)
                self.assertEqual(code, 0, stderr)

                reference = build_module.parse_gs_pass_vcf(vcf)
                with gzip.open(directory / OUTPUTS["matrix"], "rt") as handle:
                    self.assertEqual(
                        handle.read(),
                        "\n".join(
                            [
                                "\t".join(["variant_key", *reference.sample_names]),
                                *(
                                    "\t".join(row)
                                    for row in build_module.build_matrix_rows(reference, policy)
                                ),
                            ]
                        )
                        + "\n",
                    )
                for name, header, rows in (
                    (
                        OUTPUTS["sample_metadata"],
                        build_module.SAMPLE_METADATA_HEADER_WITH_QUALITY,
                        build_module.build_sample_metadata_rows("cohort", reference, policy),
                    ),
                    (
                        OUTPUTS["variant_metadata"],
                        build_module.VARIANT_METADATA_HEADER_WITH_QUALITY,
                        build_module.build_variant_metadata_rows("cohort", reference, policy),
                    ),
                    (
                        OUTPUTS["accounting"],
                        build_module.GENOTYPE_ACCOUNTING_HEADER,
                        build_module.build_genotype_accounting_rows("cohort", reference, policy),
                    ),
                ):
                    expected = (
                        "\n".join(["\t".join(header), *("\t".join(row) for row in rows)]) + "\n"
                    )
                    self.assertEqual((directory / name).read_text(), expected, name)
                self.assertEqual(
                    (directory / OUTPUTS["summary"]).read_text(),
                    build_module.build_genotype_accounting_summary_text(
                        "cohort", reference, policy
                    ),
                )
                self.assertGreater(int(_accounting(directory)["quality_masked_calls"]), 0)
                code, stderr = _verify(vcf, directory)
                self.assertEqual(code, 0, stderr)

    def test_the_enabled_production_path_never_materializes(self) -> None:
        helpers = (
            "parse_gs_pass_vcf",
            "build_matrix_rows",
            "build_sample_metadata_rows",
            "build_variant_metadata_rows",
            "build_genotype_accounting_rows",
            "build_genotype_accounting_summary_text",
            "write_matrix",
            "_resolve_record_cells",
        )
        originals = {name: getattr(build_module, name) for name in helpers}

        def forbidden(*_args: object, **_kwargs: object) -> None:
            raise AssertionError("the production path called a materializing helper")

        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            vcf = _write_generated(directory / "in.vcf.gz", 150, 5, 9)
            for name in helpers:
                setattr(build_module, name, forbidden)
            try:
                code, stderr = _build(vcf, directory, GENERATED_POLICY)
            finally:
                for name, original in originals.items():
                    setattr(build_module, name, original)
            self.assertEqual(code, 0, stderr)

    @unittest.skipUnless(Path("/proc/self/statm").exists(), "needs /proc to sample a child's RSS")
    def test_enabled_peak_memory_does_not_grow_with_the_variant_count(self) -> None:
        # The builder and the verifier both stream; hold samples fixed and
        # grow variants 16x. Anything retained per variant would show here.
        small, large = 2_000, 32_000
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            peaks: dict[str, dict[int, int]] = {"build": {}, "verify": {}}
            for variants in (small, large):
                run = directory / f"run_{variants}"
                run.mkdir()
                vcf = _write_generated(run / "in.vcf.gz", variants, 10, 21)
                peaks["build"][variants] = _peak_rss_kib(
                    [
                        sys.executable,
                        str(BIN_DIR / "build_gs_panel.py"),
                        *_build_argv(vcf, run, GENERATED_POLICY),
                    ]
                )
                peaks["verify"][variants] = _peak_rss_kib(
                    [
                        sys.executable,
                        str(BIN_DIR / "verify_gs_genotype_quality_mask.py"),
                        "--cohort-id",
                        "cohort",
                        "--genotype-quality-policy",
                        str(run / OUTPUTS["policy"]),
                        "--gs-pass-vcf",
                        str(vcf),
                        "--quality-masked-vcf",
                        str(run / OUTPUTS["masked_vcf"]),
                        "--matrix",
                        str(run / OUTPUTS["matrix"]),
                        "--variant-metadata",
                        str(run / OUTPUTS["variant_metadata"]),
                        "--sample-metadata",
                        str(run / OUTPUTS["sample_metadata"]),
                        "--genotype-accounting",
                        str(run / OUTPUTS["accounting"]),
                        "--output",
                        str(run / "v.tsv"),
                        "--summary-output",
                        str(run / "v.txt"),
                    ]
                )
        for tool, measured in peaks.items():
            with self.subTest(tool=tool):
                self.assertLess(
                    measured[large] - measured[small],
                    8 * 1024,
                    f"{tool} peak RSS grew {measured[small]} -> {measured[large]} KiB across a 16x variant increase",
                )


def _peak_rss_kib(command: list[str]) -> int:
    """Poll /proc/<pid>/statm; see test_build_gs_panel.py for why not ru_maxrss."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    statm = Path(f"/proc/{process.pid}/statm")
    page_kib = os.sysconf("SC_PAGE_SIZE") // 1024
    peak = 0
    while process.poll() is None:
        try:
            peak = max(peak, int(statm.read_text().split()[1]))
        except (OSError, IndexError, ValueError):
            break
        time.sleep(0.002)
    _out, err = process.communicate()
    if process.returncode != 0:
        raise AssertionError(err.decode(errors="replace"))
    return peak * page_kib


# ==========================================================================
# Verifier: every disagreement is caught
# ==========================================================================


class VerifierTamperTests(unittest.TestCase):
    """Each artifact is corrupted in one way; the verifier must refuse all of them."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self._tmp.name)
        code, stderr = _build(QUALITY_FIXTURE, self.directory, FIXTURE_POLICY)
        self.assertEqual(code, 0, stderr)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _rewrite_gz(self, name: str, old: str, new: str, bgzf: bool = False) -> None:
        path = self.directory / name
        with gzip.open(path, "rt") as handle:
            text = handle.read()
        self.assertIn(old, text)
        text = text.replace(old, new, 1)
        if bgzf:
            with build_module._StreamingBgzfWriter(path) as writer:
                writer.write(text)
        else:
            path.write_bytes(gzip.compress(text.encode(), mtime=0))

    def _rewrite_text(self, name: str, old: str, new: str) -> None:
        path = self.directory / name
        text = path.read_text()
        self.assertIn(old, text)
        path.write_text(text.replace(old, new, 1))

    def _assert_refused(self, expected: str) -> None:
        code, stderr = _verify(QUALITY_FIXTURE, self.directory)
        self.assertEqual(code, 1)
        self.assertIn(expected, stderr)
        self.assertFalse((self.directory / "verification.tsv").exists())

    def test_an_unmasked_cell_in_the_matrix(self) -> None:
        self._rewrite_gz(OUTPUTS["matrix"], "chrQ:200:C:T\tnan", "chrQ:200:C:T\t0")
        self._assert_refused("matrix token")

    def test_an_over_masked_cell_in_the_matrix(self) -> None:
        self._rewrite_gz(OUTPUTS["matrix"], "chrQ:100:A:G\t-1", "chrQ:100:A:G\tnan")
        self._assert_refused("matrix token")

    def test_a_masked_call_left_unmasked_in_the_vcf(self) -> None:
        self._rewrite_gz(OUTPUTS["masked_vcf"], "./.:3,2:5:40", "0/1:3,2:5:40", bgzf=True)
        self._assert_refused("masked VCF field differs")

    def test_a_retained_format_value_changed_in_the_vcf(self) -> None:
        self._rewrite_gz(OUTPUTS["masked_vcf"], "./.:3,2:5:40", "./.:3,2:5:41", bgzf=True)
        self._assert_refused("masked VCF field differs")

    def test_info_not_recomputed(self) -> None:
        self._rewrite_gz(OUTPUTS["masked_vcf"], "AC=0;AF=.;AN=0", "AC=4;AF=0.500;AN=8", bgzf=True)
        self._assert_refused("INFO is not the original")

    def test_info_changed_on_an_unmasked_row(self) -> None:
        self._rewrite_gz(
            OUTPUTS["masked_vcf"], "AC=4;AF=0.500;AN=8;DP=48", "AC=4;AF=0.5;AN=8;DP=48", bgzf=True
        )
        self._assert_refused("INFO changed on a row with no masked call")

    def test_provenance_line_for_another_policy(self) -> None:
        self._rewrite_gz(OUTPUTS["masked_vcf"], "MinDP=10", "MinDP=11", bgzf=True)
        self._assert_refused("header line")

    def test_accounting_reason_moved_between_pairs(self) -> None:
        self._rewrite_text(
            OUTPUTS["accounting"],
            "quality_masked_reason.dp_low.gq_none\t4",
            "quality_masked_reason.dp_low.gq_none\t3",
        )
        self._assert_refused("quality_masked_reason.dp_low.gq_none")

    def test_sample_metadata_masked_count(self) -> None:
        self._rewrite_text(
            OUTPUTS["sample_metadata"], "s3\t4\t0.500000\t1\t2", "s3\t4\t0.500000\t1\t1"
        )
        self._assert_refused("sample metadata row 2")

    def test_variant_metadata_masked_count(self) -> None:
        path = self.directory / OUTPUTS["variant_metadata"]
        lines = path.read_text().splitlines()
        fields = lines[2].split("\t")
        fields[-1] = "3"
        lines[2] = "\t".join(fields)
        path.write_text("\n".join(lines) + "\n")
        self._assert_refused("variant metadata row disagrees")

    def test_a_policy_document_that_differs_from_the_build(self) -> None:
        other = quality.GenotypeQualityPolicy(
            min_dp=5,
            min_gq=20,
            missing_format_field="unevaluated",
            missing_value="mask",
            malformed_value="mask",
        )
        (self.directory / OUTPUTS["policy"]).write_text(json.dumps(other.document()))
        self._assert_refused("")

    def test_a_truncated_masked_vcf(self) -> None:
        path = self.directory / OUTPUTS["masked_vcf"]
        with gzip.open(path, "rt") as handle:
            lines = handle.read().splitlines()
        with build_module._StreamingBgzfWriter(path) as writer:
            writer.write("\n".join(lines[:-1]) + "\n")
        self._assert_refused("ends before the original VCF")


class BgzfWriterTests(unittest.TestCase):
    def test_round_trips_across_blocks_including_incompressible_data(self) -> None:
        rng = random.Random(3)
        payload = (
            "".join(chr(rng.randint(33, 126)) for _ in range(200_000)) + "\n" + "A" * 150_000 + "\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.gz"
            with build_module._StreamingBgzfWriter(path) as writer:
                for start in range(0, len(payload), 7_777):
                    writer.write(payload[start : start + 7_777])
            raw = path.read_bytes()
            self.assertEqual(gzip.decompress(raw).decode(), payload)
            offset = blocks = 0
            while offset < len(raw):
                size = struct.unpack("<H", raw[offset + 16 : offset + 18])[0] + 1
                self.assertLessEqual(size, 0x10000)
                member = raw[offset : offset + size]
                self.assertEqual(
                    struct.unpack("<I", member[-8:-4])[0],
                    zlib.crc32(zlib.decompress(member[18:-8], -15)),
                )
                offset += size
                blocks += 1
            self.assertGreater(blocks, 5)


if __name__ == "__main__":
    unittest.main()


class CommittedModuleFixtureTests(unittest.TestCase):
    """The nf-test fixtures must be what the current builder produces, not a stale copy."""

    def test_masked_module_fixtures_match_a_fresh_build(self) -> None:
        fixtures = REPO_ROOT / "tests" / "modules" / "fixtures"
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            code, stderr = _build(fixtures / "gs_panel_quality.vcf.gz", directory, FIXTURE_POLICY)
            self.assertEqual(code, 0, stderr)
            for name in OUTPUTS.values():
                with self.subTest(output=name):
                    self.assertEqual(
                        (directory / name).read_bytes(),
                        (fixtures / "gs_panel_quality_masked" / name).read_bytes(),
                    )
