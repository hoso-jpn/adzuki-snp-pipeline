"""Unit tests for bin/build_run_manifest.py's schema v2 (DAG) mode.

Run with: python3 -m unittest discover -s tests/bin -v

Schema v2 is the document the Nextflow workflow builds for itself
(Issue #42). Its whole value is that it describes what a run *actually*
did, so these tests are mostly about the ways it must refuse to describe
something else: a container that two tasks disagreed about, accounting
that contradicts itself, a stale GS manifest from another cohort, an
unresolved commit, or any value that would publish where the run
happened.

tests/bin/test_build_run_manifest.py covers legacy v1 mode, which this
Issue deliberately left alone.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "bin" / "build_run_manifest.py"

# See tests/bin/test_build_run_manifest.py for why bin/ goes on sys.path
# here but needs no PYTHONPATH in production.
sys.path.insert(0, str(REPO_ROOT / "bin"))


def _load_module(name: str, path: Path) -> types.ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


manifest_module = _load_module("build_run_manifest_v2", SCRIPT_PATH)

FULL_SHA = "402f437b4a4f00023e0571bbe6d73dab3a83fadf"

PARAMETER_CLI_ARGS = [
    "--sample-ploidy", "2",
    "--genomicsdb-batch-size", "50",
    "--optical-duplicate-pixel-distance", "100",
    "--enable-gs-panel",
    "--snp-filter-qd-min", "2.0",
    "--snp-filter-qual-min", "30.0",
    "--snp-filter-sor-max", "3.0",
    "--snp-filter-fs-max", "60.0",
    "--snp-filter-mq-min", "40.0",
    "--snp-filter-mq-rank-sum-min", "-12.5",
    "--snp-filter-read-pos-rank-sum-min", "-8.0",
    "--indel-filter-qd-min", "2.0",
    "--indel-filter-qual-min", "30.0",
    "--indel-filter-fs-max", "200.0",
    "--indel-filter-read-pos-rank-sum-min", "-20.0",
]

REFERENCE_CLI_ARGS = [
    "--reference-id", "synthetic-adzuki-v1",
    "--reference-name", "Synthetic adzuki test reference",
    "--reference-species", "Vigna angularis",
    "--reference-cultivar", "synthetic",
    "--reference-accession", "",
]

# Two tasks of fastqc_raw reporting the same image, one multiqc task,
# and two independently-overridable aliases of the same GATK module.
RUNTIME_PROVENANCE_ROWS = [
    "fastqc_raw\tquay.io/biocontainers/fastqc:0.12.1@sha256:aa",
    "fastqc_raw\tquay.io/biocontainers/fastqc:0.12.1@sha256:aa",
    "fastqc_trimmed\tquay.io/biocontainers/fastqc:0.12.1@sha256:aa",
    "multiqc\tquay.io/biocontainers/multiqc:1.35@sha256:bb",
    "gatk_variantfiltration\tbroadinstitute/gatk:4.6.2.0@sha256:cc",
    "gatk_variantfiltration_gs\tbroadinstitute/gatk:4.6.2.0@sha256:dd",
    "build_run_manifest\tpython:3.12@sha256:ee",
]

INPUT_PROVENANCE_ROWS = [
    "00000002\tsample_b\tsample_b_L001\tlib_b\tILLUMINA\tpu_b\t"
    "b_R1.fastq.gz\tsha256:b1\tb_R2.fastq.gz\tsha256:b2",
    "00000000\tsample_a\tsample_a_L001\tlib_a\tILLUMINA\tpu_a1\t"
    "a1_R1.fastq.gz\tsha256:a11\ta1_R2.fastq.gz\tsha256:a12",
    "00000001\tsample_a\tsample_a_L002\tlib_a\tILLUMINA\tpu_a2\t"
    "a2_R1.fastq.gz\tsha256:a21\ta2_R2.fastq.gz\tsha256:a22",
]

REFERENCE_PROVENANCE_ROWS = [
    "fasta\tsynthetic.fa\tsha256:fa",
    "fai\tsynthetic.fa.fai\tsha256:fai",
    "dict\tsynthetic.dict\tsha256:dict",
    "bwa_index\tsynthetic.fa.pac\tsha256:i5",
    "bwa_index\tsynthetic.fa.0123\tsha256:i1",
    "bwa_index\tsynthetic.fa.amb\tsha256:i2",
    "bwa_index\tsynthetic.fa.ann\tsha256:i3",
    "bwa_index\tsynthetic.fa.bwt.2bit.64\tsha256:i4",
]

ARTIFACT_CHECKSUM_ROWS = [
    "sample_a.g.vcf.gz\tsha256:ga",
    "sample_b.g.vcf.gz\tsha256:gb",
    "cohort.raw.vcf.gz\tsha256:raw",
    "cohort.snp.pass.vcf.gz\tsha256:snp",
    "cohort.indel.pass.vcf.gz\tsha256:indel",
]

VARIANT_QC_ROWS = [
    "cohort_id\tstage\tvariant_type\tmetric\tvalue",
    "cohort\traw\tall\tnumber_of_samples\t2",
    "cohort\traw\tall\tsample_names\tsample_a,sample_b",
    "cohort\traw\tall\tcohort_total_genotypes\t4",
    "cohort\traw\tall\tcohort_missing_genotypes\t0",
    "cohort\traw\tall\tcohort_missingness_rate\t0.000000",
]

SAMPLE_QC_ROWS = [
    "\t".join(manifest_module.SAMPLE_QC_HEADER),
    "cohort\traw\tall\tsample_a\t1\t1\t0\t0\t0.000000\t8.0\t1",
    "cohort\traw\tall\tsample_b\t1\t1\t0\t0\t0.000000\t5.5\t1",
]

VARIANT_TYPE_ROWS = [
    "cohort_id\tmetric\tvalue",
    "cohort\traw_all_records\t2",
    "cohort\traw_snp_records\t2",
    "cohort\traw_indel_records\t0",
]

GS_PANEL_MANIFEST = {
    "schema_version": 2,
    "run_id": "20260901T000000Z-abcd1234",
    "generated_at": "2026-09-01T00:00:00Z",
    "cohort_id": "cohort",
    "panel_status": "empty",
    "manifest_hash": "sha256:" + "0" * 64,
    "containers": {"gs_normalize_variants": "bcftools:1.24"},
}


class _Fixture:
    """A complete, valid schema v2 invocation that each test perturbs."""

    def __init__(self, directory: Path) -> None:
        self.directory = directory
        self.runtime = self._write("runtime.tsv", RUNTIME_PROVENANCE_ROWS)
        self.inputs = [
            self._write(f"input_{index}.tsv", [row])
            for index, row in enumerate(INPUT_PROVENANCE_ROWS)
        ]
        self.reference = self._write("reference.tsv", REFERENCE_PROVENANCE_ROWS)
        self.artifacts = self._write("artifacts.tsv", ARTIFACT_CHECKSUM_ROWS)
        self.variant_qc = self._write("variant_qc.tsv", VARIANT_QC_ROWS)
        self.sample_qc = self._write("sample_qc.tsv", SAMPLE_QC_ROWS)
        self.variant_type = self._write("variant_type.tsv", VARIANT_TYPE_ROWS)
        self.gs_manifest = directory / "gs_panel.manifest.json"
        self.gs_manifest.write_text(json.dumps(GS_PANEL_MANIFEST), encoding="utf-8")
        self.output = directory / "cohort.run_manifest.json"

    def _write(self, name: str, rows: list[str]) -> Path:
        path = self.directory / name
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return path

    def argv(self, *, gs_panel: bool = True, git_commit: str = FULL_SHA) -> list[str]:
        argv = [
            "--mode", "dag-v2",
            "--cohort-id", "cohort",
            "--pipeline-version", "0.2.0",
            "--git-commit", git_commit,
            "--nextflow-version", "26.04.6",
            *REFERENCE_CLI_ARGS,
            *PARAMETER_CLI_ARGS,
            "--runtime-provenance", str(self.runtime),
            "--reference-provenance", str(self.reference),
            "--artifact-checksums", str(self.artifacts),
            "--variant-qc-tsv", str(self.variant_qc),
            "--sample-qc-tsv", str(self.sample_qc),
            "--variant-type-accounting-tsv", str(self.variant_type),
            "--output", str(self.output),
        ]
        for path in self.inputs:
            argv.extend(["--input-provenance", str(path)])
        argv.extend(
            ["--gs-panel-manifest", str(self.gs_manifest)] if gs_panel else ["--no-gs-panel"]
        )
        return argv


@contextlib.contextmanager
def fixture():
    with tempfile.TemporaryDirectory() as tmp:
        yield _Fixture(Path(tmp))


def run_main(argv: list[str]) -> tuple[int, str]:
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = manifest_module.main(argv)
    return exit_code, stderr.getvalue()


class DagManifestHappyPathTests(unittest.TestCase):
    def _manifest(self, **kwargs: object) -> dict[str, object]:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv(**kwargs))
            self.assertEqual(exit_code, 0, stderr)
            return json.loads(fake.output.read_text(encoding="utf-8"))

    def test_schema_version_is_two(self) -> None:
        self.assertEqual(self._manifest()["schema_version"], 2)

    def test_records_the_full_git_sha_and_runtime_nextflow_version(self) -> None:
        manifest = self._manifest()
        self.assertEqual(manifest["git_commit"], FULL_SHA)
        self.assertEqual(manifest["nextflow_version"], "26.04.6")

    def test_material_parameters_are_recorded_with_their_types(self) -> None:
        parameters = self._manifest()["parameters"]
        self.assertEqual(parameters["sample_ploidy"], 2)
        self.assertEqual(parameters["genomicsdb_batch_size"], 50)
        self.assertEqual(parameters["optical_duplicate_pixel_distance"], 100)
        self.assertIs(parameters["enable_gs_panel"], True)
        self.assertEqual(parameters["snp_filter_qd_min"], 2.0)
        self.assertEqual(parameters["indel_filter_read_pos_rank_sum_min"], -20.0)
        # Host-specific run parameters are not material provenance and
        # must not be published.
        self.assertNotIn("input", parameters)
        self.assertNotIn("outdir", parameters)

    def test_manifest_hash_covers_the_document_without_itself(self) -> None:
        manifest = self._manifest()
        without_hash = {k: v for k, v in manifest.items() if k != "manifest_hash"}
        self.assertEqual(
            manifest_module.canonical_json_hash(without_hash), manifest["manifest_hash"]
        )

    def test_reference_records_identity_and_every_file(self) -> None:
        reference = self._manifest()["reference"]
        self.assertEqual(reference["reference_id"], "synthetic-adzuki-v1")
        self.assertEqual(reference["reference_name"], "Synthetic adzuki test reference")
        self.assertEqual(reference["fasta"], {"filename": "synthetic.fa", "checksum": "sha256:fa"})
        self.assertEqual(
            [entry["filename"] for entry in reference["bwa_index"]],
            [
                "synthetic.fa.0123",
                "synthetic.fa.amb",
                "synthetic.fa.ann",
                "synthetic.fa.bwt.2bit.64",
                "synthetic.fa.pac",
            ],
        )

    def test_checksums_cover_the_historical_artifact_set(self) -> None:
        self.assertEqual(
            sorted(self._manifest()["checksums"]),
            [
                "cohort.indel.pass.vcf.gz",
                "cohort.raw.vcf.gz",
                "cohort.snp.pass.vcf.gz",
                "sample_a.g.vcf.gz",
                "sample_b.g.vcf.gz",
            ],
        )


class ContainerProvenanceTests(unittest.TestCase):
    def _containers(self, rows: list[str]) -> dict[str, str]:
        with fixture() as fake:
            fake.runtime.write_text("\n".join(rows) + "\n", encoding="utf-8")
            exit_code, stderr = run_main(fake.argv())
            self.assertEqual(exit_code, 0, stderr)
            manifest = json.loads(fake.output.read_text(encoding="utf-8"))
            return manifest["containers"]

    def _expect_failure(self, rows: list[str]) -> str:
        with fixture() as fake:
            fake.runtime.write_text("\n".join(rows) + "\n", encoding="utf-8")
            exit_code, stderr = run_main(fake.argv())
            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            return stderr

    def test_containers_are_keyed_by_process(self) -> None:
        containers = self._containers(RUNTIME_PROVENANCE_ROWS)
        self.assertEqual(
            sorted(containers),
            [
                "build_run_manifest",
                "fastqc_raw",
                "fastqc_trimmed",
                "gatk_variantfiltration",
                "gatk_variantfiltration_gs",
                "multiqc",
            ],
        )

    def test_many_tasks_of_one_process_collapse_to_one_identity(self) -> None:
        rows = [*RUNTIME_PROVENANCE_ROWS, "fastqc_raw\tquay.io/biocontainers/fastqc:0.12.1@sha256:aa"]
        self.assertEqual(
            self._containers(rows)["fastqc_raw"],
            "quay.io/biocontainers/fastqc:0.12.1@sha256:aa",
        )

    def test_aliases_of_one_module_keep_their_own_identities(self) -> None:
        containers = self._containers(RUNTIME_PROVENANCE_ROWS)
        self.assertEqual(containers["gatk_variantfiltration"], "broadinstitute/gatk:4.6.2.0@sha256:cc")
        self.assertEqual(containers["gatk_variantfiltration_gs"], "broadinstitute/gatk:4.6.2.0@sha256:dd")
        self.assertNotEqual(
            containers["gatk_variantfiltration"], containers["gatk_variantfiltration_gs"]
        )

    def test_one_process_with_two_identities_fails_rather_than_picking_one(self) -> None:
        stderr = self._expect_failure(
            [*RUNTIME_PROVENANCE_ROWS, "multiqc\tquay.io/biocontainers/multiqc:1.35@sha256:ff"]
        )
        self.assertIn("more than one container", stderr)
        self.assertIn("multiqc", stderr)

    def test_host_path_container_is_rejected(self) -> None:
        stderr = self._expect_failure([*RUNTIME_PROVENANCE_ROWS, "fastp\t/opt/images/fastp.sif"])
        self.assertIn("host filesystem path", stderr)

    def test_file_uri_container_is_rejected(self) -> None:
        stderr = self._expect_failure(
            [*RUNTIME_PROVENANCE_ROWS, "fastp\tfile:///opt/images/fastp.sif"]
        )
        self.assertIn("file://", stderr)

    def test_credential_bearing_container_is_rejected(self) -> None:
        stderr = self._expect_failure(
            [*RUNTIME_PROVENANCE_ROWS, "fastp\thttps://svc:s3cr3t@registry.internal/fastp:1.3.6"]
        )
        self.assertIn("credentials", stderr)

    def test_malformed_row_width_is_rejected(self) -> None:
        stderr = self._expect_failure([*RUNTIME_PROVENANCE_ROWS, "fastp"])
        self.assertIn("expected exactly 2", stderr)

    def test_non_canonical_process_key_is_rejected(self) -> None:
        stderr = self._expect_failure([*RUNTIME_PROVENANCE_ROWS, "FASTP\tquay.io/x/fastp:1"])
        self.assertIn("canonical process key", stderr)


class ReferenceProvenanceTests(unittest.TestCase):
    """Issue #42 review (P2): the reader rejects a truncated index bundle.

    The writer checks the files it can see; this side checks the shape of
    the TSV that reached it. A provenance file truncated in transit --
    partially written, partially staged -- would otherwise yield a
    manifest that silently understates the mapping input it claims to
    describe.
    """

    def _reference_rows(self, index_count: int) -> list[str]:
        index_rows = [
            f"bwa_index\tsynthetic.fa{suffix}\tsha256:i{number}"
            for number, suffix in enumerate(
                (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".sa")
            )
        ]
        return [*REFERENCE_PROVENANCE_ROWS[:3], *index_rows[:index_count]]

    def _run_with_index_count(self, index_count: int) -> tuple[int, str, bool]:
        with fixture() as fake:
            fake.reference.write_text(
                "\n".join(self._reference_rows(index_count)) + "\n", encoding="utf-8"
            )
            exit_code, stderr = run_main(fake.argv())
            return exit_code, stderr, fake.output.exists()

    def test_the_complete_five_file_index_is_accepted(self) -> None:
        exit_code, stderr, written = self._run_with_index_count(5)
        self.assertEqual(exit_code, 0, stderr)
        self.assertTrue(written)

    def test_no_index_files_is_rejected(self) -> None:
        exit_code, stderr, written = self._run_with_index_count(0)
        self.assertEqual(exit_code, 1)
        self.assertFalse(written)
        self.assertIn("exactly 5 files", stderr)

    def test_a_truncated_four_file_index_is_rejected(self) -> None:
        exit_code, stderr, written = self._run_with_index_count(4)
        self.assertEqual(exit_code, 1)
        self.assertFalse(written)
        self.assertIn("exactly 5 files", stderr)

    def test_a_six_file_index_is_rejected(self) -> None:
        exit_code, stderr, written = self._run_with_index_count(6)
        self.assertEqual(exit_code, 1)
        self.assertFalse(written)
        self.assertIn("exactly 5 files", stderr)

    def test_a_missing_single_file_role_is_rejected(self) -> None:
        with fixture() as fake:
            fake.reference.write_text(
                "\n".join(
                    row for row in REFERENCE_PROVENANCE_ROWS if not row.startswith("fai\t")
                )
                + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            self.assertIn("missing required reference role", stderr)


class InputProvenanceTests(unittest.TestCase):
    def test_read_groups_are_restored_to_samplesheet_order(self) -> None:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv())
            self.assertEqual(exit_code, 0, stderr)
            samples = json.loads(fake.output.read_text(encoding="utf-8"))["samples"]

        # The provenance files were written out of order on purpose.
        self.assertEqual(
            [entry["read_group_id"] for entry in samples],
            ["sample_a_L001", "sample_a_L002", "sample_b_L001"],
        )

    def test_samples_carry_metadata_and_filename_only_checksums(self) -> None:
        with fixture() as fake:
            run_main(fake.argv())
            samples = json.loads(fake.output.read_text(encoding="utf-8"))["samples"]

        self.assertEqual(
            samples[0],
            {
                "sample_id": "sample_a",
                "read_group_id": "sample_a_L001",
                "library_id": "lib_a",
                "platform": "ILLUMINA",
                "platform_unit": "pu_a1",
                "fastq_1": {"filename": "a1_R1.fastq.gz", "checksum": "sha256:a11"},
                "fastq_2": {"filename": "a1_R2.fastq.gz", "checksum": "sha256:a12"},
            },
        )
        # The samplesheet rank is a transport detail, never published.
        self.assertNotIn("rank", samples[0])

    def test_duplicate_samplesheet_rank_is_rejected(self) -> None:
        with fixture() as fake:
            fake.inputs[0].write_text(
                INPUT_PROVENANCE_ROWS[1].replace("sample_a_L001", "duplicate_rg") + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            self.assertIn("samplesheet position", stderr)

    def test_malformed_input_row_is_rejected(self) -> None:
        with fixture() as fake:
            fake.inputs[0].write_text("00000000\tsample_a\n", encoding="utf-8")
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("expected exactly 10", stderr)

    def test_non_integer_rank_is_rejected(self) -> None:
        with fixture() as fake:
            fake.inputs[0].write_text(
                INPUT_PROVENANCE_ROWS[0].replace("00000002", "second", 1) + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("rank is not an integer", stderr)


class AccountingTests(unittest.TestCase):
    def test_cohort_and_sample_accounting_are_recorded(self) -> None:
        with fixture() as fake:
            run_main(fake.argv())
            manifest = json.loads(fake.output.read_text(encoding="utf-8"))

        self.assertEqual(manifest["cohort_accounting"]["number_of_samples"], "2")
        self.assertEqual(
            [entry["sample"] for entry in manifest["sample_accounting"]],
            ["sample_a", "sample_b"],
        )
        self.assertEqual(
            manifest["sample_accounting"][0],
            {
                "sample": "sample_a",
                "reference_homozygous": "1",
                "non_reference_homozygous": "1",
                "heterozygous": "0",
                "missing": "0",
                "missingness_rate": "0.000000",
                "average_depth": "8.0",
                "singletons": "1",
            },
        )
        self.assertEqual(manifest["variant_type_accounting"]["raw_all_records"], "2")

    def test_sample_count_disagreement_fails(self) -> None:
        with fixture() as fake:
            fake.variant_qc.write_text(
                "\n".join(
                    row.replace("number_of_samples\t2", "number_of_samples\t3")
                    for row in VARIANT_QC_ROWS
                )
                + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            self.assertIn("reports 3 samples", stderr)

    def test_sample_name_order_disagreement_fails(self) -> None:
        with fixture() as fake:
            fake.variant_qc.write_text(
                "\n".join(
                    row.replace("sample_a,sample_b", "sample_b,sample_a")
                    for row in VARIANT_QC_ROWS
                )
                + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("do not match", stderr)

    def test_cohort_id_mismatch_in_sample_accounting_fails(self) -> None:
        with fixture() as fake:
            fake.sample_qc.write_text(
                "\n".join(SAMPLE_QC_ROWS).replace("cohort\traw", "other\traw") + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("is for cohort 'other'", stderr)

    def test_wrong_qc_stage_is_rejected(self) -> None:
        # Recording the filtered-SNP QC while claiming raw/all would be
        # undetectable downstream.
        with fixture() as fake:
            fake.variant_qc.write_text(
                "\n".join(VARIANT_QC_ROWS).replace("\traw\tall\t", "\tfiltered\tsnp\t")
                + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("must come from the raw/all", stderr)

    def test_unexpected_sample_qc_header_is_rejected(self) -> None:
        with fixture() as fake:
            rows = list(SAMPLE_QC_ROWS)
            rows[0] = rows[0].replace("singletons", "singleton_count")
            fake.sample_qc.write_text("\n".join(rows) + "\n", encoding="utf-8")
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("unexpected header", stderr)

    def test_missing_variant_type_metric_is_rejected(self) -> None:
        with fixture() as fake:
            fake.variant_type.write_text(
                "\n".join(row for row in VARIANT_TYPE_ROWS if "raw_snp_records" not in row)
                + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("missing required metric", stderr)


class GsPanelTests(unittest.TestCase):
    def test_enabled_run_embeds_the_gs_manifest_summary(self) -> None:
        with fixture() as fake:
            run_main(fake.argv(gs_panel=True))
            gs_panel = json.loads(fake.output.read_text(encoding="utf-8"))["gs_panel"]

        self.assertEqual(
            gs_panel,
            {
                "schema_version": 2,
                "run_id": GS_PANEL_MANIFEST["run_id"],
                "cohort_id": "cohort",
                "panel_status": "empty",
                "manifest_hash": GS_PANEL_MANIFEST["manifest_hash"],
            },
        )
        # A pointer, not a copy: the GS manifest's own containers and
        # checksums are not duplicated into the run manifest.
        self.assertNotIn("containers", gs_panel)

    def test_disabled_run_records_null(self) -> None:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv(gs_panel=False))
            self.assertEqual(exit_code, 0, stderr)
            self.assertIsNone(json.loads(fake.output.read_text(encoding="utf-8"))["gs_panel"])

    def test_gs_manifest_from_another_cohort_is_rejected(self) -> None:
        with fixture() as fake:
            fake.gs_manifest.write_text(
                json.dumps({**GS_PANEL_MANIFEST, "cohort_id": "other_cohort"}),
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            self.assertIn("other_cohort", stderr)

    def test_incomplete_gs_manifest_is_rejected(self) -> None:
        with fixture() as fake:
            payload = {k: v for k, v in GS_PANEL_MANIFEST.items() if k != "manifest_hash"}
            fake.gs_manifest.write_text(json.dumps(payload), encoding="utf-8")
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 1)
            self.assertIn("missing key", stderr)

    def test_gs_panel_selection_is_required(self) -> None:
        # Neither --gs-panel-manifest nor --no-gs-panel: argparse rejects
        # the invocation rather than the script guessing.
        with fixture() as fake:
            argv = [arg for arg in fake.argv(gs_panel=False) if arg != "--no-gs-panel"]
            with self.assertRaises(SystemExit):
                with contextlib.redirect_stderr(io.StringIO()):
                    manifest_module.main(argv)


class GitCommitTests(unittest.TestCase):
    def test_abbreviated_sha_is_rejected(self) -> None:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv(git_commit=FULL_SHA[:7]))

            self.assertEqual(exit_code, 1)
            self.assertFalse(fake.output.exists())
            self.assertIn("full 40-character", stderr)

    def test_empty_commit_is_rejected(self) -> None:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv(git_commit=""))

            self.assertEqual(exit_code, 1)
            self.assertIn("full 40-character", stderr)

    def test_non_hex_commit_is_rejected(self) -> None:
        with fixture() as fake:
            exit_code, _stderr = run_main(fake.argv(git_commit="z" * 40))
            self.assertEqual(exit_code, 1)


class PrivacyRegressionTests(unittest.TestCase):
    """Nothing host-specific may reach the serialized manifest."""

    def _assert_rejected(self, fake: _Fixture, needle: str) -> None:
        exit_code, stderr = run_main(fake.argv())
        self.assertEqual(exit_code, 1)
        self.assertFalse(fake.output.exists())
        self.assertIn(needle, stderr)

    def test_home_directory_path_in_a_filename_is_rejected(self) -> None:
        with fixture() as fake:
            fake.artifacts.write_text(
                "/home/alice/runs/cohort.raw.vcf.gz\tsha256:raw\n", encoding="utf-8"
            )
            self._assert_rejected(fake, "host filesystem path")

    def test_file_uri_anywhere_in_the_document_is_rejected(self) -> None:
        with fixture() as fake:
            fake.artifacts.write_text(
                "file:///opt/results/cohort.raw.vcf.gz\tsha256:raw\n", encoding="utf-8"
            )
            self._assert_rejected(fake, "file://")

    def test_credential_bearing_url_is_rejected(self) -> None:
        with fixture() as fake:
            fake.reference.write_text(
                "\n".join(
                    [
                        "fasta\thttps://user:hunter2@mirror.internal/synthetic.fa\tsha256:fa",
                        *REFERENCE_PROVENANCE_ROWS[1:],
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            self._assert_rejected(fake, "credentials")

    def test_private_network_location_is_rejected(self) -> None:
        with fixture() as fake:
            fake.runtime.write_text(
                "\n".join([*RUNTIME_PROVENANCE_ROWS, "fastp\t10.0.0.5:5000/fastp:1.3.6"])
                + "\n",
                encoding="utf-8",
            )
            self._assert_rejected(fake, "private or loopback")

    def test_loopback_location_is_rejected(self) -> None:
        with fixture() as fake:
            fake.runtime.write_text(
                "\n".join([*RUNTIME_PROVENANCE_ROWS, "fastp\tlocalhost:5000/fastp:1.3.6"])
                + "\n",
                encoding="utf-8",
            )
            self._assert_rejected(fake, "private or loopback")

    def test_a_clean_run_publishes_no_absolute_path_at_all(self) -> None:
        with fixture() as fake:
            exit_code, stderr = run_main(fake.argv())
            self.assertEqual(exit_code, 0, stderr)
            serialized = fake.output.read_text(encoding="utf-8")

            # Neither the temp directory the inputs came from nor this
            # checkout's own location appears anywhere in the document.
            self.assertNotIn(str(fake.directory), serialized)
            self.assertNotIn(str(REPO_ROOT), serialized)
            self.assertNotIn("file://", serialized)

    def test_sample_identifiers_that_look_like_names_are_not_rejected(self) -> None:
        # Privacy must not be enforced by guessing at personal names:
        # a cohort legitimately containing a sample called "alice" is
        # science, not a leak.
        with fixture() as fake:
            fake.sample_qc.write_text(
                "\n".join(SAMPLE_QC_ROWS).replace("sample_a", "alice").replace("sample_b", "bob")
                + "\n",
                encoding="utf-8",
            )
            fake.variant_qc.write_text(
                "\n".join(VARIANT_QC_ROWS).replace("sample_a,sample_b", "alice,bob") + "\n",
                encoding="utf-8",
            )
            exit_code, stderr = run_main(fake.argv())

            self.assertEqual(exit_code, 0, stderr)


class LegacyModeIsolationTests(unittest.TestCase):
    """v1 must stay reachable, and stay v1."""

    def test_default_mode_is_legacy(self) -> None:
        self.assertEqual(manifest_module._peek_mode([]), manifest_module.LEGACY_MODE)

    def test_dag_mode_flags_are_not_accepted_in_legacy_mode(self) -> None:
        # --runtime-provenance is a v2 concept; a legacy invocation that
        # passed it would otherwise silently ignore it.
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stderr(io.StringIO()):
                manifest_module.parse_args(
                    ["--mode", "legacy-v1", "--runtime-provenance", "x.tsv"]
                )

    def test_legacy_and_dag_schema_versions_differ(self) -> None:
        self.assertEqual(manifest_module.SCHEMA_VERSION, 1)
        self.assertEqual(manifest_module.DAG_SCHEMA_VERSION, 2)


if __name__ == "__main__":
    unittest.main()
