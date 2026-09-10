"""Synthetic preparation fixtures; these never authorize a real benchmark."""

import contextlib
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

HELPERS = Path(__file__).resolve().parents[2] / "benchmarks/issue45"
sys.path.insert(0, str(HELPERS))
import stage_generation as stage  # noqa: E402


def build_fixture(root):
    root.mkdir(parents=True, exist_ok=True)
    production = root / "production"
    production.mkdir()
    (production / "nextflow.config").write_text("params.outdir = 'results'\n")
    data = root / "data"
    (data / "audit").mkdir(parents=True)
    (data / "input").mkdir()
    samples, verified = [], []
    for index in range(51):
        sid, bio = f"synthetic_{index:03d}", f"SYNTHETIC_BIO_{index:03d}"
        md5s, sizes = [], []
        for mate in (1, 2):
            name = f"{sid}_{mate}.fastq.gz"
            payload = f"synthetic checksum fixture {sid} {mate}".encode()
            (data / "input" / name).write_bytes(payload)
            md5, sha = hashlib.md5(payload).hexdigest(), hashlib.sha256(payload).hexdigest()
            md5s.append(md5)
            sizes.append(str(len(payload)))
            record = {
                "name": name,
                "run_accession": sid,
                "biosample": bio,
                "published_md5": md5,
                "local_md5": md5,
                "sha256": sha,
            }
            (data / "audit" / f"{name}.verified.json").write_text(json.dumps(record))
            verified.append(record)
        samples.append(
            {
                "run_accession": sid,
                "sample_accession": bio,
                "study_accession": "PRJNA1138464",
                "library_strategy": "WGS",
                "library_layout": "PAIRED",
                "library_name": sid,
                "fastq_md5": ";".join(md5s),
                "fastq_bytes": ";".join(sizes),
            }
        )
    (data / "audit/fastq-checksums.json").write_text(json.dumps(verified))
    references, ref_manifest = {}, {}
    for role in ("fasta", "fai", "dict"):
        source = root / f"synthetic.{role}"
        source.write_text(f"synthetic {role}")
        references[role] = str(source)
        ref_manifest[role] = {"checksum": "sha256:" + stage.sha256(source)}
    (data / "audit/remote-inventory.json").write_text(
        json.dumps(
            {"references": references, "manifest": {"reference": ref_manifest}, "fastqs": {}}
        )
    )
    cohort = root / "synthetic-cohort.json"
    cohort.write_text(json.dumps({"production_sha": "a" * 40, "samples": samples}))
    return root


def stage_fixture(root):
    arguments = [
        "stage_generation.py",
        "--run-root",
        str(root / "data"),
        "--production-checkout",
        str(root / "production"),
        "--cohort-selection",
        str(root / "synthetic-cohort.json"),
        "--run-dir",
        str(root / "generated"),
    ]
    with (
        patch.object(sys, "argv", arguments),
        patch.object(stage.subprocess, "check_output", side_effect=["a" * 40 + "\n", ""]),
        patch.object(stage.shutil, "disk_usage", return_value=Mock(free=4_000_000_000_000)),
        contextlib.redirect_stdout(io.StringIO()),
    ):
        stage.main()
    return root / "generated"


class PreparationTests(unittest.TestCase):
    def test_changed_fastq_is_rejected_before_creating_run_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = build_fixture(Path(tmp))
            (root / "data/input/synthetic_000_1.fastq.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "FASTQ changed"):
                stage_fixture(root)
            self.assertFalse((root / "generated").exists())

    def test_partial_download_inventory_is_rejected_before_staging(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = build_fixture(Path(tmp))
            path = root / "data/audit/fastq-checksums.json"
            path.write_text(json.dumps(json.loads(path.read_text())[:-1]))
            with self.assertRaisesRegex(ValueError, "complete 102-file"):
                stage_fixture(root)
            self.assertFalse((root / "generated").exists())

    def test_archived_fastq_is_copied_without_changing_original_and_gate_stays_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = build_fixture(Path(tmp))
            link = root / "data/input/synthetic_000_1.fastq.gz"
            original_bytes = link.read_bytes()
            archived = root / "archived.fastq.gz"
            link.rename(archived)
            link.symlink_to(archived)
            inventory_path = root / "data/audit/remote-inventory.json"
            inventory = json.loads(inventory_path.read_text())
            inventory["fastqs"][link.name] = str(archived)
            inventory_path.write_text(json.dumps(inventory))
            output = stage_fixture(root)
            self.assertFalse(link.is_symlink())
            self.assertEqual(original_bytes, link.read_bytes())
            self.assertEqual(original_bytes, archived.read_bytes())
            manifest = json.loads((output / "production_lineage_manifest.json").read_text())
            self.assertFalse(manifest["benchmark_ready"])
            self.assertFalse(manifest["lineage_verified"])


if __name__ == "__main__":
    unittest.main()
