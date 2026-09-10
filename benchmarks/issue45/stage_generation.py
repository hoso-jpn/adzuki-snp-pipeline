#!/usr/bin/env python3
"""Stage the fixed Issue45 cohort against an unchanged production checkout.

Run only after all FASTQs have passed published checksum validation. This creates
an isolated Nextflow launch directory; it does not start or stop any workload.
"""

import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run-root", "production-checkout", "cohort-selection", "run-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    root, production = args.run_root.resolve(), args.production_checkout.resolve()
    output = args.run_dir.resolve()
    for path in (root, production, output):
        if not re.fullmatch(r"[A-Za-z0-9_./-]+", str(path)):
            raise ValueError("Launch paths must contain only shell-safe pathname characters")
    cohort = json.loads(args.cohort_selection.read_text())
    expected_sha = cohort["production_sha"]
    if (
        subprocess.check_output(
            ["git", "-C", str(production), "rev-parse", "HEAD"], text=True
        ).strip()
        != expected_sha
    ):
        raise ValueError("Production SHA differs from the fixed cohort")
    if subprocess.check_output(
        ["git", "-C", str(production), "status", "--porcelain"], text=True
    ).strip():
        raise ValueError("Production checkout is dirty")
    samples = cohort["samples"]
    if len(samples) != 51 or len({s["sample_accession"] for s in samples}) != 51:
        raise ValueError("Exactly 51 unique biological samples are required")
    inventory = json.loads((root / "audit/remote-inventory.json").read_text())
    completed = json.loads((root / "audit/fastq-checksums.json").read_text())
    expected_names = {
        f"{sample['run_accession']}_{mate}.fastq.gz" for sample in samples for mate in (1, 2)
    }
    if len(completed) != 102 or {row["name"] for row in completed} != expected_names:
        raise ValueError("The complete 102-file download verification is required")
    frozen_samples = []
    sheet = []
    for sample in sorted(samples, key=lambda row: row["run_accession"]):
        sid = sample["run_accession"]
        if (
            sample["study_accession"] != "PRJNA1138464"
            or sample["library_strategy"] != "WGS"
            or sample["library_layout"] != "PAIRED"
        ):
            raise ValueError("Unexpected project or assay")
        reads = []
        for mate in (1, 2):
            name = f"{sid}_{mate}.fastq.gz"
            path = root / "input" / name
            verified = json.loads((root / "audit" / f"{name}.verified.json").read_text())
            expected_md5 = sample["fastq_md5"].split(";")[mate - 1]
            if (
                verified["run_accession"] != sid
                or verified["biosample"] != sample["sample_accession"]
            ):
                raise ValueError("FASTQ verification sample identity mismatch")
            if verified["published_md5"] != expected_md5 or verified["local_md5"] != expected_md5:
                raise ValueError("Published FASTQ checksum mismatch")
            if (
                path.stat().st_size != int(sample["fastq_bytes"].split(";")[mate - 1])
                or sha256(path) != verified["sha256"]
            ):
                raise ValueError("FASTQ changed after download verification")
            if path.is_symlink():
                # Reused FASTQs were found in a user-trash archive. Preserve
                # the archive and freeze independent, checksum-identical copies
                # so this new run never depends on trash retention or symlinks
                # escaping Docker's automatically mounted input directory.
                if str(path.resolve(strict=True)) != inventory["fastqs"].get(name):
                    raise ValueError("FASTQ symlink differs from audited source")
                if shutil.disk_usage(root).free < 2_950_000_000_000:
                    raise ValueError("Storage gate failed before freezing reused FASTQ")
                temporary = path.with_name(path.name + f".freeze-{os.getpid()}.part")
                if temporary.exists():
                    raise ValueError("FASTQ copy staging path already exists")
                shutil.copyfile(path, temporary)
                if sha256(temporary) != verified["sha256"]:
                    raise ValueError("Reused FASTQ copy checksum mismatch")
                temporary.replace(path)
            reads.append(
                {
                    "name": name,
                    "published_md5": expected_md5,
                    "sha256": verified["sha256"],
                    "bytes": path.stat().st_size,
                }
            )
        # A run is the input read-group unit. No physical flowcell/lane is
        # inferred from the first read of an otherwise uninspected file.
        sheet.append(
            [
                sid,
                sid,
                str((root / "input" / reads[0]["name"]).resolve(strict=True)),
                str((root / "input" / reads[1]["name"]).resolve(strict=True)),
                sample["library_name"],
                "ILLUMINA",
                "",
            ]
        )
        frozen_samples.append(
            {
                "run_accession": sid,
                "biosample": sample["sample_accession"],
                "bioproject": "PRJNA1138464",
                "read_group_id": sid,
                "library_id": sample["library_name"],
                "platform": "ILLUMINA",
                "platform_unit": "",
                "fastqs": reads,
            }
        )
    references = {}
    for role, source in inventory["references"].items():
        path = Path(source)
        digest = sha256(path)
        expected = inventory["manifest"]["reference"][role]["checksum"].removeprefix("sha256:")
        if digest != expected:
            raise ValueError("Reference changed after source-run audit")
        references[role] = {"name": path.name, "sha256": digest}
    if shutil.disk_usage(root).free < 2_950_000_000_000:
        raise ValueError("Insufficient storage headroom before generation")
    output.mkdir(parents=True, exist_ok=False)
    reference_dir = output / "reference"
    reference_dir.mkdir()
    for role, source in inventory["references"].items():
        destination = reference_dir / references[role]["name"]
        shutil.copyfile(source, destination)
        if sha256(destination) != references[role]["sha256"]:
            raise ValueError("Staged reference checksum mismatch")
    (output / "bin").symlink_to(production / "bin", target_is_directory=True)
    template_path = Path(__file__).with_name("generate_gvcfs.nf.template")
    staged = output / "generate.nf"
    staged.write_text(template_path.read_text().replace("@PRODUCTION@", str(production)))
    with (output / "samples.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "sample_id",
                "read_group_id",
                "fastq_1",
                "fastq_2",
                "library_id",
                "platform",
                "platform_unit",
            )
        )
        writer.writerows(sheet)
    manifest = {
        "schema_version": 1,
        "cohort_id": output.name,
        "sample_count": 51,
        "production_sha": expected_sha,
        "production_checkout_clean": True,
        "wrapper_template_sha256": sha256(template_path),
        "preparation_helper_sha256": sha256(Path(__file__)),
        "execution_controller_sha256": sha256(Path(__file__).with_name("run_generation.py")),
        "staged_wrapper_sha256": sha256(staged),
        "samplesheet_sha256": sha256(output / "samples.csv"),
        "nextflow_version": "26.04.6",
        "gatk_version": "4.6.2.0",
        "gatk_container": "broadinstitute/gatk:4.6.2.0@sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb",
        "sample_ploidy": 2,
        "optical_duplicate_pixel_distance": 100,
        "haplotypecaller": {
            "emit_ref_confidence": "GVCF",
            "native_pair_hmm_threads": 4,
            "initial_java_heap_gib": 15,
            "create_output_variant_index": True,
        },
        "reference_accession": "GCF_016808095.1",
        "reference_id": "GCF_016808095.1",
        "reference": references,
        "samples": frozen_samples,
        "staged_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "benchmark_ready": False,
        "lineage_verified": False,
    }
    (output / "production_lineage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "nextflow.config").write_text(f"""includeConfig '{production}/nextflow.config'
params {{
    production_sha = '{expected_sha}'
    cohort_manifest = '{output}/production_lineage_manifest.json'
    input = '{output}/samples.csv'
    outdir = '{output}/results'
    reference_id = 'GCF_016808095.1'
    reference_fasta = '{reference_dir / references["fasta"]["name"]}'
    reference_fai = '{reference_dir / references["fai"]["name"]}'
    reference_dict = '{reference_dir / references["dict"]["name"]}'
}}
executor.cpus = 24
executor.memory = '96 GB'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g) --label adzuki.issue45_run={output.name}'
process {{
    withName: FASTP {{
        publishDir = [path: '${{params.outdir}}/qc/fastp', mode: 'copy', pattern: '*.{{json,html}}']
    }}
    withName: SAMTOOLS_INDEX {{ publishDir = [] }}
    withName: GATK_HAPLOTYPECALLER {{
        publishDir = [path: '${{params.outdir}}/variants/gvcf', mode: 'link', pattern: '*.g.vcf.gz*']
    }}
}}
trace {{
    enabled = true
    file = '{output}/trace.tsv'
    fields = 'task_id,hash,native_id,name,status,exit,attempt,cpus,memory,container,submit,duration,realtime,%cpu,peak_rss,peak_vmem,rchar,wchar,workdir'
}}
""")
    manifest["nextflow_config_sha256"] = sha256(output / "nextflow.config")
    (output / "production_lineage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
