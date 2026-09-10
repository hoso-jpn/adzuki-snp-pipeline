"""Compare isolated generation with the full workflow on synthetic input only."""

import argparse
import csv
import gzip
import itertools
import json
import os
import subprocess
import sys
from pathlib import Path


def execute(command, directory, name, success=True):
    with (directory / f"{name}.log").open("w") as log:
        result = subprocess.run(
            command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    if success and result.returncode != 0:
        raise RuntimeError((directory / f"{name}.log").read_text())
    if not success and result.returncode == 0:
        raise RuntimeError(f"Invalid input was accepted: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--production-checkout", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    production = args.production_checkout.resolve()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    os.environ["NXF_VER"] = "26.04.6"
    # Exercise the real stager's generated configuration, not just a hand-written
    # test config. An earlier single-quoted publishDir kept ${params.outdir}
    # literal and would have published real outputs into the wrong directory.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
    from test_issue45_generation_preparation import build_fixture, stage_fixture

    staged = stage_fixture(build_fixture(root / "synthetic-staging-contract"))
    flat = subprocess.check_output(
        ["nextflow", "config", str(staged), "-flat"], cwd=staged, text=True
    )
    for suffix in ("qc/fastp", "variants/gvcf"):
        if str(staged / "results" / suffix) not in flat:
            raise RuntimeError(
                f"Staged publishDir did not resolve to the expected location: {suffix}"
            )
    if "${params.outdir}" in flat:
        raise RuntimeError("Unexpanded publishDir in staged configuration")
    commit = subprocess.check_output(
        ["git", "-C", str(production), "rev-parse", "HEAD"], text=True
    ).strip()
    template = Path(__file__).resolve().parents[2] / "benchmarks/issue45/generate_gvcfs.nf.template"
    (root / "generate.nf").write_text(template.read_text().replace("@PRODUCTION@", str(production)))
    (root / "bin").symlink_to(production / "bin", target_is_directory=True)
    rows = list(csv.DictReader((production / "tests/data/samplesheet.csv").open()))
    for row in rows:
        for key in ("fastq_1", "fastq_2"):
            row[key] = str(production / row[key])
    with (root / "samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    expected = sorted({row["sample_id"] for row in rows})
    manifest = {
        "production_sha": commit,
        "reference_id": "synthetic-adzuki-v1",
        "sample_ploidy": 2,
        "optical_duplicate_pixel_distance": 100,
        "samples": [{"run_accession": sid, "biosample": "SYNTHETIC_" + sid} for sid in expected],
    }
    (root / "cohort.json").write_text(json.dumps(manifest))
    (root / "nextflow.config").write_text(f"""includeConfig '{production}/nextflow.config'
includeConfig '{production}/conf/test.config'
params {{
    input = '{root}/samples.csv'
    outdir = '{root}/generated'
    production_sha = '{commit}'
    cohort_manifest = '{root}/cohort.json'
    reference_id = 'synthetic-adzuki-v1'
    reference_fasta = '{production}/tests/data/reference/synthetic.fa'
    reference_fai = '{production}/tests/data/reference/synthetic.fa.fai'
    reference_dict = '{production}/tests/data/reference/prebuilt_dict/synthetic.dict'
}}
docker.enabled = true
""")
    execute(["nextflow", "lint", "generate.nf"], root, "lint")
    execute(
        [
            "nextflow",
            "run",
            "generate.nf",
            "--outdir",
            str(root / "generated"),
            "-w",
            str(root / "generation-work"),
        ],
        root,
        "generation",
    )
    full = root / "full"
    full.mkdir()
    execute(
        [
            "nextflow",
            "run",
            str(production),
            "-profile",
            "test,docker",
            "--input",
            str(root / "samples.csv"),
            "--outdir",
            str(full / "results"),
            "-w",
            str(full / "work"),
        ],
        full,
        "pipeline",
    )
    generated = root / "generated/variants/gvcf"
    if sorted(p.name for p in generated.glob("*.g.vcf.gz")) != [
        sid + ".g.vcf.gz" for sid in expected
    ]:
        raise RuntimeError("Unexpected generated sample set")
    for sid in expected:
        with (
            gzip.open(generated / f"{sid}.g.vcf.gz", "rt") as first,
            gzip.open(full / f"results/variants/gvcf/{sid}.g.vcf.gz", "rt") as second,
        ):
            left = (line for line in first if not line.startswith("#"))
            right = (line for line in second if not line.startswith("#"))
            for a, b in itertools.zip_longest(left, right):
                if a != b:
                    raise RuntimeError(f"gVCF record differs from full workflow: {sid}")
    # These invalid launches must stop before the first scientific task.
    execute(
        [
            "nextflow",
            "run",
            "generate.nf",
            "--reference_id",
            "",
            "-w",
            str(root / "bad-reference-work"),
        ],
        root,
        "missing-reference",
        success=False,
    )
    # Nextflow 26.04 coerces an empty CLI value to the string "true";
    # equality with the frozen identity must reject that value as well.
    if (
        "reference_id must match the frozen cohort manifest"
        not in (root / "missing-reference.log").read_text()
    ):
        raise RuntimeError("Missing reference failed for an unexpected reason")
    manifest["samples"][1]["biosample"] = manifest["samples"][0]["biosample"]
    (root / "cohort.json").write_text(json.dumps(manifest))
    execute(
        ["nextflow", "run", "generate.nf", "-w", str(root / "duplicate-biosample-work")],
        root,
        "duplicate-biosample",
        success=False,
    )
    if "BioSamples must be unique" not in (root / "duplicate-biosample.log").read_text():
        raise RuntimeError("Duplicate BioSample failed for an unexpected reason")
    for name in ("bad-reference-work", "duplicate-biosample-work"):
        if any((root / name).rglob(".command.run")):
            raise RuntimeError("Invalid input started a scientific task")
    print(
        "PASS: synthetic gVCF streams match full workflow; invalid lineage inputs fail before tasks"
    )


if __name__ == "__main__":
    main()
