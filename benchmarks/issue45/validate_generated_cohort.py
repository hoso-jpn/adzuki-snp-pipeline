#!/usr/bin/env python3
"""Verify every generated public gVCF against its frozen run and usable index.

Run after generation has finished. The detailed manifest contains local paths and
must remain outside Git; publish only its explicitly sanitized evidence member.
"""

import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

from stage_generation import sha256

BCFTOOLS = (
    "quay.io/biocontainers/bcftools:1.24--h118bc1c_2@sha256:"
    "a3e0d3007ffe325c409b398f660840a3e7574d076219c6e82fc994ced87d47c3"
)
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
CONTIG = re.compile(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,|>)")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def read_header(path, expected_sample, expected_contigs):
    contigs, definitions, hc = [], [], []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##contig="):
                match = CONTIG.match(line)
                if match is None:
                    raise ValueError("Malformed contig declaration")
                contigs.append((match[1], int(match[2])))
            elif line.startswith(("##INFO=", "##FORMAT=", "##FILTER=", "##ALT=")):
                definitions.append(line.rstrip("\n"))
            elif line.startswith("##GATKCommandLine=<ID=HaplotypeCaller,"):
                hc.append(line.rstrip("\n"))
            elif line.startswith("#CHROM\t"):
                if line.rstrip("\n").split("\t")[9:] != [expected_sample]:
                    raise ValueError("gVCF must contain exactly its expected sample")
                break
            elif not line.startswith("#"):
                raise ValueError("Record before sample header")
        else:
            raise ValueError("Missing sample header")
    if contigs != expected_contigs:
        raise ValueError("Ordered gVCF contigs differ from reference")
    if len(hc) != 1 or ',Version="4.6.2.0"' not in hc[0]:
        raise ValueError("Expected exactly one pinned HaplotypeCaller provenance header")
    match = re.search(r'CommandLine="(.*)",Version=', hc[0])
    if match is None:
        raise ValueError("Missing HaplotypeCaller command")
    tokens = shlex.split(match[1])
    for option, expected in (
        ("--sample-ploidy", "2"),
        ("--native-pair-hmm-threads", "4"),
        ("--emit-ref-confidence", "GVCF"),
        ("--create-output-variant-index", "true"),
    ):
        if tokens.count(option) != 1 or tokens[tokens.index(option) + 1] != expected:
            raise ValueError(f"Unexpected HaplotypeCaller condition: {option}")
    # These file identities are independently tied to the recorded task/run.
    # Keep all other expanded defaults, not just the four explicit options.
    for option in ("--input", "--output", "--reference"):
        if tokens.count(option) != 1:
            raise ValueError(f"Expected one HaplotypeCaller {option}")
        tokens[tokens.index(option) + 1] = "VERIFIED_FILE_IDENTITY"
    return {
        "sample": expected_sample,
        "ordered_contigs_sha256": digest_json(contigs),
        "shared_header_definitions_sha256": digest_json(sorted(definitions)),
        "haplotypecaller_parameters_sha256": digest_json(tokens),
    }


def verify_bgzf(path):
    if path.stat().st_size < len(BGZF_EOF):
        raise ValueError("Empty or truncated BGZF")
    with path.open("rb") as handle:
        handle.seek(-len(BGZF_EOF), os.SEEK_END)
        if handle.read() != BGZF_EOF:
            raise ValueError("Missing BGZF EOF marker")
    # Read every gzip member to check CRC/size, including data beyond the header.
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            pass


def stream_records(lines, contigs, ploidy=2):
    """Constant-memory validation and digest of bcftools' canonical records."""
    rank = {name: (index, length) for index, (name, length) in enumerate(contigs)}
    digest = hashlib.sha256()
    count = 0
    previous = (-1, 0)
    for raw in lines:
        fields = raw.decode().rstrip("\n").split("\t")
        if len(fields) != 10 or fields[0] not in rank:
            raise ValueError("Malformed single-sample gVCF record")
        index, length = rank[fields[0]]
        position = int(fields[1])
        current = (index, position)
        if not 1 <= position <= length or current < previous:
            raise ValueError("gVCF record order or coordinate violation")
        formats = fields[8].split(":")
        if "GT" not in formats:
            raise ValueError("gVCF record has no GT")
        gt = fields[9].split(":")[formats.index("GT")]
        alleles = re.split(r"[/|]", gt)
        if len(alleles) != ploidy:
            raise ValueError("gVCF genotype ploidy differs from frozen lineage")
        allele_count = 1 + len(fields[4].split(","))
        if any(a != "." and (not a.isdigit() or int(a) >= allele_count) for a in alleles):
            raise ValueError("Invalid genotype allele index")
        for info in fields[7].split(";"):
            if info.startswith("END=") and not position <= int(info[4:]) <= length:
                raise ValueError("Invalid gVCF END")
        previous = current
        digest.update(raw)
        count += 1
    if count == 0:
        raise ValueError("Empty production gVCF")
    return {"record_count": count, "record_sha256": digest.hexdigest()}


def bcftools_records(path, contigs, indexed, log_dir):
    name = f"issue45-index-check-{os.getpid()}"
    command = [
        "docker",
        "run",
        "--rm",
        "--name",
        name,
        "--network",
        "none",
        "--cpus",
        "2",
        "--memory",
        "2g",
        "--memory-swap",
        "2g",
        "-v",
        f"{path.parent}:/input:ro",
        BCFTOOLS,
        "bcftools",
        "view",
        "-H",
    ]
    if indexed:
        command += ["--regions-overlap", "pos", "-r", ",".join(name for name, _ in contigs)]
    command += ["/input/" + path.name]
    with (log_dir / f"{path.name}.{'indexed' if indexed else 'sequential'}.log").open("x") as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
        try:
            result = stream_records(process.stdout, contigs)
            if process.wait() != 0:
                raise ValueError("bcftools rejected gVCF/index; consult retained log")
        finally:
            process.stdout.close()
            if process.poll() is None:
                subprocess.run(
                    ["docker", "stop", "--time", "10", name],
                    check=False,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                process.wait()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--production-checkout", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, production, output = (
        p.resolve() for p in (args.run_dir, args.production_checkout, args.output_dir)
    )
    manifest_path = root / "production_lineage_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["sample_count"] != 51 or len(manifest["samples"]) != 51:
        raise ValueError("Expected the frozen 51-sample cohort")
    samples = manifest["samples"]
    ids = [sample["run_accession"] for sample in samples]
    biological_ids = {sample["run_accession"]: sample["biosample"] for sample in samples}
    if len(set(ids)) != 51 or len({s["biosample"] for s in samples}) != 51:
        raise ValueError("Duplicate run or biological sample")
    if any(s["bioproject"] != "PRJNA1138464" for s in samples):
        raise ValueError("Unexpected public cohort")
    if (
        subprocess.check_output(
            ["git", "-C", str(production), "rev-parse", "HEAD"], text=True
        ).strip()
        != manifest["production_sha"]
    ):
        raise ValueError("Production SHA changed")
    if subprocess.check_output(
        ["git", "-C", str(production), "status", "--porcelain"], text=True
    ).strip():
        raise ValueError("Production checkout changed")
    operations = json.loads((root / "resource-preparation.private.json").read_text())
    if operations.get("nextflow_exit_code") != 0 or not operations.get("restored"):
        raise ValueError("Generation or service restoration did not finish successfully")
    if not re.search(
        r"N E X T F L O W\s+~\s+version 26\.04\.6\b", (root / "nextflow.stdout.log").read_text()
    ):
        raise ValueError("Actual Nextflow version does not match frozen lineage")
    for name, key in (
        ("generate.nf", "staged_wrapper_sha256"),
        ("nextflow.config", "nextflow_config_sha256"),
        ("samples.csv", "samplesheet_sha256"),
    ):
        if sha256(root / name) != manifest[key]:
            raise ValueError("Frozen launch file changed")
    reference = {}
    for role, record in manifest["reference"].items():
        path = root / "reference" / record["name"]
        if sha256(path) != record["sha256"]:
            raise ValueError("Reference changed after freeze")
        reference[role] = str(path)
    contigs = [
        (fields[0], int(fields[1]))
        for line in Path(reference["fai"]).read_text().splitlines()
        if (fields := line.split("\t"))
    ]
    dictionary = []
    for line in Path(reference["dict"]).read_text().splitlines():
        if line.startswith("@SQ\t"):
            row = dict(field.split(":", 1) for field in line.split("\t")[1:])
            dictionary.append((row["SN"], int(row["LN"])))
    if dictionary != contigs:
        raise ValueError("Reference FAI/dictionary order mismatch")
    trace = list(csv.DictReader((root / "trace.tsv").open(), delimiter="\t"))
    tasks = [row for row in trace if "GATK_HAPLOTYPECALLER (" in row["name"]]
    successful = [row for row in tasks if row["status"] == "COMPLETED" and row["exit"] == "0"]
    by_sample = {}
    for task in successful:
        sid = task["name"].rsplit(" (", 1)[1].removesuffix(")")
        if (
            sid in by_sample
            or sid not in ids
            or task["container"] != manifest["gatk_container"]
            or task["cpus"] != "4"
        ):
            raise ValueError("Unexpected successful HaplotypeCaller task identity")
        by_sample[sid] = task
    if set(by_sample) != set(ids):
        raise ValueError("Not all 51 samples have a successful recorded HaplotypeCaller task")
    if output.exists() or output == root or output == production:
        raise ValueError("Validation output must be a new directory")
    output.mkdir(parents=True)
    records = []
    try:
        for sid in ids:
            task = by_sample[sid]
            work = Path(task["workdir"]).resolve()
            if not work.is_relative_to(root / "work"):
                raise ValueError("HaplotypeCaller work is outside this generation")
            gvcf = root / "results/variants/gvcf" / f"{sid}.g.vcf.gz"
            index = Path(str(gvcf) + ".tbi")
            for path in (gvcf, index):
                if (
                    not path.is_file()
                    or path.stat().st_size == 0
                    or not path.samefile(work / path.name)
                ):
                    raise ValueError("Published gVCF/index is not this recorded task's output")
            command = (work / ".command.sh").read_text()
            if "HaplotypeCaller" not in command or f"--output {sid}.g.vcf.gz" not in command:
                raise ValueError("HaplotypeCaller task command mismatch")
            if manifest["gatk_container"] not in (work / ".command.run").read_text():
                raise ValueError("Actual task launcher used a different container")
            before = (sha256(gvcf), sha256(index))
            header = read_header(gvcf, sid, contigs)
            if records:
                for key in (
                    "shared_header_definitions_sha256",
                    "haplotypecaller_parameters_sha256",
                ):
                    if header[key] != records[0][key]:
                        raise ValueError(
                            "gVCF headers or HaplotypeCaller semantics are heterogeneous"
                        )
            verify_bgzf(gvcf)
            sequential = bcftools_records(gvcf, contigs, False, output)
            indexed = bcftools_records(gvcf, contigs, True, output)
            if sequential != indexed:
                raise ValueError(
                    "Index traversal does not reproduce the complete gVCF record stream"
                )
            if before != (sha256(gvcf), sha256(index)):
                raise ValueError("gVCF/index changed during validation")
            records.append(
                {
                    **header,
                    **sequential,
                    "biosample": biological_ids[sid],
                    "bioproject": "PRJNA1138464",
                    "gvcf_name": gvcf.name,
                    "gvcf_sha256": before[0],
                    "gvcf_index_sha256": before[1],
                    "gvcf_bytes": gvcf.stat().st_size,
                    "index_bytes": index.stat().st_size,
                    "index_matches_full_record_stream": True,
                    "bgzf_crc_and_eof_valid": True,
                    "task_attempt": int(task["attempt"]),
                    "task_peak_rss": task["peak_rss"],
                    "task_wall": task["realtime"],
                    "task_memory": task["memory"],
                }
            )
            (output / "progress.json").write_text(
                json.dumps({"validated": len(records), "sample": sid}) + "\n"
            )
        evidence = {
            "schema_version": 1,
            "cohort_id": manifest["cohort_id"],
            "sample_count": 51,
            "unique_biosamples": 51,
            "production_sha": manifest["production_sha"],
            "source_manifest_sha256": sha256(manifest_path),
            "reference_accession": manifest["reference_accession"],
            "reference": manifest["reference"],
            "nextflow_version": manifest["nextflow_version"],
            "gatk_version": manifest["gatk_version"],
            "gatk_container": manifest["gatk_container"],
            "bcftools_container": BCFTOOLS,
            "sample_ploidy": manifest["sample_ploidy"],
            "samples": records,
            "benchmark_ready": True,
            "lineage_verified": True,
            "validated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        private = {
            "evidence": evidence,
            "reference_paths": reference,
            "gvcf_directory": str(root / "results/variants/gvcf"),
            "source_run": str(root),
        }
        with tempfile.NamedTemporaryFile(mode="w", dir=output, delete=False) as handle:
            json.dump(private, handle, indent=2)
            handle.write("\n")
            temporary = Path(handle.name)
        temporary.replace(output / "validated_cohort.private.json")
        print(
            "Validated 51 unique public samples: lineage and complete indexed record streams agree"
        )
    except BaseException as error:
        (output / "failure.private.json").write_text(
            json.dumps({"error": str(error), "validated_samples": len(records)}) + "\n"
        )
        raise


if __name__ == "__main__":
    main()
