"""Pinned tool execution and constant-memory joint-callset checks for Issue45."""

import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import time
from collections import Counter
from pathlib import Path

from stage_generation import sha256

GATK = (
    "broadinstitute/gatk:4.6.2.0@sha256:"
    "71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb"
)


def raw_qd(fields):
    """GATK 4.6.2 QualByDepth.getDepth with no read likelihoods at genotyping."""
    formats = fields[8].split(":")
    depth = restricted = 0
    for sample in fields[9:]:
        values = dict(zip(formats, sample.split(":")))
        alleles = re.split(r"[/|]", values["GT"])
        if "." in alleles or all(value == "0" for value in alleles):
            continue
        if values.get("AD", ".") != ".":
            ad = [int(value) for value in values["AD"].split(",")]
            total = sum(ad)
            if total:
                depth += total
                if total - ad[0] > 1:
                    restricted += total
                continue
        if values.get("DP", ".") != ".":
            depth += int(values["DP"])
    denominator = restricted or depth
    return float(fields[5]) / denominator if denominator and fields[5] != "." else None


def compare_callsets(first, second, contigs):
    """Match complete variant keys in dictionary order and classify every difference.

    QD jitter is reported separately, never silently removed from a digest. The
    high-QD branch is established from identical QUAL/AD/GT, with a conservative
    margin for QUAL's printed rounding. Current QD<2 filter membership must agree.
    """
    rank = {name: i for i, (name, _) in enumerate(contigs)}

    def positions(handle):
        current, group = None, {}
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            position = rank[fields[0]], int(fields[1])
            if current is not None and position != current:
                yield current, group
                group = {}
            current = position
            key = fields[3], fields[4]
            if key in group or len(group) >= 10000:
                raise ValueError("Duplicate/excessive variants at one position")
            group[key] = fields
        if current is not None:
            yield current, group

    counts = Counter()
    info_changes, format_changes, column_changes = Counter(), Counter(), Counter()
    with gzip.open(first, "rt") as a, gzip.open(second, "rt") as b:
        left, right = positions(a), positions(b)
        x, y = next(left, None), next(right, None)
        while x is not None or y is not None:
            if y is None or (x is not None and x[0] < y[0]):
                counts["left_only_variants"] += len(x[1])
                x = next(left, None)
                continue
            if x is None or y[0] < x[0]:
                counts["right_only_variants"] += len(y[1])
                y = next(right, None)
                continue
            counts["left_only_variants"] += len(x[1].keys() - y[1].keys())
            counts["right_only_variants"] += len(y[1].keys() - x[1].keys())
            for key in x[1].keys() & y[1].keys():
                first_fields, second_fields = x[1][key], y[1][key]
                counts["shared_variants"] += 1
                if first_fields == second_fields:
                    counts["identical_records"] += 1
                    continue
                counts["different_records"] += 1
                left_info = dict(
                    (item.split("=", 1) if "=" in item else (item, None))
                    for item in first_fields[7].split(";")
                )
                right_info = dict(
                    (item.split("=", 1) if "=" in item else (item, None))
                    for item in second_fields[7].split(";")
                )
                differing_info = {
                    key
                    for key in left_info.keys() | right_info.keys()
                    if key not in left_info
                    or key not in right_info
                    or left_info[key] != right_info[key]
                }
                info_changes.update(differing_info)
                for index, name in ((2, "ID"), (5, "QUAL"), (6, "FILTER"), (8, "FORMAT_order")):
                    if first_fields[index] != second_fields[index]:
                        column_changes[name] += 1
                for first_sample, second_sample in zip(first_fields[9:], second_fields[9:]):
                    lv = dict(zip(first_fields[8].split(":"), first_sample.split(":")))
                    rv = dict(zip(second_fields[8].split(":"), second_sample.split(":")))
                    format_changes.update(
                        {
                            key
                            for key in lv.keys() | rv.keys()
                            if key not in lv or key not in rv or lv[key] != rv[key]
                        }
                    )
                only_qd = (
                    differing_info == {"QD"}
                    and first_fields[:7] == second_fields[:7]
                    and first_fields[8:] == second_fields[8:]
                )
                if only_qd:
                    counts["qd_only_difference"] += 1
                    try:
                        raw = raw_qd(first_fields)
                        left_qd, right_qd = float(left_info["QD"]), float(right_info["QD"])
                    except (KeyError, ValueError, TypeError, IndexError):
                        # Missing annotation or unevaluable AD is an unexplained
                        # difference, never evidence of harmless QD jitter.
                        continue
                    if (
                        raw is not None
                        and math.isfinite(raw)
                        and raw >= 35.01
                        and math.isfinite(left_qd)
                        and math.isfinite(right_qd)
                    ):
                        counts["qd_only_high_qd_jitter_branch"] += 1
                        if (left_qd < 2) == (right_qd < 2):
                            counts["qd_jitter_current_filter_unchanged"] += 1
            x, y = next(left, None), next(right, None)
    result = {
        key: counts[key]
        for key in (
            "shared_variants",
            "left_only_variants",
            "right_only_variants",
            "identical_records",
            "different_records",
            "qd_only_difference",
            "qd_only_high_qd_jitter_branch",
            "qd_jitter_current_filter_unchanged",
        )
    }
    result.update(
        info_site_differences=dict(info_changes),
        format_sample_differences=dict(format_changes),
        column_site_differences=dict(column_changes),
    )
    result["all_differences_explained_by_high_qd_jitter_with_unchanged_current_filter"] = (
        counts["left_only_variants"] == counts["right_only_variants"] == 0
        and counts["different_records"] == counts["qd_jitter_current_filter_unchanged"]
    )
    return result


def filesystem_metrics(path):
    files, size = 0, 0
    for parent, _, names in os.walk(path):
        for name in names:
            entry = Path(parent) / name
            if entry.is_file() and not entry.is_symlink():
                files += 1
                size += entry.stat().st_size
    return {"bytes": size, "file_count": files}


def batches_from_log(text, sample_count, batch_size):
    observed = [
        (int(a), int(b)) for a, b in re.findall(r"Importing batch (\d+) with (\d+) samples", text)
    ]
    expected = [
        (i + 1, min(batch_size, sample_count - i * batch_size))
        for i in range((sample_count + batch_size - 1) // batch_size)
    ]
    # Multi-interval import may initialize the same batches for each partition.
    if (
        not observed
        or any(pair not in expected for pair in observed)
        or set(observed) != set(expected)
    ):
        raise ValueError("GATK logs do not prove all expected sample batches")
    completions = [
        (int(a), int(b)) for a, b in re.findall(r"Done importing batch (\d+)/(\d+)", text)
    ]
    if set(completions) != {(i, len(expected)) for i, _ in expected}:
        raise ValueError("GATK did not log completion of all batches")
    return {
        "expected": expected,
        "observed": observed,
        "completed": completions,
        "reader_initialization_fell_back_to_serial": "Falling back to serial VCF reader initialization"
        in text,
    }


def joint_integrity(path, expected_samples, contigs):
    rank = {name: (i, length) for i, (name, length) in enumerate(contigs)}
    definitions, actual_contigs = [], []
    records, genotypes, accounting = (hashlib.sha256() for _ in range(3))
    counts = {name: 0 for name, _ in contigs}
    count = called = 0
    previous = (-1, 0)
    keys_at_position = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##contig="):
                match = re.match(r"##contig=<ID=([^,>]+),length=(\d+)(?:,|>)", line)
                if not match:
                    raise ValueError("Invalid callset contig header")
                actual_contigs.append((match[1], int(match[2])))
            elif line.startswith(("##INFO=", "##FORMAT=", "##FILTER=", "##ALT=")):
                definitions.append(line.strip())
            elif line.startswith("#CHROM\t"):
                if line.rstrip("\n").split("\t")[9:] != expected_samples:
                    raise ValueError("Callset sample count/order differs from expected input")
                if actual_contigs != contigs:
                    raise ValueError("Callset ordered dictionary differs from reference")
                break
            elif not line.startswith("#"):
                raise ValueError("Record appeared before callset header")
        else:
            raise ValueError("Missing callset header")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(expected_samples) or fields[0] not in rank:
                raise ValueError("Malformed joint record")
            index, length = rank[fields[0]]
            current = index, int(fields[1])
            if current < previous or not 1 <= current[1] <= length:
                raise ValueError("Joint record order/coordinate violation")
            if current != previous:
                keys_at_position.clear()
            key = tuple(fields[:2] + fields[3:5])
            if key in keys_at_position:
                raise ValueError("Duplicate variant key, possibly overlapping gather intervals")
            keys_at_position.add(key)
            previous = current
            alternatives = fields[4].split(",")
            if "<NON_REF>" in alternatives or fields[4] == ".":
                raise ValueError("Unexpected non-variant/reference-confidence output")
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            formats = fields[8].split(":")
            if "GT" not in formats or not {"AC", "AN"} <= info.keys():
                raise ValueError("Missing GT/AC/AN accounting")
            gt_index = formats.index("GT")
            ac = [0] * len(alternatives)
            an, gts = 0, []
            for sample in fields[9:]:
                gt = sample.split(":")[gt_index]
                alleles = re.split(r"[/|]", gt)
                if len(alleles) != 2:
                    raise ValueError("Unexpected joint genotype ploidy")
                for allele in alleles:
                    if allele == ".":
                        continue
                    value = int(allele)
                    if not 0 <= value <= len(ac):
                        raise ValueError("Joint genotype allele outside ALT")
                    an += 1
                    if value:
                        ac[value - 1] += 1
                gts.append(gt)
            if [int(x) for x in info["AC"].split(",")] != ac or int(info["AN"]) != an:
                raise ValueError("INFO allele accounting disagrees with sample genotypes")
            counts[fields[0]] += 1
            count += 1
            called += an
            records.update(line.encode())
            genotypes.update(("\t".join([*key, *gts]) + "\n").encode())
            accounting.update(("\t".join([*key, info["AC"], info["AN"]]) + "\n").encode())
    return {
        "sample_count": len(expected_samples),
        "sample_order": expected_samples,
        "sample_order_sha256": hashlib.sha256("\n".join(expected_samples).encode()).hexdigest(),
        "ordered_contigs": contigs,
        "variant_count": count,
        "per_contig_variant_counts": counts,
        "called_alleles": called,
        "record_sha256": records.hexdigest(),
        "variant_gt_sha256": genotypes.hexdigest(),
        "accounting_sha256": accounting.hexdigest(),
        "shared_header_sha256": hashlib.sha256("\n".join(sorted(definitions)).encode()).hexdigest(),
        "vcf_sha256": sha256(path),
        "vcf_index_sha256": sha256(Path(str(path) + ".tbi")),
        "sample_order_valid": True,
        "contig_order_valid": True,
        "allele_accounting_valid": True,
    }


class ToolRunner:
    def __init__(
        self, root, input_dir, reference_dir, stop_event, cpu_limit=None, memory_limit=None
    ):
        self.root = root.resolve()
        self.input_dir = input_dir.resolve()
        self.reference_dir = reference_dir.resolve()
        self.stop_event = stop_event
        self.cpu_limit = cpu_limit
        self.memory_limit = memory_limit
        self.helper = Path(__file__).resolve().parent

    def effective_memory_gib(self, memory_gib):
        """Honour a synthetic-test cap without altering any scientific CLI argument."""
        return min(memory_gib, self.memory_limit) if self.memory_limit is not None else memory_gib

    def java_heap_gib(self, memory_gib, native_reserve_gib):
        """Size a JVM heap strictly below the container limit it has to live inside.

        Issue #45: a 15 GiB heap inside a 16 GiB container OOM-killed
        GenotypeGVCFs, because GenomicsDB's TileDB buffers, the JVM's own
        metaspace, thread stacks and code cache all live outside the heap.
        Deriving the heap from the effective limit keeps that reserve present
        at every allocation, including a capped synthetic one.
        """
        return max(1, self.effective_memory_gib(memory_gib) - native_reserve_gib)

    def run(self, directory, process, arguments, cpus=8, memory_gib=16):
        if self.stop_event.is_set():
            raise RuntimeError("Benchmark has been stopped by resource/failure guard")
        directory.mkdir(parents=True, exist_ok=True)
        measurement = directory / f"{process}.metrics.json"
        log_path = directory / f"{process}.log"
        if measurement.exists() or log_path.exists():
            raise ValueError("Tool output already exists; retain it and use a new run")
        container = f"issue45-{os.getpid()}-{directory.name}-{process}".lower()
        container = re.sub(r"[^a-z0-9_.-]", "-", container)
        relative = directory.relative_to(self.root)
        effective_cpus = min(cpus, self.cpu_limit) if self.cpu_limit is not None else cpus
        effective_memory = self.effective_memory_gib(memory_gib)
        command = [
            "docker",
            "run",
            "--name",
            container,
            "--network",
            "none",
            "--label",
            f"adzuki.issue45_benchmark={self.root.name}",
            "--cpus",
            str(effective_cpus),
            "--memory",
            f"{effective_memory}g",
            "--memory-swap",
            f"{effective_memory}g",
            "--user",
            f"{os.getuid()}:{os.getgid()}",
            "-v",
            f"{self.root}:/bench",
            "-v",
            f"{self.input_dir}:/input:ro",
            "-v",
            f"{self.reference_dir}:/reference:ro",
            "-v",
            f"{self.helper}:/helpers:ro",
            "-w",
            str(Path("/bench") / relative),
            GATK,
            "python3",
            "/helpers/measure_process.py",
            measurement.name,
            *arguments,
        ]
        (directory / f"{process}.command.private.json").write_text(
            json.dumps(command, indent=2) + "\n"
        )
        start = time.monotonic()
        with log_path.open("x") as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                exit_code = child.wait(timeout=24 * 3600)
            except BaseException:
                subprocess.run(
                    ["docker", "stop", "--time", "30", container],
                    check=False,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                child.wait()
                raise
        state_result = subprocess.run(
            ["docker", "inspect", "--format", "{{json .State}}", container],
            capture_output=True,
            text=True,
            check=False,
        )
        state = json.loads(state_result.stdout) if state_result.returncode == 0 else None
        # Remove only this finished measurement container; all data/logs are bind-mounted and retained.
        if state is not None and not state["Running"]:
            subprocess.run(["docker", "rm", container], check=True, stdout=subprocess.DEVNULL)
        result = (
            json.loads(measurement.read_text())
            if measurement.exists()
            else {"peak_rss_bytes": None}
        )
        result.update(
            docker_exit_code=exit_code,
            container_state=state,
            elapsed_including_container_seconds=time.monotonic() - start,
            allocated_cpus=effective_cpus,
            requested_cpus=cpus,
            allocated_memory_bytes=memory_gib * 1024**3,
            attempt=1,
            retry_count=0,
            process=process,
        )
        (directory / f"{process}.execution.json").write_text(json.dumps(result, indent=2) + "\n")
        if exit_code != 0 or result.get("exit_code") != 0 or (state and state.get("OOMKilled")):
            raise RuntimeError(f"{process} failed; retained measurement and container state")
        return result
