#!/usr/bin/env python3
"""Run the validated E0-E4 suite with an authorized local trial pause and restore."""

import argparse
import datetime
import json
import shutil
import signal
import stat
import subprocess
import threading
import time
from pathlib import Path

import execute_experiments
from execute_experiments import execute_suite
from run_generation import inspect, memory, serving, snapshot
from stage_generation import sha256

# Helpers whose content can change a measurement, and which a resumed run must
# therefore find unchanged. Reporting-only tools are deliberately excluded: a
# post-hoc summarizer edit must never block resuming a multi-day campaign, and
# recording it would make the lineage claim more than it means.
EXECUTING_HELPERS = (
    "benchmark_tools.py",
    "execute_experiments.py",
    "measure_process.py",
    "run_benchmarks.py",
    "run_generation.py",
    "stage_generation.py",
    "validate_generated_cohort.py",
)
REPORTING_ONLY_HELPERS = ("summarize_evidence.py",)

LAUNCH_MEMORY_GATE_BYTES = 110 * 1024**3
LAUNCH_STORAGE_GATE_BYTES = 2_000_000_000_000
RUNTIME_STORAGE_FLOOR_BYTES = 1_200_000_000_000


def occupied(root):
    """Bytes this run directory already holds, counting each file once.

    Regular files only, identified by inode so a hard-linked gVCF is not
    double counted, and never following a symlink out of the run directory.
    """
    seen, total = set(), 0
    for path in root.rglob("*"):
        try:
            info = path.lstat()
        except OSError:
            continue
        if not stat.S_ISREG(info.st_mode) or info.st_ino in seen:
            continue
        seen.add(info.st_ino)
        total += info.st_size
    return total


class ResourceGuard:
    def __init__(self, root, expected_trial_id):
        self.root = root
        self.expected_trial_id = expected_trial_id
        self.stop = threading.Event()
        self.monitor_done = threading.Event()
        self.monitor = None
        self.stopped_trial = False
        self.record = {}

    def save(self):
        (self.root / "resource-preparation.private.json").write_text(
            json.dumps(self.record, indent=2) + "\n"
        )

    def stop_owned(self):
        containers = subprocess.check_output(
            ["docker", "ps", "-q", "--filter", f"label=adzuki.issue45_benchmark={self.root.name}"],
            text=True,
        ).split()
        if containers:
            subprocess.run(
                ["docker", "stop", "--time", "30", *containers],
                check=False,
                stdout=subprocess.DEVNULL,
            )
        remaining = subprocess.check_output(
            ["docker", "ps", "-q", "--filter", f"label=adzuki.issue45_benchmark={self.root.name}"],
            text=True,
        ).split()
        if remaining:
            subprocess.run(["docker", "kill", *remaining], check=False, stdout=subprocess.DEVNULL)
            remaining = subprocess.check_output(
                [
                    "docker",
                    "ps",
                    "-q",
                    "--filter",
                    f"label=adzuki.issue45_benchmark={self.root.name}",
                ],
                text=True,
            ).split()
        if remaining:
            raise RuntimeError(
                "Owned benchmark containers remain active; restoring large trial is unsafe"
            )

    def monitor_host(self):
        bad = 0
        try:
            path = self.root / "host.tsv"
            new_file = not path.exists()
            with path.open("a") as handle:
                if new_file:
                    handle.write("timestamp\tmem_available_bytes\tswap_bytes\tstorage_free_bytes\n")
                while not self.monitor_done.is_set():
                    available, swap = memory()
                    free = shutil.disk_usage(self.root).free
                    handle.write(
                        f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{available}\t{swap}\t{free}\n"
                    )
                    handle.flush()
                    bad = (
                        bad + 1
                        if available < 8 * 1024**3
                        or swap - self.record["benchmark_swap_start_bytes"] > 512 * 1024**2
                        else 0
                    )
                    if bad >= 3 or free < RUNTIME_STORAGE_FLOOR_BYTES:
                        raise RuntimeError(
                            "Host available memory, swap delta, or storage guard failed"
                        )
                    self.monitor_done.wait(5)
        except BaseException as error:
            self.record["monitor_failure"] = str(error)
            self.stop.set()
            self.stop_owned()

    def __enter__(self):
        if self.expected_trial_id is None:
            return self.enter_without_pausing_any_workload()
        before = inspect("freetoken-qwen38-flash-next-trial")
        if len(self.expected_trial_id) < 12 or not before["Id"].startswith(self.expected_trial_id):
            raise ValueError("Inventoried local trial identity changed")
        if not before["State"]["Running"] or before["State"]["Paused"]:
            raise ValueError("Local trial is not in its inventoried running state")
        if before["HostConfig"]["PortBindings"] != {
            "1919/tcp": [{"HostIp": "127.0.0.1", "HostPort": "1919"}]
        }:
            raise ValueError("Trial exposure changed")
        if (
            subprocess.check_output(
                ["ss", "-Htn", "state", "established", "(", "sport", "=", ":1919", ")"], text=True
            ).strip()
            or not serving()
        ):
            raise ValueError("Trial is busy or its restoration baseline cannot be verified")
        self.before = before
        available, swap = memory()
        self.record = {
            "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "before": snapshot(before),
            "before_models_endpoint_200": True,
            "before_mem_available_bytes": available,
            "before_swap_bytes": swap,
            "storage_before_bytes": shutil.disk_usage(self.root).free,
            "stop_command": ["docker", "stop", "--time", "30", before["Id"]],
            "restore_command": ["docker", "start", before["Id"]],
            "reason": "Owner-authorized reversible local inference trial pause for Issue45 E0-E4 RAM headroom",
            "restored": False,
        }
        self.save()
        try:
            self.stopped_trial = True
            subprocess.run(self.record["stop_command"], check=True)
            available, swap = memory()
            self.record.update(
                after_stop_mem_available_bytes=available, benchmark_swap_start_bytes=swap
            )
            self.save()
            self.check_launch_headroom(available)
            self.monitor = threading.Thread(target=self.monitor_host, daemon=True)
            self.monitor.start()
            return self
        except BaseException:
            self.restore()
            raise

    def check_launch_headroom(self, available):
        """Gate the whole suite's headroom, crediting space this run already holds.

        A resumed run has spent part of its budget on experiments that are
        already finished, so comparing raw free space against the suite-wide
        gate would refuse to continue exactly when continuing is cheapest.
        """
        free = shutil.disk_usage(self.root).free
        budget = free + occupied(self.root)
        self.record.update(storage_free_at_launch_bytes=free, storage_budget_bytes=budget)
        if available < LAUNCH_MEMORY_GATE_BYTES or budget < LAUNCH_STORAGE_GATE_BYTES:
            raise ValueError("Benchmark launch headroom gate failed")

    def enter_without_pausing_any_workload(self):
        """Start measuring without stopping anything, when the host already has headroom."""
        self.before = None
        available, swap = memory()
        self.record = {
            "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "paused_workloads": [],
            "reason": "Host already had the required headroom; no workload was stopped",
            "before_mem_available_bytes": available,
            "before_swap_bytes": swap,
            "storage_before_bytes": shutil.disk_usage(self.root).free,
            "after_stop_mem_available_bytes": available,
            "benchmark_swap_start_bytes": swap,
            "restored": True,
        }
        self.save()
        try:
            self.check_launch_headroom(available)
            self.monitor = threading.Thread(target=self.monitor_host, daemon=True)
            self.monitor.start()
            return self
        except BaseException:
            self.restore()
            raise

    def restore(self):
        self.monitor_done.set()
        if self.monitor is not None:
            self.monitor.join(timeout=60)
        try:
            self.stop_owned()
            available, swap = memory()
            self.record.update(
                benchmark_end_mem_available_bytes=available, benchmark_swap_end_bytes=swap
            )
            if self.stopped_trial:
                result = subprocess.run(self.record["restore_command"], check=False)
                self.record["restore_exit_code"] = result.returncode
                for _ in range(120):
                    if serving():
                        break
                    time.sleep(5)
                after = inspect(self.before["Id"])
                self.record["after"] = snapshot(after)
                self.record["restored"] = (
                    result.returncode == 0
                    and after["State"]["Running"]
                    and serving()
                    and after["Config"] == self.before["Config"]
                    and after["HostConfig"] == self.before["HostConfig"]
                )
        except BaseException as error:
            self.record["restoration_failure"] = str(error)
            raise
        finally:
            self.record["finished_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            self.save()
        if self.stopped_trial and not self.record["restored"]:
            raise RuntimeError("Local trial restoration failed; explicit follow-up required")

    def __exit__(self, exc_type, error, traceback):
        if error is not None:
            self.record["benchmark_failure"] = str(error)
        self.restore()
        if self.stop.is_set() and error is None:
            raise RuntimeError("Benchmark failed its host resource guard")
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validated-cohort", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--expected-trial-id")
    group.add_argument(
        "--no-trial-pause",
        action="store_true",
        help="Measure without stopping any workload, when the host already has RAM headroom",
    )
    args = parser.parse_args()
    validated = json.loads(args.validated_cohort.read_text())
    evidence = validated["evidence"]
    if (
        evidence.get("sample_count") != 51
        or not evidence.get("lineage_verified")
        or not evidence.get("benchmark_ready")
    ):
        raise ValueError(
            "51-sample independent lineage validation is required before pausing any workload"
        )
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    helper = Path(__file__).resolve().parent
    lineage = json.dumps(
        {
            "validated_cohort_sha256": sha256(args.validated_cohort),
            "helper_sha256": {name: sha256(helper / name) for name in EXECUTING_HELPERS},
            "production_sha": evidence["production_sha"],
            "maximum_concurrent_tasks": execute_experiments.MAXIMUM_CONCURRENT_TASKS,
            "genomicsdb_memory_gib": execute_experiments.GENOMICSDB_MEMORY_GIB,
            "genotype_memory_gib": execute_experiments.GENOTYPE_MEMORY_GIB,
            "genotype_native_reserve_gib": execute_experiments.GENOTYPE_NATIVE_RESERVE_GIB,
            "batch_size": 50,
            "window_size_bp": 20_000_000,
            "small_scaffold_max_bp": 1_000_000,
            "benchmark_ready": True,
            "lineage_verified": True,
        },
        indent=2,
    )
    lineage_path = root / "execution_lineage.json"
    if lineage_path.exists():
        # Resuming: the cohort, helper code and resource policy must be unchanged,
        # or the completed experiments in this directory are no longer comparable.
        if lineage_path.read_text() != lineage + "\n":
            raise ValueError(
                "Cohort, helper checksums or resource policy changed since this run started; "
                "use a new output directory"
            )
    else:
        lineage_path.write_text(lineage + "\n")
    guard = ResourceGuard(root, None if args.no_trial_pause else args.expected_trial_id)

    def interrupted(signum, _frame):
        guard.stop.set()
        guard.stop_owned()
        raise KeyboardInterrupt(f"Received signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    with guard:
        execute_suite(root, validated, guard.stop)


if __name__ == "__main__":
    main()
