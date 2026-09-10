#!/usr/bin/env python3
"""Run one staged generation, monitor headroom, restore an authorized trial.

The caller must already have permission to temporarily stop the explicitly named
local experimental inference container. No other service is stopped or changed.
Private operational records stay in the run directory; sanitize before publishing.
"""

import argparse
import datetime
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time
import urllib.request
from pathlib import Path

from stage_generation import sha256


def memory():
    rows = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    values = {key: int(value.split()[0]) * 1024 for key, value in rows.items()}
    return values["MemAvailable"], values["SwapTotal"] - values["SwapFree"]


def inspect(container):
    return json.loads(subprocess.check_output(["docker", "inspect", container]))[0]


def snapshot(state):
    # Keep equality fingerprints, not Config.Env (which may contain secrets).
    return {
        "id": state["Id"],
        "state": state["State"],
        "image": state["Image"],
        "config_sha256": hashlib.sha256(
            json.dumps(state["Config"], sort_keys=True).encode()
        ).hexdigest(),
        "host_config_sha256": hashlib.sha256(
            json.dumps(state["HostConfig"], sort_keys=True).encode()
        ).hexdigest(),
    }


def interrupted(signum, _frame):
    raise KeyboardInterrupt(f"Received signal {signum}")


def serving():
    try:
        with urllib.request.urlopen("http://127.0.0.1:1919/v1/models", timeout=5) as response:
            return response.status == 200
    except OSError:
        return False


def main():
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--expected-trial-id", required=True)
    args = parser.parse_args()
    root = args.run_dir.resolve()
    manifest = json.loads((root / "production_lineage_manifest.json").read_text())
    for name, key in (
        ("generate.nf", "staged_wrapper_sha256"),
        ("samples.csv", "samplesheet_sha256"),
        ("nextflow.config", "nextflow_config_sha256"),
    ):
        if sha256(root / name) != manifest[key]:
            raise ValueError(f"Staged {name} changed after lineage freeze")
    record_path = root / "resource-preparation.private.json"
    if record_path.exists():
        raise ValueError("This launch already has an operational record")
    trial = "freetoken-qwen38-flash-next-trial"
    before = inspect(trial)
    if not before["Id"].startswith(args.expected_trial_id) or len(args.expected_trial_id) < 12:
        raise ValueError("Trial container identity changed since inventory")
    if not before["State"]["Running"] or before["State"]["Paused"]:
        raise ValueError("Trial is not in its inventoried running state")
    bindings = before["HostConfig"]["PortBindings"]
    if bindings != {"1919/tcp": [{"HostIp": "127.0.0.1", "HostPort": "1919"}]}:
        raise ValueError("Trial exposure differs from inventory")
    connections = subprocess.check_output(
        ["ss", "-Htn", "state", "established", "(", "sport", "=", ":1919", ")"], text=True
    ).strip()
    if connections:
        raise ValueError("Trial has an active connection; do not interrupt a request")
    was_serving = serving()
    if not was_serving:
        raise ValueError("Trial is not serving normally; restoration baseline is uncertain")
    available, swap = memory()
    record = {
        "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "before": snapshot(before),
        "before_models_endpoint_200": was_serving,
        "before_mem_available_bytes": available,
        "before_swap_bytes": swap,
        "stop_command": ["docker", "stop", "--time", "30", before["Id"]],
        "restore_command": ["docker", "start", before["Id"]],
        "reason": "Explicitly authorized temporary stop of local stateless LLM trial to provide gVCF generation RAM headroom",
        "restored": False,
    }
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    stopped = False
    process = None
    error = None
    try:
        stopped = True
        subprocess.run(record["stop_command"], check=True)
        available, swap_start = memory()
        record["after_stop_mem_available_bytes"] = available
        record["generation_swap_start_bytes"] = swap_start
        record_path.write_text(json.dumps(record, indent=2) + "\n")
        if available < 110 * 1024**3:
            raise ValueError("Less than 110 GiB RAM available after trial stop")
        if shutil.disk_usage(root).free < 2_950_000_000_000:
            raise ValueError("Storage launch gate failed")
        environment = os.environ.copy()
        environment["NXF_VER"] = "26.04.6"
        environment["PATH"] = str(Path.home() / ".local/bin") + ":" + environment["PATH"]
        with (
            (root / "nextflow.stdout.log").open("x") as log,
            (root / "host.tsv").open("x") as monitor,
        ):
            monitor.write("timestamp\tmem_available_bytes\tswap_bytes\tstorage_free_bytes\n")
            process = subprocess.Popen(
                [
                    "nextflow",
                    "run",
                    "generate.nf",
                    "--outdir",
                    str(root / "results"),
                    "-w",
                    str(root / "work"),
                ],
                cwd=root,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                start_new_session=True,
            )
            bad_memory_polls = 0
            while process.poll() is None:
                available, swap_current = memory()
                free = shutil.disk_usage(root).free
                monitor.write(
                    f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{available}\t{swap_current}\t{free}\n"
                )
                monitor.flush()
                bad_memory_polls = (
                    bad_memory_polls + 1
                    if available < 8 * 1024**3 or swap_current - swap_start > 512 * 1024**2
                    else 0
                )
                if bad_memory_polls >= 3 or free < 1_200_000_000_000:
                    raise RuntimeError("Host memory/swap/storage safety gate failed")
                time.sleep(5)
            record["nextflow_exit_code"] = process.returncode
            if process.returncode:
                raise RuntimeError(f"Generation exited {process.returncode}")
    except BaseException as exception:
        error = exception
        record["failure"] = str(exception)
    finally:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=45)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        # Stop only containers created by this exact labeled generation, so
        # restoring the large inference model cannot compete with orphan tasks.
        owned = subprocess.check_output(
            ["docker", "ps", "-q", "--filter", f"label=adzuki.issue45_run={root.name}"], text=True
        ).split()
        if owned:
            stop = subprocess.run(["docker", "stop", "--time", "30", *owned], check=False)
            if stop.returncode:
                subprocess.run(["docker", "kill", *owned], check=False)
            record["orphan_task_containers_stopped"] = owned
            remaining = subprocess.check_output(
                ["docker", "ps", "-q", "--filter", f"label=adzuki.issue45_run={root.name}"],
                text=True,
            ).split()
            if remaining:
                record["restoration_failure"] = (
                    "Own generation containers could not be stopped; restoring the large model would be unsafe"
                )
                record_path.write_text(json.dumps(record, indent=2) + "\n")
                raise RuntimeError(record["restoration_failure"])
        record["generation_end_mem_available_bytes"], record["generation_swap_end_bytes"] = memory()
        if stopped:
            result = subprocess.run(record["restore_command"], check=False)
            record["restore_exit_code"] = result.returncode
            for _ in range(120):
                if serving():
                    break
                time.sleep(5)
            after = inspect(before["Id"])
            record["after"] = snapshot(after)
            record["restored"] = (
                result.returncode == 0
                and after["State"]["Running"]
                and serving()
                and after["Config"] == before["Config"]
                and after["HostConfig"] == before["HostConfig"]
            )
        record["finished_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        record_path.write_text(json.dumps(record, indent=2) + "\n")
    if error is not None:
        raise error
    if not record["restored"]:
        raise RuntimeError("Trial service restoration needs attention")


if __name__ == "__main__":
    main()
