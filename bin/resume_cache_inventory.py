#!/usr/bin/env python3
"""Inventory completed local tasks after a run stops; verify before -resume.

This does not delete, repair, or silently bless changed cache files. Snapshots
must be taken while Nextflow is stopped, before the cache can be modified.
Use a new work directory after corruption, or explicitly remove/rebuild the
affected task. Input symlinks and command/log files are not cached outputs.
"""

import argparse
import json
import os
from pathlib import Path

from manifest_utils import sha256_file


def snapshot(root: Path) -> dict:
    root = root.resolve(strict=True)
    files, tasks = [], 0
    for exit_file in sorted(root.glob("*/*/.exitcode")):
        if not exit_file.parent.resolve().is_relative_to(root):
            raise ValueError("task directory resolves outside work-dir")
        if exit_file.is_symlink() or exit_file.read_text().strip() != "0":
            continue
        task = exit_file.parent
        tasks += 1
        for directory, subdirs, names in os.walk(task, followlinks=False):
            subdirs[:] = sorted(
                name
                for name in subdirs
                if not (Path(directory) == task and name.startswith("."))
                and not (Path(directory) / name).is_symlink()
            )
            for name in sorted(names):
                path = Path(directory) / name
                if (
                    (Path(directory) == task and name.startswith("."))
                    or path.is_symlink()
                    or not path.is_file()
                ):
                    continue
                files.append(
                    {
                        "path": str(path.relative_to(root)),
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                )
    return {
        "schema_version": 1,
        "contract": "stopped_local_cache_v1",
        "successful_tasks": tasks,
        "files": files,
    }


def verify(root: Path, manifest: dict) -> list[dict]:
    root = root.resolve(strict=True)
    if manifest.get("schema_version") != 1 or manifest.get("contract") != "stopped_local_cache_v1":
        raise ValueError("unsupported cache inventory")
    issues, seen, task_paths = [], set(), set()
    for row in manifest["files"]:
        relative = Path(row["path"])
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or len(relative.parts) < 3
            or str(relative) in seen
        ):
            raise ValueError("unsafe or duplicate cache inventory path")
        seen.add(str(relative))
        task_paths.add(Path(*relative.parts[:2]))
        path = root / relative
        if not path.exists():
            # Nextflow may cache a declared directory while an internal fragment
            # is missing. Only absent top-level files have the existence-based
            # re-execution behavior exercised by the synthetic resume test.
            reason = "missing" if len(relative.parts) == 3 else "missing_directory_member"
        elif path.is_symlink() or not path.resolve().is_relative_to(root):
            reason = "symlink_or_outside_workdir"
        elif (
            not path.is_file()
            or path.stat().st_size != row["bytes"]
            or sha256_file(path) != row["sha256"]
        ):
            reason = "changed"
        else:
            continue
        issues.append({"path": str(relative), "reason": reason})
    # A new GenomicsDB fragment can change a cached directory just as a
    # modified/deleted one can. Verify membership without hashing every large
    # cached output a second time. Unrelated tasks are outside this inventory.
    for task_relative in sorted(task_paths):
        task = root / task_relative
        if not task.is_dir() or task.is_symlink() or not task.resolve().is_relative_to(root):
            continue  # The recorded members already report loss or unsafe paths.
        for directory, subdirs, names in os.walk(task, followlinks=False):
            subdirs[:] = sorted(
                name for name in subdirs if not (Path(directory) == task and name.startswith("."))
            )
            for name in sorted(
                names + [name for name in subdirs if (Path(directory) / name).is_symlink()]
            ):
                path = Path(directory) / name
                if Path(directory) == task and (name.startswith(".") or path.is_symlink()):
                    continue  # Top-level staging symlinks are inputs, not outputs.
                relative = str(path.relative_to(root))
                if relative not in seen:
                    issues.append({"path": relative, "reason": "unexpected_directory_member"})
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["snapshot", "verify"])
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument(
        "--permit-missing",
        action="store_true",
        help="Permit missing top-level files only; directory-member loss remains blocked.",
    )
    args = parser.parse_args()
    try:
        if args.action == "snapshot":
            if args.inventory.resolve().is_relative_to(args.work_dir.resolve()):
                raise ValueError("write the inventory outside work-dir")
            args.inventory.write_text(
                json.dumps(snapshot(args.work_dir), indent=2, sort_keys=True) + "\n"
            )
            return 0
        issues = verify(args.work_dir, json.loads(args.inventory.read_text()))
        print(
            json.dumps(
                {
                    "issues": issues,
                    "safe_to_reuse_unchanged_files": not issues,
                    "missing_tasks_require_reexecution": any(
                        x["reason"] == "missing" for x in issues
                    ),
                },
                indent=2,
            )
        )
        return (
            0
            if not issues or (args.permit_missing and all(x["reason"] == "missing" for x in issues))
            else 2
        )
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(2, f"cache inventory rejected: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
