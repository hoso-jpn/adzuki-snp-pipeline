#!/usr/bin/env python3
"""Measure one tool and its descendants inside its resource-limited container."""

import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path


def main():
    destination = Path(sys.argv[1])
    if destination.exists():
        raise SystemExit("Measurement destination already exists")
    started = time.monotonic()
    error = None
    try:
        result = subprocess.run(sys.argv[2:], check=False)
        exit_code = result.returncode
    except OSError as exception:
        exit_code = 127
        error = str(exception)
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    record = {
        "wall_seconds": time.monotonic() - started,
        "peak_rss_bytes": usage.ru_maxrss * 1024,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "exit_code": exit_code,
        "launch_error": error,
        "major_page_faults": usage.ru_majflt,
        "filesystem_input_blocks": usage.ru_inblock,
        "filesystem_output_blocks": usage.ru_oublock,
        "measurement_scope": "Linux RUSAGE_CHILDREN: tool and waited-for descendants",
    }
    temporary = destination.with_name(destination.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(destination)
    return exit_code if exit_code >= 0 else 128 - exit_code


if __name__ == "__main__":
    sys.exit(main())
