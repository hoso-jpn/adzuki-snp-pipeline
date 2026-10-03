#!/usr/bin/env python3
"""Validate paired gVCF/index identities and one-sample headers before import."""

import argparse
import gzip
import json
from pathlib import Path


def build_map(entries: list[dict], directory: Path) -> tuple[str, dict]:
    if not entries:
        raise ValueError("empty gVCF cohort")
    seen_ids, seen_files, rows = set(), set(), []
    for entry in sorted(entries, key=lambda item: item["sample_id"]):
        sample, gvcf, index = (entry[key] for key in ("sample_id", "gvcf", "index"))
        if (
            not sample
            or any(c in sample for c in "\t\r\n")
            or sample in seen_ids
            or gvcf in seen_files
            or index != gvcf + ".tbi"
            or any(Path(name).name != name for name in (gvcf, index))
        ):
            raise ValueError("duplicate sample/file or mismatched gVCF/index tuple")
        if not (directory / index).is_file() or (directory / index).stat().st_size == 0:
            raise ValueError("missing/empty gVCF index")
        header = None
        with gzip.open(directory / gvcf, "rt") as handle:
            for line in handle:
                if line.startswith("#CHROM\t"):
                    header = line.rstrip("\n").split("\t")[9:]
                    break
                if not line.startswith("#"):
                    break
        if header != [sample]:
            raise ValueError(f"gVCF header does not match declared sample {sample}")
        seen_ids.add(sample)
        seen_files.add(gvcf)
        rows.append(f"{sample}\t{gvcf}\t{index}\n")
    return "".join(rows), {
        "schema_version": 1,
        "sample_order": sorted(seen_ids),
        "sample_count": len(seen_ids),
        "pairing": "declared_tuple_and_header",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    try:
        text, summary = build_map(json.loads(args.entries.read_text()), Path.cwd())
        args.output.write_text(text)
        args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"sample map rejected: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
