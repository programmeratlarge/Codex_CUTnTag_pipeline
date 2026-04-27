#!/usr/bin/env python3
"""Merge per-stage read counts and write MultiQC custom content."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


STAGES = [
    "raw_read_pairs",
    "trimmed_read_pairs",
    "aligned_fragments",
    "properly_paired_fragments",
    "mapq_proper_fragments",
    "mitochondrial_removed_fragments",
    "blacklist_filtered_fragments",
    "pre_dedup_fragments",
    "final_usable_fragments",
]


def read_rows(paths):
    for path in paths:
        with Path(path).open(newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def merge(files_by_type):
    data = {}
    for kind, paths in files_by_type.items():
        for row in read_rows(paths):
            sid = row["sample_id"]
            data.setdefault(sid, {"sample_id": sid})
            for key, value in row.items():
                if key == "sample_id" or value == "":
                    continue
                data[sid][key] = value
    return data


def as_int(row, key):
    try:
        return int(float(row.get(key, 0) or 0))
    except ValueError:
        return 0


def add_percentages(row):
    raw = as_int(row, "raw_read_pairs")
    final = as_int(row, "final_usable_fragments")
    pre_dedup = as_int(row, "pre_dedup_fragments")
    mito_before = as_int(row, "mapq_proper_fragments")
    mito_after = as_int(row, "mitochondrial_removed_fragments")
    bl_after = as_int(row, "blacklist_filtered_fragments")
    row["final_retention_fraction"] = f"{final / raw:.6f}" if raw else "0"
    row["mitochondrial_fraction_removed"] = f"{(mito_before - mito_after) / mito_before:.6f}" if mito_before else "0"
    row["blacklist_fraction_removed"] = f"{(mito_after - bl_after) / mito_after:.6f}" if mito_after else "0"
    row["duplicate_fraction_removed"] = f"{(pre_dedup - final) / pre_dedup:.6f}" if pre_dedup else row.get("duplicate_rate", "0")
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--multiqc", required=True)
    parser.add_argument("--raw", nargs="+", required=True)
    parser.add_argument("--trimmed", nargs="+", required=True)
    parser.add_argument("--filter", nargs="+", required=True)
    parser.add_argument("--dedup", nargs="+", required=True)
    args = parser.parse_args()

    data = merge({"raw": args.raw, "trimmed": args.trimmed, "filter": args.filter, "dedup": args.dedup})
    rows = [add_percentages(data[sid]) for sid in sorted(data)]
    fields = ["sample_id"] + STAGES + [
        "duplicate_removed_fragments", "duplicate_rate", "final_retention_fraction",
        "mitochondrial_fraction_removed", "blacklist_fraction_removed", "duplicate_fraction_removed",
    ]
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    with open(args.multiqc, "w", newline="") as handle:
        handle.write("# plot_type: linegraph\n")
        handle.write("# section_name: CUT&Tag read retention\n")
        handle.write("# description: Paired-end fragments retained at major processing stages.\n")
        writer = csv.DictWriter(handle, fieldnames=["sample_id"] + STAGES, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

