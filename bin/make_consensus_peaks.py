#!/usr/bin/env python3
"""Build union/consensus peak sets across comparable merged groups."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_table(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_peak(group_peaks: Path, macs_name: str, mode: str) -> Path | None:
    suffix = "broadPeak" if mode == "broad" else "narrowPeak"
    candidates = list(group_peaks.glob(f"**/{macs_name}_peaks.{suffix}"))
    return candidates[0] if candidates else None


def read_bed3(path: Path):
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                yield parts[0], int(parts[1]), int(parts[2])


def merge_intervals(records):
    by_chrom = defaultdict(list)
    for chrom, start, end in records:
        by_chrom[chrom].append((start, end))
    merged = []
    for chrom in sorted(by_chrom):
        intervals = sorted(by_chrom[chrom])
        if not intervals:
            continue
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chrom, cur_start, cur_end))
    return merged


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group-peaks", required=True)
    parser.add_argument("--consensus-jobs", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    group_peaks = Path(args.group_peaks)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    for job in read_table(Path(args.consensus_jobs)):
        macs_names = [x for x in job["macs_names"].split(";") if x]
        mode = job["peak_calling_mode"]
        records = []
        sources = []
        for macs_name in macs_names:
            peak = find_peak(group_peaks, macs_name, mode)
            if peak:
                sources.append(str(peak))
                records.extend(read_bed3(peak))
        merged = merge_intervals(records)
        out = outdir / f"{job['comparable_set_id']}.consensus.bed"
        with out.open("w") as handle:
            for i, (chrom, start, end) in enumerate(merged, start=1):
                handle.write(f"{chrom}\t{start}\t{end}\t{job['comparable_set_id']}_peak_{i}\t0\t.\n")
        rows.append({
            "comparable_set_id": job["comparable_set_id"],
            "peak_count": len(merged),
            "source_peak_files": ";".join(sources),
            "consensus_bed": str(out),
        })

    with (outdir / "consensus_peak_counts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["comparable_set_id", "peak_count", "source_peak_files", "consensus_bed"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with (outdir / "consensus_peak_counts_mqc.tsv").open("w", newline="") as handle:
        handle.write("# plot_type: bargraph\n")
        handle.write("# section_name: Consensus peak counts\n")
        handle.write("# description: Union peak counts for comparable non-control merged groups.\n")
        writer = csv.DictWriter(handle, fieldnames=["comparable_set_id", "peak_count"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

