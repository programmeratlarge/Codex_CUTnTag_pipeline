#!/usr/bin/env python3
"""Estimate per-TSS signal enrichment from BigWig tracks."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import pyBigWig


def read_tss(path: Path):
    records = []
    with path.open() as handle:
        for i, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[:3]
            name = parts[3] if len(parts) > 3 else f"TSS_{i}"
            center = (int(start) + int(end)) // 2
            records.append((chrom, center, name))
    return records


def mean_signal(bw, chrom, start, end):
    chroms = bw.chroms()
    if chrom not in chroms:
        return 0.0
    start = max(0, start)
    end = min(chroms[chrom], end)
    if end <= start:
        return 0.0
    values = bw.values(chrom, start, end, numpy=False)
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return sum(vals) / len(vals) if vals else 0.0


def write_plot(summary, score_map, outdir: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    labels = [r["sample_id"] for r in summary]
    medians = [float(r["median_tss_enrichment"]) for r in summary]
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.5), 4))
    ax.bar(labels, medians)
    ax.set_ylabel("Median TSS enrichment")
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "tss_enrichment_median.png", dpi=200)
    if score_map:
        fig, ax = plt.subplots(figsize=(max(6, len(score_map) * 0.5), 4))
        labels = list(score_map)
        ax.boxplot([score_map[label] for label in labels], labels=labels, showfliers=False)
        ax.set_ylabel("Per-TSS enrichment")
        ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / "tss_enrichment_distribution.png", dpi=200)


def median(values):
    values = sorted(values)
    if not values:
        return 0.0
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bigwigs", nargs="+", required=True)
    parser.add_argument("--tss-bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=3000)
    parser.add_argument("--center-window", type=int, default=100)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tss_records = read_tss(Path(args.tss_bed))
    rows = []
    summary = []
    score_map = {}
    flank_half = max(args.center_window, args.window // 2)

    for bw_path in map(Path, args.bigwigs):
        sample = bw_path.name.replace(".bw", "").replace(".bigWig", "").replace(".bigwig", "")
        scores = []
        bw = pyBigWig.open(str(bw_path))
        try:
            for chrom, center, name in tss_records:
                center_signal = mean_signal(bw, chrom, center - args.center_window // 2, center + args.center_window // 2)
                left = mean_signal(bw, chrom, center - args.window, center - flank_half)
                right = mean_signal(bw, chrom, center + flank_half, center + args.window)
                flank = (left + right) / 2
                score = center_signal / (flank + 1e-9)
                scores.append(score)
                rows.append({
                    "sample_id": sample,
                    "tss_id": name,
                    "center_signal": f"{center_signal:.6f}",
                    "flank_signal": f"{flank:.6f}",
                    "tss_enrichment": f"{score:.6f}",
                })
        finally:
            bw.close()
        score_map[sample] = scores
        summary.append({
            "sample_id": sample,
            "tss_count": len(scores),
            "mean_tss_enrichment": f"{sum(scores) / len(scores):.6f}" if scores else "0",
            "median_tss_enrichment": f"{median(scores):.6f}",
        })

    with (outdir / "tss_enrichment.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "tss_id", "center_signal", "flank_signal", "tss_enrichment"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (outdir / "tss_enrichment_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "tss_count", "mean_tss_enrichment", "median_tss_enrichment"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    with (outdir / "tss_enrichment_mqc.tsv").open("w", newline="") as handle:
        handle.write("# plot_type: bargraph\n")
        handle.write("# section_name: TSS enrichment\n")
        handle.write("# description: Median BigWig signal enrichment at TSS centers over flanking signal.\n")
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "median_tss_enrichment"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    write_plot(summary, score_map, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
