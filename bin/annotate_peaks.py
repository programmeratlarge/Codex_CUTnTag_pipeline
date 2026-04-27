#!/usr/bin/env python3
"""Lightweight peak annotation from GTF/GFF plus optional TSS BED."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


BIN_SIZE = 100_000
PEAK_EXTS = (".narrowPeak", ".broadPeak", ".bed")


def attrs(text: str):
    values = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, val = parts
        values[key] = val.strip().strip('"')
    return values


def bins_for(start: int, end: int):
    return range(max(0, start) // BIN_SIZE, max(start, end - 1) // BIN_SIZE + 1)


def add(index, chrom, start, end, payload):
    if end <= start:
        return
    for b in bins_for(start, end):
        index[chrom][b].append((start, end, payload))


def query(index, chrom, start, end):
    seen = set()
    for b in bins_for(start, end):
        for i, rec in enumerate(index.get(chrom, {}).get(b, [])):
            key = (b, i)
            if key in seen:
                continue
            seen.add(key)
            rstart, rend, payload = rec
            if rstart < end and rend > start:
                yield rec


def load_gtf(path: Path, downstream_window: int):
    exons = defaultdict(lambda: defaultdict(list))
    genes = defaultdict(lambda: defaultdict(list))
    downstream = defaultdict(lambda: defaultdict(list))
    tss_records = []
    gene_spans = {}

    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr_text = parts[:9]
            start0, end1 = int(start) - 1, int(end)
            info = attrs(attr_text)
            gene_id = info.get("gene_id") or info.get("ID") or info.get("Parent") or "unknown"
            gene_name = info.get("gene_name") or info.get("Name") or gene_id
            if feature == "gene":
                gene_spans[gene_id] = (chrom, start0, end1, strand, gene_name)
            if feature == "exon":
                add(exons, chrom, start0, end1, gene_name)

    for gene_id, (chrom, start, end, strand, gene_name) in gene_spans.items():
        add(genes, chrom, start, end, gene_name)
        if strand == "-":
            tss = end
            dstart, dend = max(0, start - downstream_window), start
        else:
            tss = start
            dstart, dend = end, end + downstream_window
        tss_records.append((chrom, tss, tss + 1, gene_name, strand))
        add(downstream, chrom, dstart, dend, gene_name)
    return exons, genes, downstream, tss_records


def load_tss_bed(path: str, derived_tss):
    if not path:
        return derived_tss
    bed = Path(path)
    if not bed.exists():
        return derived_tss
    records = []
    with bed.open() as handle:
        for i, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[:3]
            name = parts[3] if len(parts) > 3 else f"TSS_{i}"
            strand = parts[5] if len(parts) > 5 else "."
            center = (int(start) + int(end)) // 2
            records.append((chrom, center, center + 1, name, strand))
    return records or derived_tss


def promoter_index(tss_records, upstream: int, downstream: int):
    index = defaultdict(lambda: defaultdict(list))
    for chrom, start, end, name, strand in tss_records:
        tss = start
        if strand == "-":
            pstart, pend = max(0, tss - downstream), tss + upstream
        else:
            pstart, pend = max(0, tss - upstream), tss + downstream
        add(index, chrom, pstart, pend, name)
    return index


def classify(chrom, start, end, promoter, exons, genes, downstream):
    if any(query(promoter, chrom, start, end)):
        return "promoter_or_TSS_proximal"
    if any(query(exons, chrom, start, end)):
        return "exon"
    if any(query(genes, chrom, start, end)):
        return "intron"
    if any(query(downstream, chrom, start, end)):
        return "downstream"
    return "intergenic"


def peak_files(root: Path):
    return sorted(p for p in root.rglob("*") if p.is_file() and p.name.endswith(PEAK_EXTS))


def safe_name(path: Path):
    stem = path.name
    for ext in PEAK_EXTS:
        stem = stem.replace(ext, "")
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", stem)


def write_plot(summary_rows, outdir: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    by_peak = defaultdict(dict)
    for row in summary_rows:
        by_peak[row["peak_set"]][row["category"]] = int(row["count"])
    cats = sorted({r["category"] for r in summary_rows})
    labels = list(by_peak)
    bottoms = [0] * len(labels)
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.45), 4))
    for cat in cats:
        values = [by_peak[label].get(cat, 0) for label in labels]
        ax.bar(labels, values, bottom=bottoms, label=cat)
        bottoms = [a + b for a, b in zip(bottoms, values)]
    ax.set_ylabel("Peak count")
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "annotation_category_distribution.png", dpi=200)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--tss-bed", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--promoter-upstream", type=int, default=2000)
    parser.add_argument("--promoter-downstream", type=int, default=500)
    parser.add_argument("--downstream-window", type=int, default=3000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    exons, genes, downstream, derived_tss = load_gtf(Path(args.gtf), args.downstream_window)
    promoters = promoter_index(load_tss_bed(args.tss_bed, derived_tss), args.promoter_upstream, args.promoter_downstream)
    summary_rows = []

    for peak_path in peak_files(Path(args.peaks)):
        peak_set = safe_name(peak_path)
        counts = Counter()
        annotated = outdir / f"{peak_set}.annotated.tsv"
        with peak_path.open() as inp, annotated.open("w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["chrom", "start", "end", "peak_name", "annotation"])
            for i, line in enumerate(inp, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name = parts[3] if len(parts) > 3 else f"{peak_set}_{i}"
                category = classify(chrom, start, end, promoters, exons, genes, downstream)
                counts[category] += 1
                writer.writerow([chrom, start, end, name, category])
        total = sum(counts.values()) or 1
        for category, count in sorted(counts.items()):
            summary_rows.append({
                "peak_set": peak_set,
                "category": category,
                "count": count,
                "fraction": f"{count / total:.6f}",
            })

    with (outdir / "annotation_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["peak_set", "category", "count", "fraction"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    with (outdir / "annotation_category_mqc.tsv").open("w", newline="") as handle:
        handle.write("# plot_type: bargraph\n")
        handle.write("# section_name: Peak annotation categories\n")
        handle.write("# description: Genomic annotation distribution of called peaks.\n")
        writer = csv.DictWriter(handle, fieldnames=["peak_set", "category", "count"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    write_plot(summary_rows, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

