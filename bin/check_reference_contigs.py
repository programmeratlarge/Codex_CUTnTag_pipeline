#!/usr/bin/env python3
"""Check whether reference-side files share BAM contig names."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def bam_contigs(path: Path) -> set[str]:
    proc = subprocess.run(
        ["samtools", "view", "-H", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    contigs = set()
    for line in proc.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                contigs.add(field[3:])
                break
    return contigs


def first_column_contigs(path: Path) -> set[str]:
    contigs = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            contigs.add(line.split("\t", 1)[0])
    return contigs


def gtf_contigs(path: Path) -> set[str]:
    contigs = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            contigs.add(line.split("\t", 1)[0])
    return contigs


def chr_style(contigs: set[str]) -> str:
    if not contigs:
        return "empty"
    prefixed = sum(1 for c in contigs if c.startswith("chr"))
    if prefixed == len(contigs):
        return "chr_prefixed"
    if prefixed == 0:
        return "unprefixed"
    return "mixed"


def summarize(label: str, path: str, contigs: set[str], bam: set[str], required: bool):
    if not path:
        return {
            "reference": label,
            "path": "",
            "status": "SKIPPED",
            "shared_contigs": 0,
            "reference_contigs": 0,
            "bam_style": chr_style(bam),
            "reference_style": "",
            "message": "not provided",
        }
    shared = bam & contigs
    status = "OK" if shared else ("ERROR" if required else "WARNING")
    message = "shared contig names detected"
    if not shared:
        message = "no shared contig names with BAM header"
        if chr_style(bam) != chr_style(contigs):
            message += f"; likely naming mismatch ({chr_style(bam)} BAM vs {chr_style(contigs)} reference)"
    return {
        "reference": label,
        "path": path,
        "status": status,
        "shared_contigs": len(shared),
        "reference_contigs": len(contigs),
        "bam_style": chr_style(bam),
        "reference_style": chr_style(contigs),
        "message": message,
    }


def maybe_load(path: str, loader):
    return loader(Path(path)) if path else set()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--blacklist", default="")
    parser.add_argument("--annotation-gtf", default="")
    parser.add_argument("--tss-bed", default="")
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    bam = bam_contigs(Path(args.bam))
    if not bam:
        raise SystemExit(f"No @SQ contigs found in BAM header: {args.bam}")

    rows = [
        summarize("chrom_sizes", args.chrom_sizes, first_column_contigs(Path(args.chrom_sizes)), bam, True),
        summarize("blacklist", args.blacklist, maybe_load(args.blacklist, first_column_contigs), bam, False),
        summarize("annotation_gtf", args.annotation_gtf, maybe_load(args.annotation_gtf, gtf_contigs), bam, False),
        summarize("tss_bed", args.tss_bed, maybe_load(args.tss_bed, first_column_contigs), bam, False),
    ]

    fields = [
        "reference", "path", "status", "shared_contigs", "reference_contigs",
        "bam_style", "reference_style", "message",
    ]
    with open(args.tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "REFERENCE CONTIG COMPATIBILITY",
        f"BAM: {args.bam}",
        f"BAM contigs: {len(bam)} ({chr_style(bam)})",
        "",
    ]
    for row in rows:
        lines.append(
            f"{row['status']}: {row['reference']} - {row['message']} "
            f"[shared={row['shared_contigs']}, reference_contigs={row['reference_contigs']}]"
        )
    Path(args.report).write_text("\n".join(lines) + "\n", encoding="utf-8")

    errors = [r for r in rows if r["status"] == "ERROR"]
    warnings = [r for r in rows if r["status"] == "WARNING"]
    for row in warnings:
        print(f"WARNING: {row['reference']}: {row['message']}", file=sys.stderr)
    if errors:
        for row in errors:
            print(f"ERROR: {row['reference']}: {row['message']}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

