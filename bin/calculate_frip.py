#!/usr/bin/env python3
"""Calculate FRiP against sample, group, and consensus peak sets."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def read_table(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(cmd: str):
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def peak_path(root: Path, macs_name: str, mode: str):
    suffix = "broadPeak" if mode == "broad" else "narrowPeak"
    matches = list(root.glob(f"**/{macs_name}_peaks.{suffix}"))
    return matches[0] if matches else None


def consensus_path(root: Path, comparable_set_id: str):
    matches = list(root.glob(f"**/{comparable_set_id}.consensus.bed"))
    return matches[0] if matches else None


def make_fragments(bam: Path, outbed: Path):
    cmd = (
        f"bedtools bamtobed -bedpe -i {bam} | "
        "awk 'BEGIN{OFS=\"\\t\"} $1==$4 && $2>=0 && $5>=0 "
        "{s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e,$7,0,$9}' | "
        f"sort -k1,1 -k2,2n > {outbed}"
    )
    run(cmd)


def count_lines(path: Path) -> int:
    if not path or not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def overlap_count(fragments: Path, peaks: Path) -> int:
    if not peaks or not peaks.exists() or count_lines(peaks) == 0:
        return 0
    proc = subprocess.run(
        f"bedtools intersect -u -a {fragments} -b {peaks} | wc -l",
        shell=True, executable="/bin/bash", check=True, text=True, stdout=subprocess.PIPE,
    )
    return int(proc.stdout.strip() or 0)


def add_frip(rows, entity_type, entity_id, bam, peaks, peak_set_type, peak_set_id, fragment_dir):
    frag = fragment_dir / f"{entity_type}.{entity_id}.fragments.bed"
    if not frag.exists():
        make_fragments(bam, frag)
    total = count_lines(frag)
    in_peaks = overlap_count(frag, peaks) if peaks else 0
    rows.append({
        "entity_type": entity_type,
        "entity_id": entity_id,
        "peak_set_type": peak_set_type,
        "peak_set_id": peak_set_id,
        "total_fragments": total,
        "fragments_in_peaks": in_peaks,
        "frip": f"{in_peaks / total:.6f}" if total else "0",
        "peak_file": str(peaks) if peaks else "",
    })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--sample-peak-jobs", required=True)
    parser.add_argument("--group-peak-jobs", required=True)
    parser.add_argument("--sample-bam-dir", required=True)
    parser.add_argument("--group-bam-dir", required=True)
    parser.add_argument("--sample-peaks", required=True)
    parser.add_argument("--group-peaks", required=True)
    parser.add_argument("--consensus-peaks", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    fragment_dir = outdir / "fragments"
    fragment_dir.mkdir(parents=True, exist_ok=True)
    sample_bams = {p.name.replace(".final.bam", ""): p for p in Path(args.sample_bam_dir).glob("*.final.bam")}
    group_bams = {p.name.replace(".merged.bam", ""): p for p in Path(args.group_bam_dir).glob("*.merged.bam")}
    sample_jobs = {r["sample_id"]: r for r in read_table(Path(args.sample_peak_jobs))}
    group_jobs = {r["merge_group_id"]: r for r in read_table(Path(args.group_peak_jobs))}
    sample_peaks = Path(args.sample_peaks)
    group_peaks = Path(args.group_peaks)
    consensus_peaks = Path(args.consensus_peaks)
    rows = []

    for sample in read_table(Path(args.samplesheet)):
        sid = sample["sample_id"]
        bam = sample_bams.get(sid)
        if not bam:
            continue
        sjob = sample_jobs.get(sid)
        if sjob:
            peaks = peak_path(sample_peaks, sjob["macs_name"], sjob["peak_calling_mode"])
            add_frip(rows, "sample", sid, bam, peaks, "own_sample", sid, fragment_dir)
        gjob = group_jobs.get(sample["merge_group_id"])
        if gjob:
            peaks = peak_path(group_peaks, gjob["macs_name"], gjob["peak_calling_mode"])
            add_frip(rows, "sample", sid, bam, peaks, "matched_group", sample["merge_group_id"], fragment_dir)
            cpeaks = consensus_path(consensus_peaks, gjob["comparable_set_id"])
            add_frip(rows, "sample", sid, bam, cpeaks, "consensus", gjob["comparable_set_id"], fragment_dir)

    for gid, gjob in group_jobs.items():
        if gjob["is_control"] == "true":
            continue
        bam = group_bams.get(gid)
        if not bam:
            continue
        peaks = peak_path(group_peaks, gjob["macs_name"], gjob["peak_calling_mode"])
        add_frip(rows, "merged_group", gid, bam, peaks, "own_group", gid, fragment_dir)

    fields = ["entity_type", "entity_id", "peak_set_type", "peak_set_id", "total_fragments", "fragments_in_peaks", "frip", "peak_file"]
    with (outdir / "frip.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (outdir / "frip_mqc.tsv").open("w", newline="") as handle:
        handle.write("# plot_type: bargraph\n")
        handle.write("# section_name: FRiP scores\n")
        handle.write("# description: Fraction of fragments in called peak sets.\n")
        writer = csv.DictWriter(handle, fieldnames=["entity_id", "peak_set_type", "frip"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({"entity_id": row["entity_id"], "peak_set_type": row["peak_set_type"], "frip": row["frip"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

