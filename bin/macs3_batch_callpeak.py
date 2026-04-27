#!/usr/bin/env python3
"""Run MACS3 callpeak jobs from explicit sample/group job tables."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path


GENOME_ALIASES = {
    "hg18": "hs", "hg19": "hs", "hg38": "hs", "human": "hs",
    "mm9": "mm", "mm10": "mm", "mm39": "mm", "mouse": "mm",
    "rn6": "2.9e9", "rn7": "2.9e9", "rat": "2.9e9",
}


def read_table(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def infer_genome_size(value: str, genome: str) -> str:
    if value and value != "auto":
        return value
    return GENOME_ALIASES.get((genome or "").lower(), "hs")


def peak_file(outdir: Path, macs_name: str, mode: str) -> Path:
    suffix = "broadPeak" if mode == "broad" else "narrowPeak"
    return outdir / f"{macs_name}_peaks.{suffix}"


def count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def prepare_one(row, args, bam_map, script_dir: Path):
    is_control = row.get("is_control", "").lower() == "true"
    if args.level == "sample":
        job_id = row["sample_id"]
        treatment = bam_map.get(f"{job_id}.final.bam")
        controls = [bam_map.get(f"{sid}.final.bam") for sid in row.get("control_sample_ids", "").split(";") if sid]
    else:
        job_id = row["merge_group_id"]
        treatment = bam_map.get(f"{job_id}.merged.bam")
        ctrl = row.get("resolved_control_merge_group_id", "")
        controls = [bam_map.get(f"{ctrl}.merged.bam")] if ctrl else []

    controls = [c for c in controls if c]
    macs_name = row["macs_name"]
    mode = row.get("peak_calling_mode", "narrow")
    job_out = Path(args.outdir) / macs_name
    job_out.mkdir(parents=True, exist_ok=True)

    if is_control and not args.call_control_peaks:
        (job_out / "SKIPPED_CONTROL.txt").write_text("Control peak calling skipped by default.\n", encoding="utf-8")
        return {
            "id": job_id, "macs_name": macs_name, "level": args.level, "mode": mode,
            "status": "skipped_control", "peak_count": 0, "peak_file": "", "script": "",
        }, None
    if not treatment:
        raise RuntimeError(f"No BAM found for {args.level} job {job_id}")

    cmd = [
        "macs3", "callpeak",
        "-t", str(treatment),
        "-f", "BAMPE",
        "-n", macs_name,
        "-g", infer_genome_size(args.genome_size, row.get("genome", "")),
        "-q", str(args.qvalue),
        "--outdir", str(job_out),
    ]
    if controls:
        cmd += ["-c"] + [str(c) for c in controls]
    if mode == "broad":
        cmd += ["--broad", "--broad-cutoff", str(args.broad_cutoff)]
    if args.extra:
        cmd += shlex.split(args.extra)

    peaks = peak_file(job_out, macs_name, mode)
    log = job_out / f"{macs_name}.macs3.log"
    script = script_dir / f"{macs_name}.sh"
    quoted = " ".join(shlex.quote(str(x)) for x in cmd)
    script.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"echo {shlex.quote(quoted)} > {shlex.quote(str(log))}\n"
        f"{quoted} >> {shlex.quote(str(log))} 2>&1\n",
        encoding="utf-8",
    )
    return {
        "id": job_id, "macs_name": macs_name, "level": args.level, "mode": mode,
        "status": "pending", "peak_count": 0, "peak_file": str(peaks), "script": str(script),
    }, script


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", choices=["sample", "group"], required=True)
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--qvalue", type=float, required=True)
    parser.add_argument("--broad-cutoff", type=float, required=True)
    parser.add_argument("--genome-size", default="auto")
    parser.add_argument("--extra", default="")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--call-control-peaks", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bam_dir = Path(args.bam_dir)
    bam_map = {p.name: p.resolve() for p in bam_dir.glob("*.bam")}
    jobs = read_table(Path(args.jobs))
    results = []
    scripts = []
    script_dir = outdir / "parallel_jobs"
    script_dir.mkdir(parents=True, exist_ok=True)

    for row in jobs:
        result, script = prepare_one(row, args, bam_map, script_dir)
        results.append(result)
        if script:
            scripts.append(script)

    if scripts:
        script_list = outdir / "macs3_parallel_jobs.txt"
        script_list.write_text("\n".join(str(s) for s in scripts) + "\n", encoding="utf-8")
        cmd = ["parallel", "--will-cite", "--halt", "soon,fail=1", "-j", str(max(1, args.threads)), "bash", "::::", str(script_list)]
        completed = subprocess.run(cmd, text=True)
        if completed.returncode:
            raise RuntimeError("One or more MACS3 jobs failed; inspect per-job logs.")

    for row in results:
        if row["status"] == "pending":
            row["status"] = "ok"
            row["peak_count"] = count_lines(Path(row["peak_file"]))

    fields = ["level", "id", "macs_name", "mode", "status", "peak_count", "peak_file"]
    with (outdir / "peak_counts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda r: (r["level"], r["id"])))

    with (outdir / "peak_counts_mqc.tsv").open("w", newline="") as handle:
        handle.write("# plot_type: bargraph\n")
        handle.write(f"# section_name: MACS3 {args.level} peak counts\n")
        handle.write("# description: Number of called peaks after MACS3 callpeak.\n")
        writer = csv.DictWriter(handle, fieldnames=["id", "peak_count"], delimiter="\t")
        writer.writeheader()
        for row in sorted(results, key=lambda r: r["id"]):
            writer.writerow({"id": row["id"], "peak_count": row["peak_count"]})
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
