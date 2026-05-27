#!/usr/bin/env python3
"""Validate CUT&Tag association CSV and create explicit Nextflow job tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id", "species", "genome", "antibody", "condition", "replicate",
    "group_id", "is_control", "control_group_id", "merge_group_id",
    "peak_calling_mode", "notes",
]
BROAD_HINTS = ("h3k27me3", "h3k36me3", "h3k9me3", "h3k79me2", "h4k20me1", "broad")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.+-]+$")


def norm_bool(value: str) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n", ""}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def sanitize(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(text).strip())
    return text.strip("_") or "unnamed"


def resolved_mode(mode: str, antibody: str) -> str:
    mode = (mode or "auto").strip().lower()
    if mode not in {"auto", "narrow", "broad"}:
        raise ValueError(f"peak_calling_mode must be narrow, broad, or auto; got {mode}")
    if mode != "auto":
        return mode
    antibody_l = antibody.lower()
    return "broad" if any(hint in antibody_l for hint in BROAD_HINTS) else "narrow"


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"association_csv is missing required columns: {', '.join(missing)}")
        return list(reader)


def write_table(path: Path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-pairs", required=True)
    parser.add_argument("--association-csv", required=True)
    parser.add_argument("--target-genome", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--allow-extra-association-rows", action="store_true")
    parser.add_argument("--allow-control-free-peak-calling", action="store_true")
    parser.add_argument("--call-control-peaks", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    errors, warnings = [], []
    target_genome = args.target_genome.strip()
    if not target_genome:
        errors.append("--target-genome must not be blank")

    pairs = read_tsv(Path(args.fastq_pairs))
    pair_rows = [r for r in pairs if r.get("sample_id") and r.get("status") == "OK"]
    pair_ids = {r["sample_id"] for r in pair_rows}
    detected_pair_ids = {r["sample_id"] for r in pairs if r.get("sample_id")}

    assoc = read_csv(Path(args.association_csv))
    assoc_by_id = {}
    for row in assoc:
        sid = row["sample_id"].strip()
        if not sid:
            errors.append("association_csv contains a blank sample_id")
            continue
        if not SAMPLE_ID_RE.match(sid):
            errors.append(f"sample_id contains unsupported characters: {sid}")
        if sid in assoc_by_id:
            errors.append(f"Duplicate sample_id in association_csv: {sid}")
        assoc_by_id[sid] = row
        try:
            row["is_control_norm"] = norm_bool(row["is_control"])
            row["peak_calling_mode_norm"] = resolved_mode(row["peak_calling_mode"], row["antibody"])
        except ValueError as exc:
            errors.append(f"{sid}: {exc}")

    assoc_ids = set(assoc_by_id)
    missing = sorted(pair_ids - assoc_ids)
    target_assoc_ids = {
        sid for sid, row in assoc_by_id.items()
        if str(row.get("genome", "")).strip() == target_genome
    }
    extra = sorted(target_assoc_ids - detected_pair_ids)
    if missing:
        errors.append("FASTQ sample_id values missing from association_csv: " + ", ".join(missing))
    if extra and not args.allow_extra_association_rows:
        errors.append(
            f"association_csv contains --genome {target_genome} sample_id rows without FASTQ pairs: "
            + ", ".join(extra)
        )
    elif extra:
        warnings.append(
            f"Ignoring extra --genome {target_genome} association rows without FASTQ pairs: "
            + ", ".join(extra)
        )

    selected_rows = []
    excluded_by_genome = []
    for pair in pair_rows:
        sid = pair["sample_id"]
        if sid not in assoc_by_id:
            continue
        row = dict(assoc_by_id[sid])
        row.update(pair)
        if str(row.get("genome", "")).strip() != target_genome:
            excluded_by_genome.append({
                "sample_id": sid,
                "association_genome": row.get("genome", ""),
                "target_genome": target_genome,
                "reason": "association genome does not match --genome",
            })
            continue
        for col in ["species", "genome", "antibody", "condition", "replicate", "group_id", "merge_group_id"]:
            if not str(row.get(col, "")).strip():
                errors.append(f"{sid}: required association field {col} is blank")
        selected_rows.append(row)

    pair_errors = [r for r in pairs if r.get("status") == "ERROR"]
    for pair_error in pair_errors:
        sid = pair_error.get("sample_id", "")
        assoc_row = assoc_by_id.get(sid)
        if assoc_row and str(assoc_row.get("genome", "")).strip() != target_genome:
            excluded_by_genome.append({
                "sample_id": sid,
                "association_genome": assoc_row.get("genome", ""),
                "target_genome": target_genome,
                "reason": f"non-target genome FASTQ pair has scanner error and was skipped: {pair_error.get('message', '')}",
            })
            continue
        errors.append(f"FASTQ pairing error for {sid or pair_error.get('r1')}: {pair_error.get('message')}")

    if not selected_rows:
        errors.append(f"No valid FASTQ pairs matched --genome {target_genome}.")

    merge_groups = defaultdict(list)
    group_ids_to_control_samples = defaultdict(list)
    merge_ids_to_control_samples = defaultdict(list)
    for row in selected_rows:
        merge_groups[row["merge_group_id"]].append(row)
        if row.get("is_control_norm"):
            group_ids_to_control_samples[row["group_id"]].append(row)
            merge_ids_to_control_samples[row["merge_group_id"]].append(row)

    for merge_group_id, members in merge_groups.items():
        for col in ["species", "genome", "antibody", "condition"]:
            values = {m[col] for m in members}
            if len(values) > 1:
                errors.append(f"merge_group_id {merge_group_id} has incompatible {col}: {', '.join(sorted(values))}")
        modes = {m.get("peak_calling_mode_norm", "") for m in members}
        if len(modes) > 1:
            errors.append(f"merge_group_id {merge_group_id} mixes peak_calling_mode values: {', '.join(sorted(modes))}")
        control_flags = {bool(m.get("is_control_norm")) for m in members}
        if len(control_flags) > 1:
            errors.append(f"merge_group_id {merge_group_id} mixes control and treatment samples")

    def resolve_control(control_group_id: str):
        cid = (control_group_id or "").strip()
        if not cid:
            return "", []
        if cid in merge_ids_to_control_samples:
            return cid, merge_ids_to_control_samples[cid]
        if cid in group_ids_to_control_samples:
            samples = group_ids_to_control_samples[cid]
            merge_ids = sorted({s["merge_group_id"] for s in samples})
            if len(merge_ids) > 1:
                errors.append(f"control_group_id {cid} maps to multiple control merge_group_id values: {', '.join(merge_ids)}")
            return merge_ids[0], samples
        return "", []

    for merge_group_id, members in merge_groups.items():
        first = members[0]
        if first.get("is_control_norm"):
            for m in members:
                m["resolved_control_merge_group_id"] = ""
            continue
        control_ids = {m.get("control_group_id", "").strip() for m in members}
        if len(control_ids) > 1:
            errors.append(f"merge_group_id {merge_group_id} has inconsistent control_group_id values: {', '.join(sorted(control_ids))}")
        control_id = next(iter(control_ids)) if control_ids else ""
        resolved, control_samples = resolve_control(control_id)
        if not resolved and not args.allow_control_free_peak_calling:
            errors.append(f"merge_group_id {merge_group_id} has no valid control_group_id; set --allow_control_free_peak_calling to override")
        if not resolved:
            warnings.append(f"merge_group_id {merge_group_id} will be peak-called without a control")
        for m in members:
            m["resolved_control_merge_group_id"] = resolved
            m["control_sample_ids"] = ";".join(sorted(s["sample_id"] for s in control_samples))

    if errors:
        report = outdir / "association_validation_report.txt"
        report.write_text("VALIDATION FAILED\n\nErrors:\n- " + "\n- ".join(errors) + "\n", encoding="utf-8")
        for err in errors:
            print(err, file=sys.stderr)
        return 2

    pair_by_id = {r["sample_id"]: r for r in pair_rows}
    samplesheet = []
    group_members = []
    sample_peak_jobs = []
    for row in selected_rows:
        sid = row["sample_id"]
        pair = pair_by_id[sid]
        is_control = bool(row["is_control_norm"])
        mode = row["peak_calling_mode_norm"]
        macs_name = sanitize(f"{sid}_{row['antibody']}_{row['condition']}")
        samplesheet.append({
            "sample_id": sid,
            "r1": pair["r1"],
            "r2": pair["r2"],
            "input_dir": pair["input_dir"],
            "detected_pattern": pair["detected_pattern"],
            "species": row["species"],
            "genome": row["genome"],
            "antibody": row["antibody"],
            "condition": row["condition"],
            "replicate": row["replicate"],
            "group_id": row["group_id"],
            "is_control": str(is_control).lower(),
            "control_group_id": row.get("control_group_id", ""),
            "resolved_control_merge_group_id": row.get("resolved_control_merge_group_id", ""),
            "merge_group_id": row["merge_group_id"],
            "peak_calling_mode": mode,
            "notes": row.get("notes", ""),
        })
        group_members.append({
            "sample_id": sid,
            "merge_group_id": row["merge_group_id"],
            "is_control": str(is_control).lower(),
            "control_group_id": row.get("control_group_id", ""),
            "resolved_control_merge_group_id": row.get("resolved_control_merge_group_id", ""),
            "peak_calling_mode": mode,
            "genome": row["genome"],
            "species": row["species"],
            "antibody": row["antibody"],
            "condition": row["condition"],
            "replicate": row["replicate"],
        })
        sample_peak_jobs.append({
            "sample_id": sid,
            "is_control": str(is_control).lower(),
            "control_group_id": row.get("control_group_id", ""),
            "resolved_control_merge_group_id": row.get("resolved_control_merge_group_id", ""),
            "control_sample_ids": row.get("control_sample_ids", ""),
            "peak_calling_mode": mode,
            "macs_name": macs_name,
            "species": row["species"],
            "genome": row["genome"],
            "antibody": row["antibody"],
            "condition": row["condition"],
            "merge_group_id": row["merge_group_id"],
        })

    group_peak_jobs = []
    for merge_group_id, members in sorted(merge_groups.items()):
        first = members[0]
        is_control = bool(first["is_control_norm"])
        comparable = sanitize(f"{first['species']}_{first['genome']}_{first['antibody']}")
        macs_name = sanitize(f"{merge_group_id}_{first['antibody']}_{first['condition']}")
        group_peak_jobs.append({
            "merge_group_id": merge_group_id,
            "is_control": str(is_control).lower(),
            "control_group_id": first.get("control_group_id", ""),
            "resolved_control_merge_group_id": first.get("resolved_control_merge_group_id", ""),
            "peak_calling_mode": first["peak_calling_mode_norm"],
            "macs_name": macs_name,
            "species": first["species"],
            "genome": first["genome"],
            "antibody": first["antibody"],
            "condition": first["condition"],
            "sample_ids": ";".join(sorted(m["sample_id"] for m in members)),
            "comparable_set_id": comparable,
        })

    consensus_rows = []
    by_comparable = defaultdict(list)
    for row in group_peak_jobs:
        if row["is_control"] != "true":
            by_comparable[row["comparable_set_id"]].append(row)
    for comparable, members in sorted(by_comparable.items()):
        if not members:
            continue
        consensus_rows.append({
            "comparable_set_id": comparable,
            "merge_group_ids": ";".join(m["merge_group_id"] for m in members),
            "macs_names": ";".join(m["macs_name"] for m in members),
            "peak_calling_mode": members[0]["peak_calling_mode"],
            "species": members[0]["species"],
            "genome": members[0]["genome"],
            "antibody": members[0]["antibody"],
        })

    write_table(outdir / "samplesheet.validated.tsv", [
        "sample_id", "r1", "r2", "input_dir", "detected_pattern", "species", "genome",
        "antibody", "condition", "replicate", "group_id", "is_control", "control_group_id",
        "resolved_control_merge_group_id", "merge_group_id", "peak_calling_mode", "notes"
    ], samplesheet)
    write_table(outdir / "group_members.tsv", [
        "sample_id", "merge_group_id", "is_control", "control_group_id",
        "resolved_control_merge_group_id", "peak_calling_mode", "genome", "species",
        "antibody", "condition", "replicate"
    ], group_members)
    write_table(outdir / "sample_peak_jobs.tsv", [
        "sample_id", "is_control", "control_group_id", "resolved_control_merge_group_id",
        "control_sample_ids", "peak_calling_mode", "macs_name", "species", "genome",
        "antibody", "condition", "merge_group_id"
    ], sample_peak_jobs)
    write_table(outdir / "group_peak_jobs.tsv", [
        "merge_group_id", "is_control", "control_group_id", "resolved_control_merge_group_id",
        "peak_calling_mode", "macs_name", "species", "genome", "antibody", "condition",
        "sample_ids", "comparable_set_id"
    ], group_peak_jobs)
    write_table(outdir / "consensus_jobs.tsv", [
        "comparable_set_id", "merge_group_ids", "macs_names", "peak_calling_mode",
        "species", "genome", "antibody"
    ], consensus_rows)
    write_table(outdir / "samples_excluded_by_genome.tsv", [
        "sample_id", "association_genome", "target_genome", "reason"
    ], excluded_by_genome)

    metadata = {
        "target_genome": target_genome,
        "sample_count": len(samplesheet),
        "excluded_by_genome_count": len(excluded_by_genome),
        "merge_group_count": len(group_peak_jobs),
        "consensus_set_count": len(consensus_rows),
        "warnings": warnings,
        "allow_control_free_peak_calling": args.allow_control_free_peak_calling,
        "call_control_peaks": args.call_control_peaks,
    }
    (outdir / "metadata.validation.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    report_lines = [
        "VALIDATION PASSED",
        "",
        f"Target genome: {target_genome}",
        f"Samples selected for processing: {len(samplesheet)}",
        f"Samples excluded by genome: {len(excluded_by_genome)}",
        f"Merge groups: {len(group_peak_jobs)}",
    ]
    if warnings:
        report_lines += ["", "Warnings:"] + [f"- {w}" for w in warnings]
    (outdir / "association_validation_report.txt").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
