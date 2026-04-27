#!/usr/bin/env python3
"""Scan an input directory for paired-end FASTQ files and emit a pairing TSV."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


FASTQ_SUFFIX_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)

AUTO_PATTERNS = [
    ("illumina_R1_001", re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]\d+)?\.(?:fastq|fq)(?:\.gz)?$", re.I)),
    ("illumina_R1_lane", re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]L\d{3})?(?:[._-]\d+)?\.(?:fastq|fq)(?:\.gz)?$", re.I)),
    ("dot_R1", re.compile(r"^(?P<sample>.+?)\.R(?P<read>[12])\.(?:fastq|fq)(?:\.gz)?$", re.I)),
    ("underscore_1", re.compile(r"^(?P<sample>.+?)_(?P<read>[12])\.(?:fastq|fq)(?:\.gz)?$", re.I)),
    ("dot_1", re.compile(r"^(?P<sample>.+?)\.(?P<read>[12])\.(?:fastq|fq)(?:\.gz)?$", re.I)),
]


def compile_user_pattern(pattern: str):
    if not pattern:
        return None
    try:
        rx = re.compile(pattern)
    except re.error as exc:
        raise SystemExit(f"Invalid --paired_pattern regular expression: {exc}") from exc
    needed = {"sample", "read"}
    if not needed.issubset(rx.groupindex):
        raise SystemExit("--paired_pattern must contain named groups (?P<sample>...) and (?P<read>[12])")
    return ("user_pattern", rx)


def detect(path: Path, user_pattern=None):
    name = path.name
    patterns = [user_pattern] if user_pattern else AUTO_PATTERNS
    for label, rx in patterns:
        match = rx.match(name)
        if match:
            read = match.group("read").upper().replace("R", "")
            if read not in {"1", "2"}:
                continue
            sample = match.group("sample")
            return sample, read, label
    return None, None, None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--paired-pattern", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist or is not a directory: {input_dir}")

    user_pattern = compile_user_pattern(args.paired_pattern)
    fastqs = sorted(p for p in input_dir.rglob("*") if p.is_file() and FASTQ_SUFFIX_RE.search(p.name))
    groups = defaultdict(lambda: {"1": [], "2": [], "patterns": set(), "unmatched": []})
    unmatched = []

    for fastq in fastqs:
        sample, read, pattern = detect(fastq, user_pattern)
        if not sample:
            unmatched.append(fastq)
            continue
        groups[sample][read].append(fastq)
        groups[sample]["patterns"].add(pattern)

    rows = []
    for sample, info in sorted(groups.items()):
        r1s = info["1"]
        r2s = info["2"]
        patterns = ",".join(sorted(info["patterns"]))
        status = "OK"
        messages = []
        if len(r1s) != 1:
            status = "ERROR"
            messages.append(f"expected 1 R1, found {len(r1s)}")
        if len(r2s) != 1:
            status = "ERROR"
            messages.append(f"expected 1 R2, found {len(r2s)}")
        rows.append({
            "sample_id": sample,
            "r1": str(r1s[0].resolve()) if r1s else "",
            "r2": str(r2s[0].resolve()) if r2s else "",
            "input_dir": str(input_dir),
            "detected_pattern": patterns,
            "status": status,
            "message": "; ".join(messages),
        })

    for path in unmatched:
        rows.append({
            "sample_id": "",
            "r1": str(path.resolve()),
            "r2": "",
            "input_dir": str(input_dir),
            "detected_pattern": "unmatched",
            "status": "WARNING",
            "message": "FASTQ did not match any supported paired-end naming pattern",
        })

    if not rows:
        rows.append({
            "sample_id": "",
            "r1": "",
            "r2": "",
            "input_dir": str(input_dir),
            "detected_pattern": "",
            "status": "ERROR",
            "message": "No FASTQ files found",
        })

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", "r1", "r2", "input_dir", "detected_pattern", "status", "message"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    errors = [r for r in rows if r["status"] == "ERROR"]
    if errors:
        print(f"Detected {len(errors)} FASTQ pairing error(s); validation will fail.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

