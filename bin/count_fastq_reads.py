#!/usr/bin/env python3
"""Count paired FASTQ records."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def opener(path: Path):
    return gzip.open if path.name.endswith(".gz") else open


def count_records(path: Path) -> int:
    n = 0
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for n, _ in enumerate(handle, start=1):
            pass
    if n % 4:
        raise SystemExit(f"FASTQ line count is not divisible by 4: {path}")
    return n // 4


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    r1 = count_records(Path(args.r1))
    r2 = count_records(Path(args.r2))
    if r1 != r2:
        raise SystemExit(f"Read count mismatch for {args.sample_id}: R1={r1}, R2={r2}")
    col = "raw_read_pairs" if args.stage == "raw" else f"{args.stage}_read_pairs"
    with open(args.output, "w", encoding="utf-8") as out:
        out.write(f"sample_id\tstage\t{col}\tread1_records\tread2_records\n")
        out.write(f"{args.sample_id}\t{args.stage}\t{r1}\t{r1}\t{r2}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

