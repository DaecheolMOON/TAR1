#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse FASTA headers like '>::chr19q:25658-27061', merge multiple ranges per arm,
and write a CSV with columns: arm, min_pos, max_pos, span (span = max_pos - min_pos).

Notes:
- Only headers are parsed; sequence lines are ignored.
- Multiple occurrences for the same arm are merged by taking global min(start) and max(end).
- Span follows the user's instruction as (max - min), not inclusive length.
- Optional filter: records where start or end exceeds a given --max-pos value can be excluded.
"""

import argparse
import csv
import os
import re
from collections import defaultdict

HDR_RE = re.compile(
    r"^>\s*[:]*\s*(chr[0-9XY]+[pq]?)\s*:\s*(\d+)\s*-\s*(\d+)\s*$"
)


def parse_headers(path, max_pos=None):
    """Yield (arm, start, end) tuples from FASTA headers. Ignore non-matching lines."""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            s = line.strip()
            m = HDR_RE.match(s)
            if not m:
                # Normalize headers of the form '>::...' to '>'
                s_norm = re.sub(r"^>\s*:+", ">", s)
                m = HDR_RE.match(s_norm)
            if m:
                arm = m.group(1)
                start = int(m.group(2))
                end = int(m.group(3))

                # Optional filter: exclude records with coordinates larger than max_pos
                if max_pos is not None and (start > max_pos or end > max_pos):
                    continue

                if start > end:
                    start, end = end, start
                yield arm, start, end


def merge_ranges(records):
    """
    records: iterable of (arm, start, end)
    returns dict arm -> (min_pos, max_pos)
    """
    merged = {}
    for arm, start, end in records:
        if arm not in merged:
            merged[arm] = [start, end]
        else:
            merged[arm][0] = min(merged[arm][0], start)
            merged[arm][1] = max(merged[arm][1], end)
    return {arm: (v[0], v[1]) for arm, v in merged.items()}


def write_csv(merged, out_csv):
    """Write CSV with columns: arm, min_pos, max_pos, span."""
    def arm_key(a):
        """
        Sorting key: numeric order for chr1..chr22, then X,Y;
        within each chromosome, 'p' before 'q', then no suffix.
        """
        m = re.match(r"chr(\d+|X|Y)([pq]?)$", a)
        if not m:
            return (999, a)
        chrom = m.group(1)
        arm_suffix = m.group(2)
        if chrom == "X":
            chrom_num = 23
        elif chrom == "Y":
            chrom_num = 24
        else:
            chrom_num = int(chrom)
        arm_rank = {"p": 0, "q": 1, "": 2}.get(arm_suffix, 2)
        return (chrom_num, arm_rank)

    rows = []
    for arm, (mn, mx) in merged.items():
        span = mx - mn  # difference, not inclusive length
        rows.append((arm, mn, mx, span))
    rows.sort(key=lambda r: arm_key(r[0]))

    with open(out_csv, "w", newline="", encoding="utf-8") as w:
        writer = csv.writer(w)
        writer.writerow(["arm", "min_pos", "max_pos", "span"])
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(
        description="Extract and merge chromosome-arm ranges from FASTA headers into a CSV file."
    )
    ap.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA file.",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Path to output CSV file.",
    )
    ap.add_argument(
        "--max-pos",
        type=int,
        default=None,
        help=(
            "Maximum allowed coordinate. Headers whose start or end position "
            "exceeds this value will be excluded."
        ),
    )
    args = ap.parse_args()

    fasta_path = os.path.abspath(args.fasta)
    out_csv = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    recs = list(parse_headers(fasta_path, args.max_pos))
    merged = merge_ranges(recs)
    write_csv(merged, out_csv)
    print(f"[Done] Wrote {len(merged)} rows to {out_csv}")


if __name__ == "__main__":
    main()
