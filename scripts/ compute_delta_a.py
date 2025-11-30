#!/usr/bin/env python3
import argparse
import os
import re
import csv
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute delta_a from Mode 2 peaks in Y002_430/mode2_peak_summary.txt "
            "under sample subdirectories (e.g. JH*)."
        )
    )
    parser.add_argument(
        "base_dir",
        help="Base directory containing sample subdirectories (e.g. JH*).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output CSV file path (default: base_dir/delta_a_summary.csv)",
        default=None,
    )
    parser.add_argument(
        "--sample-prefix",
        default="JH",
        help=(
            "Prefix of sample directories inside base_dir (default: JH). "
            "Only directories whose names start with this prefix are processed."
        ),
    )
    parser.add_argument(
        "--single-peak-max",
        type=int,
        default=800,
        help=(
            "Reference peak position (bp) used when only one Mode 2 peak is found "
            "(delta_a = single_peak_max - first_peak). Default: 800."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = args.base_dir
    sample_prefix = args.sample_prefix
    single_peak_max = args.single_peak_max

    if not os.path.isdir(base_dir):
        sys.exit(f"Error: Base directory not found: {base_dir}")

    output_path = args.output or os.path.join(base_dir, "delta_a_summary.csv")

    # Regex patterns for two-peak and single-peak cases
    pattern_two = re.compile(r"Mode 2 peaks \(bp\):\s*(\d+),\s*(\d+)")
    pattern_one = re.compile(r"Mode 2 peaks \(bp\):\s*(\d+)\b")

    results = []

    for entry in sorted(os.listdir(base_dir)):
        dir_path = os.path.join(base_dir, entry)
        # Process only sample directories with the specified prefix
        if not (os.path.isdir(dir_path) and entry.startswith(sample_prefix)):
            continue

        summary_file = os.path.join(dir_path, "Y002_430", "mode2_peak_summary.txt")
        if not os.path.isfile(summary_file):
            print(f"Warning: Not found {summary_file}", file=sys.stderr)
            continue

        with open(summary_file, "r") as f:
            text = f.read()

        match_two = pattern_two.search(text)
        if match_two:
            first, second = map(int, match_two.groups())
            delta_a = second - first
        else:
            match_one = pattern_one.search(text)
            if match_one:
                first = int(match_one.group(1))
                delta_a = single_peak_max - first
            else:
                print(
                    f"Warning: Pattern not found in {summary_file}", file=sys.stderr
                )
                continue

        results.append({"sample": entry, "delta_a": delta_a})

    if not results:
        print("No valid data found; exiting.", file=sys.stderr)
        sys.exit(1)

    with open(output_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["sample", "delta_a"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"Saved CSV to {output_path}")


if __name__ == "__main__":
    main()
