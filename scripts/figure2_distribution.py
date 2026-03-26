from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

from figure2_common import (
    CUTOFFS_DIST,
    STYLE_MAP,
    TOTAL_ARMS,
    ensure_dir,
    get_tar1_fasta,
    run_counting,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 2B cumulative TAR1 distribution plot")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--data-root", default="./data", help="Root directory for data")
    parser.add_argument(
        "--species",
        nargs="*",
        default=["human", "chimp", "bonobo", "gorilla", "borang", "sorang"],
        help="Species list for plotting",
    )
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    results_rows = []
    results_counts = {}

    for species in args.species:
        tar1_file = get_tar1_fasta(species)
        counts, num_items = run_counting(tar1_file, CUTOFFS_DIST)
        results_counts[species] = counts
        print(f"[{species.upper()}] Parsed {num_items} blocks -> counts OK")

        denom = TOTAL_ARMS.get(species, 48)
        for cutoff, count in zip(CUTOFFS_DIST, counts):
            results_rows.append({
                "species": species,
                "cutoff": int(cutoff),
                "unique_arms_with_sequence": int(count),
                "total_arms": int(denom),
                "fraction": float(count) / float(denom),
                "percent": 100.0 * float(count) / float(denom),
            })

    df = pd.DataFrame(results_rows)
    summary_csv = Path(outdir) / "tar1_counts_summary_all_species.csv"
    df.to_csv(summary_csv, index=False)

    plt.figure(figsize=(10, 6))
    for species, counts in results_counts.items():
        if species not in STYLE_MAP:
            continue
        st = STYLE_MAP[species]
        denom = TOTAL_ARMS.get(species, 48)
        perc = [100.0 * v / denom for v in counts]

        plt.plot(
            CUTOFFS_DIST,
            perc,
            color=st["color"],
            marker=st["marker"],
            linestyle=st["style"],
            linewidth=2,
            label=st["label"],
        )

    plt.xscale("log")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:g}"))
    plt.ylim(0, 105)
    plt.xlabel("Distance from chromosome end (bp, log scale)")
    plt.ylabel("Cumulative TAR1-positive chromosome arms (%)")
    plt.legend()
    plt.tight_layout()

    out_png = Path(outdir) / "Figure2B_cumulative_distribution.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] summary -> {summary_csv}")
    print(f"[OK] plot -> {out_png}")


if __name__ == "__main__":
    main()