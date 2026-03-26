from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from figure2_common import (
    STYLE_MAP,
    cumulative_to_bin_counts,
    ensure_dir,
    to_log_coord,
    unify_grid,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 2C normalized TAR1 density plot")
    parser.add_argument("--summary-csv", required=True, help="Output from figure2_distribution.py")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    df = pd.read_csv(args.summary_csv)

    species_bins = {}
    for sp in sorted(df["species"].unique()):
        d = df[df["species"] == sp].copy().sort_values("cutoff")
        x, nbin = cumulative_to_bin_counts(
            d["cutoff"].values,
            d["unique_arms_with_sequence"].values,
        )
        species_bins[sp] = (x, nbin)

    grid = unify_grid(species_bins)
    x_grid = grid["cutoff"].values.astype(float)
    z_grid = to_log_coord(x_grid, x0_ref=1.0)

    norm_df = pd.DataFrame(index=np.arange(len(x_grid)))
    norm_df["cutoff"] = x_grid
    norm_df["z_log_cutoff"] = z_grid

    for sp in [c for c in grid.columns if c != "cutoff"]:
        vec = grid[sp].values.astype(float)
        s = vec.sum()
        if s > 0:
            vec = vec / s
        norm_df[sp] = vec

    norm_csv = Path(outdir) / "Figure2C_normalized_density_table.csv"
    norm_df.to_csv(norm_csv, index=False)

    plt.figure(figsize=(10, 6))
    for sp in [c for c in norm_df.columns if c not in {"cutoff", "z_log_cutoff"}]:
        st = STYLE_MAP.get(sp, {"label": sp})
        plt.plot(
            norm_df["z_log_cutoff"].values,
            norm_df[sp].values,
            linewidth=2,
            label=st["label"],
        )

    plt.xlabel("z = log(distance from chromosome end)")
    plt.ylabel("Normalized density")
    plt.legend()
    plt.tight_layout()

    out_png = Path(outdir) / "Figure2C_normalized_density.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] table -> {norm_csv}")
    print(f"[OK] plot -> {out_png}")


if __name__ == "__main__":
    import numpy as np
    main()