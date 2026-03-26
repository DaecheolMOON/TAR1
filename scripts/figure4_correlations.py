from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from figure4_common import (
    build_total_human_tar1_vs_telomere_table,
    ensure_dir,
    get_delta_a_csv,
    load_schmidt_meta,
)


def safe_pearson(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(tmp) < 3 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, len(tmp)
    r, p = stats.pearsonr(tmp["x"], tmp["y"])
    return float(r), float(p), len(tmp)


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(tmp) < 3 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, len(tmp)
    r, p = stats.spearmanr(tmp["x"], tmp["y"])
    return float(r), float(p), len(tmp)


def build_doubling_df() -> pd.DataFrame:
    return pd.DataFrame({
        "Sample_Name": ["Calu-3", "G-292", "HOS", "HT-1080", "HT-29", "SK-LU-1", "SK-N-AS", "Saos-2", "U2-OS", "SK-N-F1"],
        "Doubling_Time": [105.83, 102.19, 41.06, 78.31, 28.35, 98.21, 46.57, 99.13, 71.88, 97.22],
    })


def panel_a_total_and_cancer(data_root: str, outdir: Path) -> None:
    df = build_total_human_tar1_vs_telomere_table(data_root=data_root)
    df.to_csv(outdir / "Figure4A_total_human_tar1_vs_telomere.csv", index=False)

    df_all = df.dropna(subset=["ground_truth_bps", "tar1_median_bps"]).copy()
    df_cancer = df_all[df_all["sample_group"] == "Cancer"].copy()

    r_all, p_all, n_all = safe_pearson(df_all["tar1_median_bps"], df_all["ground_truth_bps"])
    r_c, p_c, n_c = safe_pearson(df_cancer["tar1_median_bps"], df_cancer["ground_truth_bps"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.regplot(data=df_all, x="tar1_median_bps", y="ground_truth_bps", ax=axes[0],
                scatter_kws={"s": 50, "alpha": 0.7}, line_kws={"linewidth": 2})
    axes[0].set_title(f"Total human samples\nn={n_all}, Pearson r={r_all:.3f}, p={p_all:.3g}")
    axes[0].set_xlabel("TAR1 median position (bp)")
    axes[0].set_ylabel("Mean telomere length (bp)")
    axes[0].grid(True, linestyle="--", alpha=0.5)

    sns.regplot(data=df_cancer, x="tar1_median_bps", y="ground_truth_bps", ax=axes[1],
                scatter_kws={"s": 50, "alpha": 0.7}, line_kws={"linewidth": 2}, color="#d62728")
    axes[1].set_title(f"Cancer cell lines\nn={n_c}, Pearson r={r_c:.3f}, p={p_c:.3g}")
    axes[1].set_xlabel("TAR1 median position (bp)")
    axes[1].set_ylabel("Mean telomere length (bp)")
    axes[1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(outdir / "Figure4A_total_and_cancer_correlation.png", dpi=300, bbox_inches="tight")
    plt.close()


def panel_g_and_h(data_root: str, outdir: Path) -> None:
    delta_df = pd.read_csv(get_delta_a_csv(data_root))
    schmidt = load_schmidt_meta(data_root)
    dt_df = build_doubling_df()

    df = pd.merge(delta_df, schmidt, left_on="sample", right_on="Run_ID", how="inner")
    df = pd.merge(df, dt_df, on="Sample_Name", how="inner")

    for col in ["delta_a", "Lower_quartile_bp", "Doubling_Time"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df.to_csv(outdir / "Figure4GH_merged_deltaa_schmidt_dt.csv", index=False)

    df_no_hos = df[df["Sample_Name"] != "HOS"].copy()
    df_g = df_no_hos.groupby("Sample_Name").mean(numeric_only=True).reset_index()
    r_g, p_g, n_g = safe_pearson(df_g["Lower_quartile_bp"], df_g["Doubling_Time"])

    df_h = df.groupby("Sample_Name").mean(numeric_only=True).reset_index()
    r_h, p_h, n_h = safe_spearman(df_h["delta_a"], df_h["Lower_quartile_bp"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.regplot(data=df_g, x="Lower_quartile_bp", y="Doubling_Time", ax=axes[0],
                scatter_kws={"s": 90, "alpha": 0.8, "edgecolor": "k"},
                line_kws={"linewidth": 2})
    for _, row in df_g.iterrows():
        axes[0].text(row["Lower_quartile_bp"], row["Doubling_Time"], str(row["Sample_Name"]), fontsize=9)
    axes[0].set_title(f"Figure 4G\nn={n_g}, Pearson r={r_g:.3f}, p={p_g:.3g}")
    axes[0].set_xlabel("Lower-quartile telomere length (bp)")
    axes[0].set_ylabel("Doubling time (hours)")
    axes[0].grid(True, linestyle="--", alpha=0.5)

    sns.regplot(data=df_h, x="delta_a", y="Lower_quartile_bp", ax=axes[1],
                scatter_kws={"s": 90, "alpha": 0.8, "edgecolor": "k"},
                line_kws={"linewidth": 2}, color="#2ca02c")
    for _, row in df_h.iterrows():
        axes[1].text(row["delta_a"], row["Lower_quartile_bp"], str(row["Sample_Name"]), fontsize=9)
    axes[1].set_title(f"Figure 4H\nn={n_h}, Spearman rho={r_h:.3f}, p={p_h:.3g}")
    axes[1].set_xlabel("delta_a (bp)")
    axes[1].set_ylabel("Lower-quartile telomere length (bp)")
    axes[1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(outdir / "Figure4G_H_correlations.png", dpi=300, bbox_inches="tight")
    plt.close()

    summary = pd.DataFrame([
        {"panel": "G", "n": n_g, "method": "Pearson", "x": "Lower_quartile_bp", "y": "Doubling_Time", "r": r_g, "p": p_g, "exclude": "HOS"},
        {"panel": "H", "n": n_h, "method": "Spearman", "x": "delta_a", "y": "Lower_quartile_bp", "r": r_h, "p": p_h, "exclude": ""},
    ])
    summary.to_csv(outdir / "Figure4GH_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 4 correlation panels")
    parser.add_argument("--data-root", default="./data")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    panel_a_total_and_cancer(args.data_root, outdir)
    panel_g_and_h(args.data_root, outdir)

    print(f"[OK] correlation panels saved under -> {outdir}")


if __name__ == "__main__":
    main()