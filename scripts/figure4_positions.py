from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure4_common import (
    ensure_dir,
    build_final_position_table,
    get_cluster_master_csv,
    load_donor_meta,
)

# =========================================================
# Local constants
# =========================================================

FINAL_GROUP_ORDER = [
    "Cancer (Fast Growth)",
    "Cancer (Slow Growth)",
    "Short Telomere Disease",
    "Age 0",
    "Age 1-59",
    "Age 60+",
]

GROUP_ORDER_RATIO = [
    "Newborn (Age 0)",
    "Adult (Age 1+)",
    "Cancer",
]

AGE_GROUP_ORDER = [
    "Age 0",
    "Age 1-59",
    "Age 60+",
]

DENSITY_GROUP_ORDER = [
    "Age 0",
    "Age 1-59",
    "Age 60+",
    "Cancer (SRR)",
]

COLOR_DICT = {
    "Cancer (Fast Growth)": "#6baed6",
    "Cancer (Slow Growth)": "#2171b5",
    "Short Telomere Disease": "#2ca25f",
    "Age 0": "#cb181d",
    "Age 1-59": "#fb6a4a",
    "Age 60+": "#fcae91",
}

MARKER_DICT = {
    "TERT": "s",
    "ALT": "^",
    "Normal": "o",
    "Unknown": "X",
}

TMM_ORDER = ["TERT", "ALT"]

MIN_POS_K = 500
MAX_POS_K = 15000


# =========================================================
# Helpers
# =========================================================

def get_star_from_pvalue(p: float) -> str:
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def mannwhitney_safe(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    from scipy.stats import mannwhitneyu

    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()

    if len(a) < 1 or len(b) < 1:
        return np.nan, np.nan

    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(stat), float(p)


def scatter_with_jitter(ax, x_center: float, y: pd.Series, width: float = 0.14, size: float = 36.0) -> None:
    vals = pd.to_numeric(y, errors="coerce").dropna().values
    if len(vals) == 0:
        return
    rng = np.random.default_rng(42)
    xs = x_center + rng.uniform(-width, width, size=len(vals))
    ax.scatter(xs, vals, s=size, color="black", alpha=0.6, zorder=3)


def style_boxplot(bp, facecolors: list[str]) -> None:
    for patch, fc in zip(bp["boxes"], facecolors):
        patch.set_facecolor(fc)
        patch.set_alpha(0.65)
        patch.set_linewidth(2.0)
    for med in bp["medians"]:
        med.set_color("black")
        med.set_linewidth(2.0)
    for whisk in bp["whiskers"]:
        whisk.set_linewidth(1.8)
    for cap in bp["caps"]:
        cap.set_linewidth(1.8)


def map_ratio_group(final_group: str) -> str | None:
    if final_group == "Age 0":
        return "Newborn (Age 0)"
    if final_group in {"Age 1-59", "Age 60+"}:
        return "Adult (Age 1+)"
    if isinstance(final_group, str) and "Cancer" in final_group:
        return "Cancer"
    return None


def map_density_group(sample_id: str, donor_age_map: dict[str, float]) -> str | None:
    sid = str(sample_id).strip()

    if sid.startswith("SRR2684") or sid.startswith("SRR"):
        return "Cancer (SRR)"

    donor_id = sid.split(".", 1)[0] if "." in sid else sid
    age = donor_age_map.get(donor_id, np.nan)
    if pd.isna(age):
        return None

    age = int(age)
    if age == 0:
        return "Age 0"
    if 1 <= age <= 59:
        return "Age 1-59"
    if age >= 60:
        return "Age 60+"
    return None


def parse_interval_from_header(header: str) -> tuple[int | None, int | None]:
    if header is None or (isinstance(header, float) and np.isnan(header)):
        return None, None
    m = re.search(r"(\d+)-(\d+)\s*$", str(header).strip())
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


# =========================================================
# Panel B
# =========================================================

def build_panel_b(df_final: pd.DataFrame, outdir: Path) -> None:
    outdir = ensure_dir(outdir)

    fig, ax = plt.subplots(figsize=(8.5, 7.5))

    used_groups = [g for g in FINAL_GROUP_ORDER if g in df_final["Final_Group"].dropna().unique()]

    for fg in used_groups:
        sub = df_final[df_final["Final_Group"] == fg].copy()
        for _, row in sub.iterrows():
            mk = MARKER_DICT.get(str(row["Shape_Type"]), "o")
            color = COLOR_DICT.get(fg, "gray")
            ax.scatter(
                float(row["Cluster_0_Pos_Median"]),
                float(row["Cluster_1_Pos_Median"]),
                s=120,
                marker=mk,
                color=color,
                edgecolors="black",
                linewidths=0.8,
                alpha=0.85,
                zorder=3,
            )

    x_min = float(pd.concat([
        pd.to_numeric(df_final["Cluster_0_Pos_Median"], errors="coerce"),
        pd.to_numeric(df_final["Cluster_1_Pos_Median"], errors="coerce"),
    ]).min())
    x_max = float(pd.concat([
        pd.to_numeric(df_final["Cluster_0_Pos_Median"], errors="coerce"),
        pd.to_numeric(df_final["Cluster_1_Pos_Median"], errors="coerce"),
    ]).max())

    ax.plot([x_min, x_max], [x_min, x_max], linestyle="--", color="black", alpha=0.6, zorder=1)

    ax.set_title("Figure 4B: Cluster 0 vs Cluster 1 positions")
    ax.set_xlabel("Cluster 0 median position (bp)")
    ax.set_ylabel("Cluster 1 median position (bp)")
    ax.grid(True, linestyle="--", alpha=0.35)

    out_png = Path(outdir) / "Figure4B_cluster0_vs_cluster1.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    df_final.to_csv(Path(outdir) / "Figure4B_source_table.csv", index=False)
    print(f"[OK] Figure4B -> {out_png}")


# =========================================================
# Panel C
# =========================================================

def build_panel_c(df_final: pd.DataFrame, outdir: Path) -> None:
    outdir = ensure_dir(outdir)

    df_ratio = df_final.copy()
    df_ratio["Ratio_Group"] = df_ratio["Final_Group"].map(map_ratio_group)
    df_ratio = df_ratio[df_ratio["Ratio_Group"].isin(GROUP_ORDER_RATIO)].copy()
    df_ratio["Pos_Ratio"] = pd.to_numeric(df_ratio["Cluster_1_Pos_Median"], errors="coerce") / pd.to_numeric(df_ratio["Cluster_0_Pos_Median"], errors="coerce")
    df_ratio = df_ratio.dropna(subset=["Pos_Ratio"])

    fig, ax = plt.subplots(figsize=(8, 7))

    data_list = [
        df_ratio.loc[df_ratio["Ratio_Group"] == g, "Pos_Ratio"].astype(float).values
        for g in GROUP_ORDER_RATIO
    ]
    facecolors = ["#cb181d", "#3182bd", "#4d4d4d"]

    bp = ax.boxplot(data_list, tick_labels=GROUP_ORDER_RATIO, showfliers=False, patch_artist=True)
    style_boxplot(bp, facecolors)

    for i, g in enumerate(GROUP_ORDER_RATIO, start=1):
        scatter_with_jitter(ax, i, df_ratio.loc[df_ratio["Ratio_Group"] == g, "Pos_Ratio"])

    pairs = [
        ("Newborn (Age 0)", "Adult (Age 1+)"),
        ("Adult (Age 1+)", "Cancer"),
        ("Newborn (Age 0)", "Cancer"),
    ]

    y_max = float(df_ratio["Pos_Ratio"].max())
    y_min = float(df_ratio["Pos_Ratio"].min())
    y_range = max(0.05, y_max - y_min)
    base = y_max + y_range * 0.06

    for idx, (g1, g2) in enumerate(pairs):
        _, p = mannwhitney_safe(
            df_ratio.loc[df_ratio["Ratio_Group"] == g1, "Pos_Ratio"],
            df_ratio.loc[df_ratio["Ratio_Group"] == g2, "Pos_Ratio"],
        )
        x1 = GROUP_ORDER_RATIO.index(g1) + 1
        x2 = GROUP_ORDER_RATIO.index(g2) + 1
        y = base + idx * y_range * 0.10
        h = y_range * 0.02
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", linewidth=1.8)
        ax.text((x1 + x2) / 2.0, y + h, get_star_from_pvalue(p), ha="center", va="bottom", fontsize=14)

    ax.set_title("Figure 4C: Cluster 1 / Cluster 0 ratio")
    ax.set_ylabel("Position ratio")
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)

    out_png = Path(outdir) / "Figure4C_ratio_boxplot.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    df_ratio.to_csv(Path(outdir) / "Figure4C_ratio_source_table.csv", index=False)
    print(f"[OK] Figure4C -> {out_png}")


# =========================================================
# Panel D
# =========================================================

def build_panel_d(df_final: pd.DataFrame, outdir: Path) -> None:
    outdir = ensure_dir(outdir)

    df_tmm = df_final[
        df_final["Final_Group"].astype(str).str.contains("Cancer", na=False)
        & df_final["Shape_Type"].isin(TMM_ORDER)
    ].copy()

    fig, ax = plt.subplots(figsize=(6.5, 7.0))

    data_list = [
        df_tmm.loc[df_tmm["Shape_Type"] == g, "Cluster_0_Pos_Median"].astype(float).values
        for g in TMM_ORDER
    ]
    facecolors = ["#e67e22", "#8e44ad"]

    bp = ax.boxplot(data_list, tick_labels=TMM_ORDER, showfliers=False, patch_artist=True)
    style_boxplot(bp, facecolors)

    for i, g in enumerate(TMM_ORDER, start=1):
        scatter_with_jitter(ax, i, df_tmm.loc[df_tmm["Shape_Type"] == g, "Cluster_0_Pos_Median"])

    _, p = mannwhitney_safe(
        df_tmm.loc[df_tmm["Shape_Type"] == "TERT", "Cluster_0_Pos_Median"],
        df_tmm.loc[df_tmm["Shape_Type"] == "ALT", "Cluster_0_Pos_Median"],
    )

    if len(df_tmm) > 0:
        y_max = float(pd.to_numeric(df_tmm["Cluster_0_Pos_Median"], errors="coerce").max())
        y_min = float(pd.to_numeric(df_tmm["Cluster_0_Pos_Median"], errors="coerce").min())
        y_range = max(1.0, y_max - y_min)
        y = y_max + y_range * 0.08
        h = y_range * 0.02
        ax.plot([1, 1, 2, 2], [y, y + h, y + h, y], color="black", linewidth=1.8)
        ax.text(1.5, y + h, get_star_from_pvalue(p), ha="center", va="bottom", fontsize=14)

    ax.set_title("Figure 4D: Cluster 0 position (TERT vs ALT)")
    ax.set_ylabel("Cluster 0 median position (bp)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)

    out_png = Path(outdir) / "Figure4D_TERT_vs_ALT_cluster0.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    df_tmm.to_csv(Path(outdir) / "Figure4D_TERT_vs_ALT_source_table.csv", index=False)
    print(f"[OK] Figure4D -> {out_png}")


# =========================================================
# Panels E/F
# =========================================================

def build_age_boxplot(df_final: pd.DataFrame, outdir: Path, cluster_col: str, panel_label: str) -> None:
    outdir = ensure_dir(outdir)

    df_age = df_final[df_final["Final_Group"].isin(AGE_GROUP_ORDER)].copy()
    df_age = df_age.dropna(subset=[cluster_col])

    fig, ax = plt.subplots(figsize=(7.5, 7.0))

    data_list = [
        df_age.loc[df_age["Final_Group"] == g, cluster_col].astype(float).values
        for g in AGE_GROUP_ORDER
    ]
    facecolors = ["#cb181d", "#fb6a4a", "#fcae91"]

    bp = ax.boxplot(data_list, tick_labels=AGE_GROUP_ORDER, showfliers=False, patch_artist=True)
    style_boxplot(bp, facecolors)

    for i, g in enumerate(AGE_GROUP_ORDER, start=1):
        scatter_with_jitter(ax, i, df_age.loc[df_age["Final_Group"] == g, cluster_col])

    pairs = [
        ("Age 0", "Age 1-59"),
        ("Age 1-59", "Age 60+"),
        ("Age 0", "Age 60+"),
    ]

    if len(df_age) > 0:
        y_max = float(pd.to_numeric(df_age[cluster_col], errors="coerce").max())
        y_min = float(pd.to_numeric(df_age[cluster_col], errors="coerce").min())
        y_range = max(1.0, y_max - y_min)
        base = y_max + y_range * 0.08

        idx_sig = 0
        for g1, g2 in pairs:
            _, p = mannwhitney_safe(
                df_age.loc[df_age["Final_Group"] == g1, cluster_col],
                df_age.loc[df_age["Final_Group"] == g2, cluster_col],
            )
            if pd.isna(p) or p >= 0.05:
                continue
            x1 = AGE_GROUP_ORDER.index(g1) + 1
            x2 = AGE_GROUP_ORDER.index(g2) + 1
            y = base + idx_sig * y_range * 0.10
            h = y_range * 0.02
            ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", linewidth=1.8)
            ax.text((x1 + x2) / 2.0, y + h, get_star_from_pvalue(p), ha="center", va="bottom", fontsize=14)
            idx_sig += 1

    ax.set_title(f"{panel_label}: {cluster_col}")
    ax.set_ylabel("Median position (bp)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)

    out_png = Path(outdir) / f"{panel_label}_{cluster_col}.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    df_age.to_csv(Path(outdir) / f"{panel_label}_{cluster_col}_source_table.csv", index=False)
    print(f"[OK] {panel_label} -> {out_png}")


# =========================================================
# Panel K
# Original-dump-like logic:
# use cluster master + header parsing + coverage accumulation
# =========================================================

def build_density_panel_k(data_root: str | None, outdir: Path) -> None:
    outdir = ensure_dir(outdir)

    cluster_csv = get_cluster_master_csv(data_root)
    df = pd.read_csv(cluster_csv)

    donor_meta = load_donor_meta(data_root)
    donor_meta["donor_id"] = donor_meta["donor_id"].astype(str).str.strip()
    donor_meta["age_years"] = pd.to_numeric(donor_meta["age_years"], errors="coerce")
    donor_age_map = dict(zip(donor_meta["donor_id"], donor_meta["age_years"]))

    group_coverage = {
        g: {
            0: np.zeros(MAX_POS_K + 1, dtype=int),
            1: np.zeros(MAX_POS_K + 1, dtype=int),
        }
        for g in DENSITY_GROUP_ORDER
    }
    group_counts = {g: 0 for g in DENSITY_GROUP_ORDER}

    for _, row in df.iterrows():
        sample_id = str(row["Sample_ID"]).strip()
        cluster_code = int(row["Cluster_Code"])
        full_header = row["Full_Header"]

        group_name = map_density_group(sample_id, donor_age_map)
        if group_name is None:
            continue

        start, end = parse_interval_from_header(full_header)
        if start is None or end is None:
            continue

        group_counts[group_name] += 1

        s = min(start, end)
        e = max(start, end)
        if s > MAX_POS_K:
            continue

        s_clamp = max(s, MIN_POS_K)
        e_clamp = min(e, MAX_POS_K)
        if s_clamp <= e_clamp and cluster_code in [0, 1]:
            group_coverage[group_name][cluster_code][s_clamp:e_clamp + 1] += 1

    x_positions = np.arange(MIN_POS_K, MAX_POS_K + 1)

    fig, axes = plt.subplots(1, 4, figsize=(22, 5.6), sharey=False)

    density_rows = []

    for i, group_name in enumerate(DENSITY_GROUP_ORDER):
        ax = axes[i]

        y_c0 = group_coverage[group_name][0][MIN_POS_K:MAX_POS_K + 1].astype(float)
        y_c1 = group_coverage[group_name][1][MIN_POS_K:MAX_POS_K + 1].astype(float)

        sum_c0 = y_c0.sum()
        sum_c1 = y_c1.sum()

        norm_c0 = y_c0 / sum_c0 if sum_c0 > 0 else y_c0
        norm_c1 = y_c1 / sum_c1 if sum_c1 > 0 else y_c1

        ax.plot(x_positions, norm_c0, linewidth=1.6, label=f"Cluster 0 (n={int(sum_c0):,})")
        ax.fill_between(x_positions, 0, norm_c0, alpha=0.10)

        ax.plot(x_positions, norm_c1, linewidth=1.6, label=f"Cluster 1 (n={int(sum_c1):,})")
        ax.fill_between(x_positions, 0, norm_c1, alpha=0.10)

        ax.set_title(group_name)
        ax.set_xlabel("Position (bp)")
        if i == 0:
            ax.set_ylabel("Normalized density")
        ax.grid(True, linestyle=":", alpha=0.45)
        ax.set_xlim(MIN_POS_K, MAX_POS_K)
        ax.legend(loc="upper right", fontsize=8)

        for pos, v0, v1 in zip(x_positions, norm_c0, norm_c1):
            density_rows.append({
                "group": group_name,
                "position_bp": int(pos),
                "cluster0_density": float(v0),
                "cluster1_density": float(v1),
                "cluster0_total": int(sum_c0),
                "cluster1_total": int(sum_c1),
                "n_sequences_in_group": int(group_counts[group_name]),
            })

    out_png = Path(outdir) / "Figure4K_cluster_density_profiles.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(density_rows).to_csv(Path(outdir) / "Figure4K_cluster_density_profiles.csv", index=False)
    print(f"[OK] Figure4K -> {out_png}")


# =========================================================
# Main
# =========================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 4 position-based panels")
    parser.add_argument("--data-root", default="./data", help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    df_final = build_final_position_table(data_root=args.data_root)
    df_final_csv = Path(outdir) / "Figure4_final_position_table.csv"
    df_final.to_csv(df_final_csv, index=False)
    print(f"[OK] final position table -> {df_final_csv}")

    build_panel_b(df_final, Path(outdir) / "panelB")
    build_panel_c(df_final, Path(outdir) / "panelC")
    build_panel_d(df_final, Path(outdir) / "panelD")
    build_age_boxplot(df_final, Path(outdir) / "panelE", "Cluster_0_Pos_Median", "Figure4E")
    build_age_boxplot(df_final, Path(outdir) / "panelF", "Cluster_1_Pos_Median", "Figure4F")
    build_density_panel_k(args.data_root, Path(outdir) / "panelK")


if __name__ == "__main__":
    main()