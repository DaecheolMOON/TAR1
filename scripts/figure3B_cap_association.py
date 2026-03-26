from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, binomtest

from figure3_common import (
    CHIMP_CAP_PLUS,
    CUTOFFS_DIST,
    build_universe,
    ensure_dir,
    get_tar1_fasta,
    normalize_arm_token,
    run_counting,
)


def parse_tar1_fasta_items(fasta_path):
    import re

    hdr_re = re.compile(r"^(chr[\w]+[pq]):(\d+)-(\d+)")
    items = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip().split()[0]
            m = hdr_re.match(header)
            if not m:
                continue
            arm = m.group(1)
            s = int(m.group(2))
            e = int(m.group(3))
            if s > e:
                s, e = e, s
            items.append((arm, s, e))
    return items


def build_tar1_counts_summary_chimp(data_root=None) -> pd.DataFrame:
    tar1_path = get_tar1_fasta("chimp", data_root=data_root)
    items = parse_tar1_fasta_items(tar1_path)

    rows = []
    for cutoff in CUTOFFS_DIST:
        found_arms = set()
        for arm, s, e in items:
            if s <= cutoff and e <= cutoff:
                found_arms.add(arm)

        arm_list = sorted(found_arms)
        rows.append({
            "cutoff": int(cutoff),
            "total_sequences": int(len(found_arms)),
            "unique_arms_with_sequence": int(len(found_arms)),
            "arms_list": ";".join(arm_list),
        })

    return pd.DataFrame(rows)


def parse_arms_list(cell) -> set:
    try:
        if pd.isna(cell):
            return set()
    except Exception:
        pass
    s = str(cell).strip()
    if s == "" or s.lower() in {"nan", "none", "null"}:
        return set()
    parts = [p for p in re_split_multi(s) if p]
    out = set()
    for p in parts:
        try:
            out.add(normalize_arm_token(p))
        except Exception:
            continue
    return out


def re_split_multi(s: str):
    import re
    return re.split(r"[;,\s]+", s)


def fisher_2x2_selected(selected: set, cap_plus: set, universe: set, alternative: str):
    a = len(selected & cap_plus)
    b = len(selected - cap_plus)
    c = len(cap_plus - selected)
    d = len(universe - (selected | cap_plus))
    table = [[a, b], [c, d]]
    oddsratio, p_value = fisher_exact(table, alternative=alternative)
    return a, b, c, d, float(oddsratio), float(p_value)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 3B chimp cap association plot")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    df = build_tar1_counts_summary_chimp(data_root=args.data_root)
    df["cutoff"] = pd.to_numeric(df["cutoff"], errors="coerce")
    df = df.dropna(subset=["cutoff"]).copy()
    df["cutoff"] = df["cutoff"].astype(int)
    df = df[(df["cutoff"] >= 30_000) & (df["cutoff"] <= 100_000_000)].copy()
    df = df.sort_values("cutoff").reset_index(drop=True)

    all_arms = build_universe()

    rows_out = []
    for _, row in df.iterrows():
        cutoff = int(row["cutoff"])
        tar1_set = parse_arms_list(row.get("arms_list", ""))
        if len(tar1_set) == 0:
            continue

        cap_minus_in = len(tar1_set - CHIMP_CAP_PLUS)
        n = len(tar1_set)
        binom = binomtest(k=cap_minus_in, n=n, p=0.5, alternative="greater")

        A, B, C, D, or_capplus, p_fisher = fisher_2x2_selected(
            selected=tar1_set,
            cap_plus=CHIMP_CAP_PLUS,
            universe=all_arms,
            alternative="less",
        )

        if np.isfinite(or_capplus) and or_capplus > 0:
            or_capminus = 1.0 / or_capplus
        elif or_capplus == 0:
            or_capminus = np.inf
        else:
            or_capminus = np.nan

        rows_out.append({
            "cutoff": cutoff,
            "tar1_plus_arms_n": n,
            "cap_plus_in_tar1_plus": A,
            "cap_minus_in_tar1_plus": B,
            "cap_plus_not_tar1_plus": C,
            "cap_minus_not_tar1_plus": D,
            "binom_p": float(binom.pvalue),
            "fisher_p": float(p_fisher),
            "odds_ratio_capplus": float(or_capplus),
            "odds_ratio_capminus": float(or_capminus),
            "cap_plus_frac_in_tar1_plus": A / n if n > 0 else np.nan,
        })

    out_df = pd.DataFrame(rows_out).sort_values("cutoff").reset_index(drop=True)
    out_csv = Path(outdir) / "Figure3B_chimp_cap_association.csv"
    out_df.to_csv(out_csv, index=False)

    tmp = out_df[["cutoff", "odds_ratio_capminus"]].copy()
    tmp["cutoff"] = pd.to_numeric(tmp["cutoff"], errors="coerce")
    tmp["odds_ratio_capminus"] = pd.to_numeric(tmp["odds_ratio_capminus"], errors="coerce")
    tmp = tmp.dropna().sort_values("cutoff")

    x = tmp["cutoff"].to_numpy(dtype=float)
    y = tmp["odds_ratio_capminus"].to_numpy(dtype=float)

    plt.figure(figsize=(8, 5))
    plt.plot(x, y, marker="o", linestyle="-", linewidth=2)

    plt.xscale("log")
    plt.xlim(30_000, 100_000_000)

    plt.axvspan(5e4, 1e5, alpha=0.30, label="50–100 kb window")
    plt.axvspan(5e5, 2e6, alpha=0.20, label="0.5–2 Mbp window")

    finite_y = y[np.isfinite(y)]
    ymax = float(np.max(finite_y)) if finite_y.size > 0 else 1.0
    plt.ylim(0, max(5, ymax * 1.05))

    plt.xlabel("Distance from telomere (bp, log scale)")
    plt.ylabel("Odds ratio (cap- enrichment)")
    plt.title("Figure 3B: chimp heterochromatin cap association")
    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}" if v >= 1 else f"{v:g}"))
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()

    out_png = Path(outdir) / "Figure3B_chimp_cap_association.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] csv -> {out_csv}")
    print(f"[OK] plot -> {out_png}")


if __name__ == "__main__":
    main()