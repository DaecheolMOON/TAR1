# =========================================================
# Cell 10: summary_runs.csv -> single PDF signal plot
# Saves next to summary_runs.csv:
#   signal_plot_{species}_{ref_cutoff}_len{len_bp}.pdf
# =========================================================

import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt

_required = ["DIR_RESULT_PLOT", "SELECTED_SPECIES", "ANALYSIS_BP"]
_missing = [k for k in _required if k not in globals()]
if _missing:
    raise RuntimeError(f"Missing globals: {_missing}")

_sp = str(SELECTED_SPECIES)
_len_bp = int(ANALYSIS_BP)
_ref_cutoff = int(REF_CUTOFF) if "REF_CUTOFF" in globals() else None

def resolve_summary_runs_csv():
    if "summary_runs_csv" in globals():
        p = globals()["summary_runs_csv"]
        if isinstance(p, str) and p and os.path.exists(p):
            return os.path.abspath(p)

    ref_cwt_root = os.path.join(DIR_RESULT_PLOT, "reference_cwt")
    if _ref_cutoff is not None:
        expected_dir = os.path.join(ref_cwt_root, f"{_sp}_{_ref_cutoff}_len{_len_bp}")
        expected_csv = os.path.join(expected_dir, "summary_runs.csv")
        if os.path.exists(expected_csv):
            return os.path.abspath(expected_csv)

    pattern = os.path.join(ref_cwt_root, "**", "summary_runs.csv")
    hits = glob.glob(pattern, recursive=True)
    if not hits:
        raise FileNotFoundError(f"summary_runs.csv not found under: {ref_cwt_root}")

    def score(path):
        s = 0
        p = path.lower()
        if _sp.lower() in p:
            s += 10
        if f"len{_len_bp}" in p:
            s += 10
        if _ref_cutoff is not None and str(_ref_cutoff) in p:
            s += 5
        return s

    hits_sorted = sorted(hits, key=lambda x: score(x), reverse=True)
    return os.path.abspath(hits_sorted[0])

INPUT_CSV = resolve_summary_runs_csv()
IN_DIR = os.path.dirname(INPUT_CSV)

out_pdf = os.path.join(
    IN_DIR,
    f"signal_plot_{_sp}" + (f"_{_ref_cutoff}" if _ref_cutoff is not None else "") + f"_len{_len_bp}.pdf"
)

print("[INFO] INPUT_CSV:", INPUT_CSV)
print("[INFO] OUT_PDF  :", out_pdf)

def plot_signals_from_summary_runs(df: pd.DataFrame, total_length: int, output_path: str):
    required_cols = {"record_id", "start_base", "end_base", "mean_atcg"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    def arm_key(a):
        a = str(a).strip()
        a2 = a.lower()
        if not a2.startswith("chr"):
            a2 = "chr" + a2
        m = re.match(r"^chr(\d+|x|y)([pq])$", a2)
        if not m:
            return (999, 9, a)
        c = m.group(1)
        pq = m.group(2)
        cn = 23 if c == "x" else (24 if c == "y" else int(c))
        pq_rank = 0 if pq == "p" else 1
        return (cn, pq_rank, a)

    arms = sorted(df["record_id"].dropna().unique().tolist(), key=arm_key)
    if not arms:
        raise ValueError("No arms found in record_id")

    df = df.copy()
    for c in ["start_base", "end_base", "mean_atcg"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if total_length is None or total_length <= 0:
        total_length = int(df["end_base"].max(skipna=True))

    n_arms = len(arms)
    fig, axes = plt.subplots(n_arms, 1, figsize=(12, 1.3 * n_arms), sharex=True)
    if n_arms == 1:
        axes = [axes]

    for ax, arm in zip(axes, arms):
        arm_df = df[df["record_id"] == arm]
        for _, row in arm_df.iterrows():
            start = row["start_base"]
            end = row["end_base"]
            if pd.isna(start) or pd.isna(end):
                continue
            start = max(0, min(int(start), total_length))
            end = max(0, min(int(end), total_length))
            if end < start:
                start, end = end, start

            height = row["mean_atcg"]
            if pd.isna(height):
                height = 0.0
            height = float(height)

            ax.fill_between([start, end], [0, 0], [height, height], color="black", linewidth=0)

        ax.set_xlim(0, total_length)
        ax.set_ylim(0, 4)
        ax.set_ylabel(str(arm), rotation=0, labelpad=30, fontsize=8)
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    plt.xlabel("Genomic position within analysis window (bp)")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
    print("[OK] Saved:", output_path)

df = pd.read_csv(INPUT_CSV)
plot_signals_from_summary_runs(df, total_length=_len_bp, output_path=out_pdf)
