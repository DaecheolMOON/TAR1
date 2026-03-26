from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure4_common import (
    compute_vmode_profile,
    ensure_dir,
    get_vmode_fasta,
)


def draw_v2_panel(species_key: str, title: str, out_png: Path, out_csv: Path, data_root: str | None) -> None:
    scales = np.arange(200, 401)
    fasta = get_vmode_fasta(species_key, data_root=data_root)
    res = compute_vmode_profile(fasta_path=fasta, cwt_scales=scales, kmer_k=20)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(res["scales"], res["v2"], marker="o", markersize=3, linewidth=1.5, color="black")

    peak_text = ""
    if len(res["major_peaks"]) >= 1:
        for p in res["major_peaks"]:
            ax.axvline(float(p), linestyle="--", color="red", alpha=0.7)
        peak_text = ", ".join(str(int(x)) for x in res["major_peaks"])

    delta_a = res["delta_a"]
    title_line = title
    if pd.notna(delta_a):
        title_line += f"\nmajor peaks: {peak_text} | delta_a={int(delta_a)} bp"

    ax.set_title(title_line)
    ax.set_xlabel("Scale (bp)")
    ax.set_ylabel("V2 loading")
    ax.grid(True, linestyle=":", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    out_df = pd.DataFrame({
        "scale_bp": res["scales"],
        "v2_loading": res["v2"],
    })
    out_df.to_csv(out_csv, index=False)

    summary_df = pd.DataFrame([{
        "species_key": species_key,
        "fasta_path": str(fasta),
        "reads_n": int(res["reads_n"]),
        "major_peaks_bp": ",".join(str(int(x)) for x in res["major_peaks"]) if len(res["major_peaks"]) > 0 else "",
        "delta_a_bp": float(delta_a) if pd.notna(delta_a) else np.nan,
        "v2_explained_ratio": float(res["var_ratio"][1]) if len(res["var_ratio"]) > 1 else float(res["var_ratio"][0]),
    }])
    summary_df.to_csv(out_png.with_suffix(".summary.csv"), index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 4I-J V2 profiles")
    parser.add_argument("--data-root", default="./data", help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    draw_v2_panel(
        species_key="human_chm13",
        title="Figure 4I: CHM13 V2 scale-mode profile",
        out_png=Path(outdir) / "Figure4I_CHM13_V2_profile.png",
        out_csv=Path(outdir) / "Figure4I_CHM13_V2_profile.csv",
        data_root=args.data_root,
    )

    draw_v2_panel(
        species_key="hg002",
        title="Figure 4J: HG002 V2 scale-mode profile",
        out_png=Path(outdir) / "Figure4J_HG002_V2_profile.png",
        out_csv=Path(outdir) / "Figure4J_HG002_V2_profile.csv",
        data_root=args.data_root,
    )

    print(f"[OK] Figure 4I, J written under -> {outdir}")


if __name__ == "__main__":
    main()