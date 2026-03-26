from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt

from figure2_common import (
    build_log_scales,
    choose_downsample_factor,
    downsample_positions_mean,
    downsample_signal_mean,
    ensure_dir,
    find_optimal_scale,
    get_centromere_bed,
    get_chromosome_sequence,
    get_repeatmasker_out,
    parse_repeatmasker_tar1,
    read_centromere_interval_from_bed,
    sampled_kmer_repetitiveness_signal,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 2D centromere-focused chr13 CWT")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--chrom", default="chr13", help="Chromosome for centromere-focused analysis")
    parser.add_argument("--window-bp-cap", type=int, default=12_000_000, help="Maximum extracted centromere window length")
    parser.add_argument("--kmer-k", type=int, default=20)
    parser.add_argument("--kmer-stride-bp", type=int, default=500)
    parser.add_argument("--target-signal-len", type=int, default=12_000)
    parser.add_argument("--scale-min", type=int, default=100)
    parser.add_argument("--scale-max", type=int, default=400)
    parser.add_argument("--n-scales", type=int, default=64)
    parser.add_argument("--wavelet", default="morl")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    chrom = str(args.chrom)
    _, full_seq = get_chromosome_sequence("human", chrom, data_root=args.data_root)
    chrom_len = len(full_seq)

    centromere_bed = get_centromere_bed(args.data_root)
    cent = read_centromere_interval_from_bed(centromere_bed, chrom)
    if cent is None:
        raise ValueError(f"No centromere interval found for {chrom} in {centromere_bed}")

    cent_start, cent_end = cent
    cent_len = cent_end - cent_start
    if cent_len <= 0:
        raise ValueError(f"Invalid centromere interval for {chrom}: {cent}")

    if cent_len > int(args.window_bp_cap):
        center = (cent_start + cent_end) // 2
        half = int(args.window_bp_cap) // 2
        win_start = max(0, center - half)
        win_end = min(chrom_len, center + half)
    else:
        win_start, win_end = cent_start, cent_end

    seq = full_seq[win_start:win_end]

    rm_out = get_repeatmasker_out(args.data_root)
    tar1_df = parse_repeatmasker_tar1(rm_out, chrom_filter=chrom)
    if tar1_df.empty:
        tar1_df = pd.DataFrame(columns=["chrom", "start", "end", "rep_name", "rep_class_family"])

    tar1_df = tar1_df[(tar1_df["end"] >= win_start) & (tar1_df["start"] <= win_end)].copy()
    if not tar1_df.empty:
        tar1_df["local_start"] = (tar1_df["start"] - win_start).clip(lower=0)
        tar1_df["local_end"] = (tar1_df["end"] - win_start).clip(upper=len(seq))
    else:
        tar1_df["local_start"] = []
        tar1_df["local_end"] = []

    raw_positions_bp, raw_signal = sampled_kmer_repetitiveness_signal(
        seq=seq,
        k=int(args.kmer_k),
        stride_bp=int(args.kmer_stride_bp),
    )
    if raw_signal.size == 0:
        raise ValueError("Empty signal for centromere window")

    ds_factor = choose_downsample_factor(len(raw_signal), target_len=int(args.target_signal_len))
    positions_bp = downsample_positions_mean(raw_positions_bp, ds_factor)
    signal = downsample_signal_mean(raw_signal, ds_factor)

    scales = build_log_scales(
        scale_min=int(args.scale_min),
        scale_max=int(args.scale_max),
        n_scales=int(args.n_scales),
    )
    coeffs, _ = pywt.cwt(signal, scales, args.wavelet)
    abs_coeffs = np.abs(coeffs)

    opt_scale = find_optimal_scale(abs_coeffs, scales, peak_prom_frac=0.05, width_for_scan=3)
    scale_idx = int(np.argmin(np.abs(scales - opt_scale)))
    primary = abs_coeffs[scale_idx, :]
    effective_stride_bp = int(args.kmer_stride_bp) * int(ds_factor)
    scale_used_bp = float(scales[scale_idx] * effective_stride_bp)

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 7),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
        constrained_layout=True,
    )

    x_min = 0.0
    x_max = float(win_end - win_start)

    im = axes[0].imshow(
        abs_coeffs,
        extent=[x_min, x_max, float(scales[0] * effective_stride_bp), float(scales[-1] * effective_stride_bp)],
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    for row in tar1_df.itertuples(index=False):
        axes[0].axvspan(float(row.local_start), float(row.local_end), color="red", alpha=0.30)

    axes[0].set_title(f"Figure 2D: human {chrom} centromeric region")
    axes[0].set_ylabel("Approx scale (bp)")
    plt.colorbar(im, ax=axes[0], fraction=0.03, pad=0.02)

    axes[1].plot(positions_bp, primary, linewidth=1.2, color="black")
    for row in tar1_df.itertuples(index=False):
        axes[1].axvspan(float(row.local_start), float(row.local_end), color="red", alpha=0.30)

    ymin, ymax = axes[1].get_ylim()
    y_tar1 = ymin + 0.03 * (ymax - ymin)
    for row in tar1_df.itertuples(index=False):
        axes[1].plot(
            [float(row.local_start), float(row.local_end)],
            [y_tar1, y_tar1],
            color="red",
            linewidth=5,
            solid_capstyle="butt",
        )

    axes[1].set_title(f"Primary scale profile (scale≈{int(scale_used_bp):,} bp)")
    axes[1].set_xlabel("Position along centromere window (bp)")
    axes[1].set_ylabel("Amplitude")

    out_png = Path(outdir) / "Figure2D_human_chr13_centromere.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    summary = pd.DataFrame([{
        "chrom": chrom,
        "window_start_bp": int(win_start),
        "window_end_bp": int(win_end),
        "window_len_bp": int(len(seq)),
        "n_tar1_fragments": int(len(tar1_df)),
        "raw_signal_len": int(len(raw_signal)),
        "downsample_factor": int(ds_factor),
        "signal_len_after_downsampling": int(len(signal)),
        "kmer_stride_bp": int(args.kmer_stride_bp),
        "effective_stride_bp": int(effective_stride_bp),
        "optimal_scale_ds": float(scales[scale_idx]),
        "optimal_scale_bp": float(scale_used_bp),
        "plot_path": str(out_png),
    }])
    summary_csv = Path(outdir) / "Figure2D_summary.csv"
    summary.to_csv(summary_csv, index=False)

    tar1_csv = Path(outdir) / "Figure2D_tar1_fragments.csv"
    tar1_df.to_csv(tar1_csv, index=False)

    print(f"[OK] plot -> {out_png}")
    print(f"[OK] summary -> {summary_csv}")
    print(f"[OK] TAR1 table -> {tar1_csv}")


if __name__ == "__main__":
    main()