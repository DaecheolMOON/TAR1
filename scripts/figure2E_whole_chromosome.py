from __future__ import annotations

import argparse
from pathlib import Path
from collections import Counter

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
    get_chromosome_sequence,
    load_manifest,
    whole_chromosome_tar1_intervals,
)


def sampled_kmer_repetitiveness_signal(seq: str, k: int = 20, stride_bp: int = 1):
    seq = seq.upper()
    if len(seq) < k:
        return np.array([], dtype=int), np.array([], dtype=float)

    starts = np.arange(0, len(seq) - k + 1, max(1, int(stride_bp)), dtype=int)
    kmers = [seq[i:i + k] for i in starts]
    freq = Counter(kmers)
    signal = np.array([freq[km] for km in kmers], dtype=float)
    return starts, signal


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 2E whole-chromosome CWT examples")
    parser.add_argument("--manifest", required=True, help="Path to figure2E_manifest.tsv")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--kmer-k", type=int, default=20)
    parser.add_argument("--kmer-stride-bp", type=int, default=2000)
    parser.add_argument("--target-signal-len", type=int, default=10_000)
    parser.add_argument("--scale-min", type=int, default=100)
    parser.add_argument("--scale-max", type=int, default=400)
    parser.add_argument("--n-scales", type=int, default=64)
    parser.add_argument("--wavelet", default="morl")
    parser.add_argument("--telomere-zone-bp", type=int, default=1_000_000)
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    plots_dir = ensure_dir(Path(outdir) / "plots")

    manifest = load_manifest(args.manifest)
    scales = build_log_scales(
        scale_min=int(args.scale_min),
        scale_max=int(args.scale_max),
        n_scales=int(args.n_scales),
    )

    summary_rows = []

    for row in manifest.itertuples(index=False):
        chrom = str(row.chrom)
        species = str(row.species)

        _, seq = get_chromosome_sequence(species, chrom, data_root=args.data_root)
        chrom_len = len(seq)

        raw_positions_bp, raw_signal = sampled_kmer_repetitiveness_signal(
            seq=seq,
            k=int(args.kmer_k),
            stride_bp=int(args.kmer_stride_bp),
        )
        if raw_signal.size == 0:
            raise ValueError(f"Empty signal for {species} {chrom}")

        ds_factor = choose_downsample_factor(len(raw_signal), target_len=int(args.target_signal_len))
        positions_bp = downsample_positions_mean(raw_positions_bp, ds_factor)
        signal = downsample_signal_mean(raw_signal, ds_factor)

        if signal.size < 16:
            raise ValueError(f"Signal too short after downsampling for {species} {chrom}")

        coeffs, _ = pywt.cwt(signal, scales, args.wavelet)
        abs_coeffs = np.abs(coeffs)

        opt_scale = find_optimal_scale(abs_coeffs, scales, peak_prom_frac=0.05, width_for_scan=3)
        scale_idx = int(np.argmin(np.abs(scales - opt_scale)))
        primary = abs_coeffs[scale_idx, :]
        effective_stride_bp = int(args.kmer_stride_bp) * int(ds_factor)
        scale_used_bp = float(scales[scale_idx] * effective_stride_bp)

        tar1_intervals = whole_chromosome_tar1_intervals(
            species=species,
            chrom=chrom,
            chrom_len=chrom_len,
            data_root=args.data_root,
        )

        fig, axes = plt.subplots(
            2, 1, figsize=(13, 7),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True,
            constrained_layout=True,
        )

        x_min = 0.0
        x_max = float(chrom_len)

        im = axes[0].imshow(
            abs_coeffs,
            extent=[x_min, x_max, float(scales[0] * effective_stride_bp), float(scales[-1] * effective_stride_bp)],
            aspect="auto",
            origin="lower",
            interpolation="nearest",
        )

        tel_zone = min(int(args.telomere_zone_bp), chrom_len // 2)
        axes[0].axvspan(0, tel_zone, color="gray", alpha=0.10)
        axes[0].axvspan(chrom_len - tel_zone, chrom_len, color="gray", alpha=0.10)

        for s, e, arm in tar1_intervals:
            axes[0].axvspan(float(s), float(e), color="orange", alpha=0.25)

        axes[0].set_title(f"{row.panel}: {species} {chrom}")
        axes[0].set_ylabel("Approx scale (bp)")
        plt.colorbar(im, ax=axes[0], fraction=0.03, pad=0.02)

        axes[1].plot(positions_bp, primary, linewidth=1.2, color="black")
        axes[1].axvspan(0, tel_zone, color="gray", alpha=0.10)
        axes[1].axvspan(chrom_len - tel_zone, chrom_len, color="gray", alpha=0.10)

        for s, e, arm in tar1_intervals:
            axes[1].axvspan(float(s), float(e), color="orange", alpha=0.25)

        ymin, ymax = axes[1].get_ylim()
        y_tar1 = ymin + 0.03 * (ymax - ymin)
        for s, e, arm in tar1_intervals:
            axes[1].plot(
                [float(s), float(e)],
                [y_tar1, y_tar1],
                color="orange",
                linewidth=5,
                solid_capstyle="butt",
            )

        axes[1].set_title(f"Primary scale profile (scale≈{int(scale_used_bp):,} bp)")
        axes[1].set_xlabel("Chromosome coordinate (bp)")
        axes[1].set_ylabel("Amplitude")

        out_png = plots_dir / f"{row.panel}_{row.species}_{chrom}_whole_chromosome.png"
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

        summary_rows.append({
            "panel": row.panel,
            "species": species,
            "chrom": chrom,
            "chrom_len_bp": int(chrom_len),
            "n_tar1_intervals": int(len(tar1_intervals)),
            "raw_signal_len": int(len(raw_signal)),
            "downsample_factor": int(ds_factor),
            "signal_len_after_downsampling": int(len(signal)),
            "kmer_stride_bp": int(args.kmer_stride_bp),
            "effective_stride_bp": int(effective_stride_bp),
            "optimal_scale_ds": float(scales[scale_idx]),
            "optimal_scale_bp": float(scale_used_bp),
            "plot_path": str(out_png),
        })

        print(f"[OK] {row.panel}: saved {out_png}")

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = Path(outdir) / "Figure2E_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"[OK] summary -> {summary_csv}")


if __name__ == "__main__":
    main()