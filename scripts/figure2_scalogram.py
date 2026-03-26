from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt

from figure2_common import (
    ensure_dir,
    load_manifest,
    read_single_fasta,
    sampled_kmer_repetitiveness_signal,
    choose_downsample_factor,
    downsample_positions_mean,
    downsample_signal_mean,
    sanitize_id_for_filename,
    tar1_intervals_within_window,
)


def build_fasta_name(panel: str, species: str, arm: str, window_bp: int) -> str:
    return f"{panel}_{species}_{sanitize_id_for_filename(arm)}_{int(window_bp)}bp.fa"


def analyze_signal_cwt(
    seq: str,
    window_bp: int,
    kmer_k: int = 20,
    stride_for_large: int = 200,
    target_signal_len_large: int = 40000,
):
    """
    Policy for Figure 2A:
      - 30 kb / 500 kb: no downsampling, scales 200-400
      - 10 Mb: sampled + adaptive downsampling, still scales 200-400
    """
    window_bp = int(window_bp)

    if window_bp in (30_000, 500_000):
        raw_positions_bp, raw_signal = sampled_kmer_repetitiveness_signal(
            seq=seq,
            k=int(kmer_k),
            stride_bp=1,
        )
        ds_factor = 1
        positions_bp = raw_positions_bp.astype(float)
        signal = raw_signal.astype(float)
        stride_bp = 1

    elif window_bp == 10_000_000:
        raw_positions_bp, raw_signal = sampled_kmer_repetitiveness_signal(
            seq=seq,
            k=int(kmer_k),
            stride_bp=int(stride_for_large),
        )
        ds_factor = choose_downsample_factor(
            len(raw_signal),
            target_len=int(target_signal_len_large),
        )
        positions_bp = downsample_positions_mean(raw_positions_bp, ds_factor)
        signal = downsample_signal_mean(raw_signal, ds_factor)
        stride_bp = int(stride_for_large)

    else:
        # fallback
        raw_positions_bp, raw_signal = sampled_kmer_repetitiveness_signal(
            seq=seq,
            k=int(kmer_k),
            stride_bp=1,
        )
        ds_factor = 1
        positions_bp = raw_positions_bp.astype(float)
        signal = raw_signal.astype(float)
        stride_bp = 1

    if signal.size == 0:
        raise ValueError("Signal is empty after preprocessing.")

    scales = np.arange(200, 401, dtype=int)
    coeffs, _ = pywt.cwt(signal, scales, "morl")
    abs_coeffs = np.abs(coeffs)

    mean_power = abs_coeffs.mean(axis=1)
    best_idx = int(np.argmax(mean_power))
    best_scale = int(scales[best_idx])

    return {
        "raw_signal_len": int(len(raw_signal)),
        "positions_bp": positions_bp,
        "signal": signal,
        "ds_factor": int(ds_factor),
        "stride_bp": int(stride_bp),
        "effective_stride_bp": int(stride_bp * ds_factor),
        "scales": scales,
        "coeffs_abs": abs_coeffs,
        "best_scale": best_scale,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 2A CWT scalograms with TAR1 overlay")
    parser.add_argument("--manifest", required=True, help="Path to figure2_manifest.tsv")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--fasta-dir", required=True, help="Prepared window FASTA directory")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--kmer-k", type=int, default=20)
    parser.add_argument("--large-window-stride-bp", type=int, default=500)
    parser.add_argument("--large-window-target-signal-len", type=int, default=25000)
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    plot_dir = ensure_dir(Path(outdir) / "plots")
    summary_rows = []

    manifest = load_manifest(args.manifest)

    for row in manifest.itertuples(index=False):
        fasta_name = build_fasta_name(row.panel, row.species, row.arm, int(row.window_bp))
        fasta_path = Path(args.fasta_dir) / fasta_name
        if not fasta_path.exists():
            raise FileNotFoundError(f"Prepared FASTA not found: {fasta_path}")

        header, seq = read_single_fasta(fasta_path)

        res = analyze_signal_cwt(
            seq=seq,
            window_bp=int(row.window_bp),
            kmer_k=int(args.kmer_k),
            stride_for_large=int(args.large_window_stride_bp),
            target_signal_len_large=int(args.large_window_target_signal_len),
        )

        positions_bp = res["positions_bp"]
        signal = res["signal"]
        scales = res["scales"]
        abs_coeffs = res["coeffs_abs"]
        best_scale = res["best_scale"]

        tar1_intervals = tar1_intervals_within_window(
            species=row.species,
            arm=row.arm,
            window_bp=int(row.window_bp),
            data_root=args.data_root,
        )

        fig, ax = plt.subplots(figsize=(10, 4))

        # IMPORTANT:
        # small scale at bottom, large scale at top
        # so use extent=[x0, x1, scales[0], scales[-1]] with origin="lower"
        x_min = 0.0
        x_max = float(positions_bp[-1]) if len(positions_bp) > 0 else float(len(seq))

        im = ax.imshow(
            abs_coeffs,
            extent=[x_min, x_max, float(scales[0]), float(scales[-1])],
            aspect="auto",
            origin="lower",
        )

        for s, e in tar1_intervals:
            ax.axvspan(float(s), float(e), alpha=0.25)

        ax.set_title(f"{row.panel}: {row.species} {row.arm} | {int(row.window_bp):,} bp")
        ax.set_xlabel("Distance from telomere (bp)")
        ax.set_ylabel("Scale")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

        out_png = plot_dir / f"{row.panel}_{row.species}_{sanitize_id_for_filename(row.arm)}_{int(row.window_bp)}bp_scalogram.png"
        plt.tight_layout()
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

        summary_rows.append({
            "panel": row.panel,
            "species": row.species,
            "arm": row.arm,
            "window_bp": int(row.window_bp),
            "raw_signal_len": int(res["raw_signal_len"]),
            "signal_len_after_processing": int(len(signal)),
            "downsample_factor": int(res["ds_factor"]),
            "kmer_stride_bp": int(res["stride_bp"]),
            "effective_stride_bp": int(res["effective_stride_bp"]),
            "n_tar1_intervals": int(len(tar1_intervals)),
            "best_scale": int(best_scale),
            "plot_path": str(out_png),
        })

        print(
            f"[OK] {row.panel}: saved {out_png} "
            f"(window_bp={int(row.window_bp)}, raw_signal_len={int(res['raw_signal_len'])}, "
            f"processed_signal_len={int(len(signal))}, stride_bp={int(res['stride_bp'])}, "
            f"downsample_factor={int(res['ds_factor'])}, best_scale={int(best_scale)})"
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = Path(outdir) / "figure2A_scalogram_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"[OK] summary -> {summary_csv}")


if __name__ == "__main__":
    main()