from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt
from scipy.signal import find_peaks

from figure1_common import ensure_dir, kmer_repetitiveness_signal, load_manifest, read_single_fasta, sanitize_id_for_filename


def find_optimal_scale(abs_coeffs_matrix: np.ndarray, scales_array: np.ndarray, peak_prom_frac: float = 0.05, width_for_scan: int = 3) -> float:
    max_avg_peak = -1.0
    best_scale = None
    for i, sc in enumerate(scales_array):
        row = abs_coeffs_matrix[i, :]
        rmin, rmax = float(np.min(row)), float(np.max(row))
        rrange = rmax - rmin
        if rrange <= 0:
            continue
        prominence = rrange * peak_prom_frac
        peaks, _ = find_peaks(row, prominence=prominence, width=width_for_scan)
        if peaks.size > 0:
            avgp = float(np.mean(row[peaks]))
            if avgp > max_avg_peak:
                max_avg_peak = avgp
                best_scale = float(sc)
    if best_scale is None and len(scales_array) > 0:
        best_scale = float(scales_array[len(scales_array) // 2])
    return best_scale


def compute_cwt_summary(seq_id: str, seq: str, k: int = 20, scale_min: int = 200, scale_max: int = 400) -> dict:
    signal = kmer_repetitiveness_signal(seq, k=k)
    if signal.size == 0:
        raise ValueError(f"Sequence too short for k={k}: {seq_id}")
    scales = np.arange(scale_min, scale_max + 1, dtype=int)
    coeffs, _ = pywt.cwt(signal, scales, "morl")
    abs_coeffs = np.abs(coeffs)
    opt_scale = find_optimal_scale(abs_coeffs, scales)
    scale_idx = int(np.argmin(np.abs(scales - opt_scale)))
    primary = abs_coeffs[scale_idx, :]
    return {
        "seq_id": seq_id,
        "sequence_length": len(seq),
        "signal_length": int(signal.size),
        "scale_min": scale_min,
        "scale_max": scale_max,
        "dominant_scale": float(scales[scale_idx]),
        "primary_mean": float(np.mean(primary)),
        "primary_max": float(np.max(primary)),
        "signal": signal,
        "scales": scales,
        "abs_coeffs": abs_coeffs,
    }


def plot_scalogram(summary: dict, out_path: str | Path) -> None:
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    signal = summary["signal"]
    scales = summary["scales"]
    abs_coeffs = summary["abs_coeffs"]
    dominant_scale = summary["dominant_scale"]

    plt.figure(figsize=(10, 5))
    plt.imshow(
        abs_coeffs,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
        extent=[0, signal.shape[0], scales[0], scales[-1]],
    )
    plt.axhline(dominant_scale, linestyle="--", linewidth=1.5)
    plt.title(f"Scalogram (CWT) — {summary['seq_id']} | dominant scale={dominant_scale:.1f}")
    plt.xlabel("Position (k-mer index)")
    plt.ylabel("Wavelet scale")
    plt.colorbar(label="|Coefficient|")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 1 CWT scalograms")
    parser.add_argument("--manifest", required=True, help="Path to figure1_manifest.tsv")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing extracted Figure 1 FASTA files")
    parser.add_argument("--outdir", required=True, help="Output directory for CWT plots and summary table")
    parser.add_argument("--k", type=int, default=20)
    parser.add_argument("--scale-min", type=int, default=200)
    parser.add_argument("--scale-max", type=int, default=400)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    manifest = manifest[manifest["use_for_cwt"].astype(int) == 1].copy()
    outdir = ensure_dir(args.outdir)

    rows = []
    for row in manifest.itertuples(index=False):
        fasta_path = Path(args.fasta_dir) / f"{row.panel}_{row.species}_{sanitize_id_for_filename(row.record_id)}.fa"
        seq_id, seq = read_single_fasta(fasta_path)
        summary = compute_cwt_summary(seq_id=seq_id, seq=seq, k=args.k, scale_min=args.scale_min, scale_max=args.scale_max)
        png_path = outdir / f"{row.panel}_{row.species}_{sanitize_id_for_filename(row.record_id)}_scalogram.png"
        plot_scalogram(summary, png_path)
        rows.append({
            "panel": row.panel,
            "species": row.species,
            "record_id": row.record_id,
            "seq_id": seq_id,
            "sequence_length": summary["sequence_length"],
            "signal_length": summary["signal_length"],
            "dominant_scale": summary["dominant_scale"],
            "primary_mean": summary["primary_mean"],
            "primary_max": summary["primary_max"],
            "plot_path": str(png_path),
        })
        print(f"[OK] {row.panel}: saved {png_path}")

    df = pd.DataFrame(rows)
    summary_csv = outdir / "figure1_cwt_summary.csv"
    df.to_csv(summary_csv, index=False)
    print(f"[OK] summary -> {summary_csv}")


if __name__ == "__main__":
    main()
