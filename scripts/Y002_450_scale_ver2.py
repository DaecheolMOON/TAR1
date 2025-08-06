#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAR1 repeat SVD analysis (read-level, Mode 2 only, adaptive scale ≤ 800 bp)

1) k-mer count signal  → CWT
2) raw energy matrix   → SVD
3) Mode 2 scale-loading plot + two leading peak positions
"""

import numpy as np
import pywt
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from collections import Counter
import argparse
from scipy.signal import find_peaks

# ───── global parameters ──────────────────────────────────────────
KMER_K          = 20
MODES           = [2]           # Mode 2 only
INIT_MIN_SCALE  = 200
INIT_MAX_SCALE  = 430
MAX_ALLOWED     = 800           # hard ceiling
STEP_EXTEND     = 20            # step for extension
PEAK_DIST       = 50            # min distance between peaks
PEAK_HEIGHT     = 0.05          # relative threshold
# ──────────────────────────────────────────────────────────────────


def compute_cwt_matrix(seq: str, scales: np.ndarray) -> np.ndarray:
    """Return |CWT| coefficients summed over positions for given scales."""
    kmers = [seq[i:i + KMER_K] for i in range(len(seq) - KMER_K + 1)]
    counts = Counter(kmers)
    signal = np.array([counts[k] for k in kmers])
    coeffs, _ = pywt.cwt(signal, scales, 'morl', method='fft')
    return np.abs(coeffs)               # shape: (n_scales, len(signal))


def ensure_two_peaks(vector: np.ndarray,
                     scales: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (peak_indices, peak_scales) ensuring ≥2 peaks if possible."""
    peaks, _ = find_peaks(np.abs(vector),
                          distance=PEAK_DIST,
                          height=PEAK_HEIGHT)
    return peaks, scales[peaks]


def analyze_single_fasta(fasta_path: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n==== {os.path.basename(fasta_path)} → {output_dir}")

    # 1) load FASTA
    reads = {rec.id: str(rec.seq).upper()
             for rec in SeqIO.parse(fasta_path, "fasta")}
    if not reads:
        print("[WARN] empty FASTA, skip")
        return

    # 2) adaptive CWT–SVD loop
    max_scale = INIT_MAX_SCALE
    peak_scales = None
    final_scales = None
    U = S = Vt = None

    while max_scale <= MAX_ALLOWED:
        scales = np.arange(INIT_MIN_SCALE, max_scale + 1)
        cwt_sum = []
        processed_ids = []

        for rid, seq in reads.items():
            if len(seq) <= max_scale + KMER_K:
                continue
            mat = compute_cwt_matrix(seq, scales)
            cwt_sum.append(mat.sum(axis=1))
            processed_ids.append(rid)

        if not processed_ids:
            print("[ERROR] no reads long enough, abort")
            return

        M = np.vstack(cwt_sum)
        try:
            U, S, Vt = np.linalg.svd(M, full_matrices=False)
        except np.linalg.LinAlgError as e:
            print(f"[ERROR] SVD failed: {e}")
            return

        vec = Vt[MODES[0] - 1]           # Mode 2 vector
        peaks, p_scales = ensure_two_peaks(vec, scales)

        if len(peaks) >= 2:
            peak_scales = p_scales[:2]   # first two peaks
            final_scales = scales
            break

        max_scale += STEP_EXTEND         # extend range and retry

    if peak_scales is None:
        print("[WARN] <2 peaks even at max range, analysis continues with "
              f"{len(peaks)} peak(s)")
        peak_scales = p_scales
        final_scales = scales

    # 3) report peak positions
    print(f"[INFO] Mode 2 peaks: {peak_scales.tolist()} bp")

    # 4) save summary file
    summary_path = os.path.join(output_dir, "mode2_peak_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Mode 2 peaks (bp): {', '.join(map(str, peak_scales))}\n")
    print(f"[OUT] summary → {summary_path}")

    # 5) save plot
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(final_scales, vec, '.', lw=1, color='royalblue', label='Mode 2')
    for idx, p in enumerate(peak_scales):
        ax.axvline(p, color='red', ls='--', alpha=0.8,
                   label=f'Peak {idx + 1}: {p} bp' if idx == 0 else None)

    ax.set_title(f"Mode 2 scale-loading\n{os.path.basename(fasta_path)}")
    ax.set_xlabel("Scale (bp)")
    ax.set_ylabel("Loading")
    ax.set_xlim(final_scales[0], final_scales[-1])
    ax.grid(True, ls=':')
    ax.legend()
    plt.tight_layout()
    plot_path = os.path.join(output_dir, "mode2_scale_loading.png")
    plt.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"[OUT] plot → {plot_path}")


def main() -> None:
    argp = argparse.ArgumentParser(
        description="Batch TAR1 SVD (read-level, Mode 2 only, adaptive scale)")
    argp.add_argument("--base_dir", required=True,
                      help="Directory containing SRRT*/tar1_blocks.fa")
    args = argp.parse_args()

    base_dir = args.base_dir
    fasta_name = "tar1_blocks.fa"
    out_dir_name = "Y002_430"

    count = 0
    for sub in sorted(os.listdir(base_dir)):
        if not sub.startswith("SRRT"):
            continue
        srr_dir = os.path.join(base_dir, sub)
        fasta_path = os.path.join(srr_dir, fasta_name)
        if not os.path.isfile(fasta_path):
            continue
        count += 1
        try:
            analyze_single_fasta(fasta_path,
                                 os.path.join(srr_dir, out_dir_name))
        except Exception as err:
            print(f"[ERROR] {fasta_path}: {err}")

    if count == 0:
        print("No matching FASTA files found.")
    else:
        print(f"Batch complete: {count} FASTA processed.")


if __name__ == "__main__":
    main()
