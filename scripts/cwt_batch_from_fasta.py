#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch CWT analysis on chromosome-arm FASTA sequences.

For each record (e.g. chr1p, chr1q):
  1) Build a k-mer repetitiveness signal (k-mer frequency per position).
  2) Compute a CWT scalogram (Morl wavelet).
  3) Select an "optimal" scale by average peak strength.
  4) Detect strong (peaks) and weak (troughs) runs in the primary scale.
  5) For each run, compute AT/CG ratio, sample sequence, and coordinates.
  6) Save a scalogram plot and append run-level information to a summary CSV.
"""

import os
import argparse
from collections import Counter

import numpy as np
import pandas as pd
import pywt
from scipy.signal import find_peaks
from Bio import SeqIO
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------- k-mer repetitiveness signal ----------
def build_kmer_signal(seq, k=20):
    """
    Build a k-mer repetitiveness signal:
    each position reflects how many times the corresponding k-mer appears in the sequence.
    """
    seq = seq.upper()
    if len(seq) < k:
        return np.array([], dtype=float)
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    freq = Counter(kmers)
    signal = np.array([freq[km] for km in kmers], dtype=float)
    return signal


# ---------- contiguous index grouping ----------
def group_contiguous_indices(indices, min_gap=10):
    """
    Group sorted indices into contiguous runs, merging gaps <= min_gap.
    Returns a list of (start_idx, end_idx) pairs.
    """
    indices = np.array(indices, dtype=int)
    if indices.size == 0:
        return []
    if indices.size == 1:
        return [(int(indices[0]), int(indices[0]))]

    sorted_indices = np.sort(indices)
    groups = []
    start = int(sorted_indices[0])
    end = int(sorted_indices[0])
    for x in sorted_indices[1:]:
        x = int(x)
        if x - end <= min_gap:
            end = x
        else:
            groups.append((start, end))
            start = x
            end = x
    groups.append((start, end))
    return groups


# ---------- AT/CG ratio in a local window ----------
def atcg_ratio_calc(sequence, indices, window_size=25):
    """
    Compute AT/CG ratio in a local window around each index.

    Returns a numpy array of ratios (float), possibly containing np.nan or np.inf.
    """
    seq = sequence.upper()
    L = len(seq)
    ratios = []
    for idx in indices:
        start = max(0, idx - window_size // 2)
        end = min(L, idx + window_size // 2 + (window_size % 2))
        window_seq = seq[start:end]
        if not window_seq:
            ratios.append(np.nan)
            continue
        at_count = window_seq.count("A") + window_seq.count("T")
        cg_count = window_seq.count("C") + window_seq.count("G")
        if cg_count == 0:
            ratios.append(np.inf if at_count > 0 else np.nan)
        else:
            ratios.append(at_count / cg_count)
    return np.array(ratios, dtype=float)


# ---------- find optimal scale by average peak strength ----------
def find_optimal_scale(abs_coeffs_matrix, scales_array, peak_prom_frac=0.05, width_for_scan=3):
    """
    Select the wavelet scale with the highest average peak strength.

    For each scale:
      - compute dynamic range (min, max)
      - set peak prominence as 'peak_prom_frac * range'
      - detect peaks, then compute the mean value at detected peaks
    The scale with the largest mean peak value is returned.
    """
    max_avg_peak = -1.0
    best_scale = None
    for i, sc in enumerate(scales_array):
        s = abs_coeffs_matrix[i, :]
        s_min, s_max = np.min(s), np.max(s)
        s_range = s_max - s_min
        if s_range <= 0:
            continue
        prominence = s_range * peak_prom_frac
        peaks, _ = find_peaks(s, prominence=prominence, width=width_for_scan)
        if peaks.size > 0:
            avg_peak_strength = float(np.mean(s[peaks]))
            if avg_peak_strength > max_avg_peak:
                max_avg_peak = avg_peak_strength
                best_scale = float(sc)
    if best_scale is None and len(scales_array) > 0:
        best_scale = float(scales_array[len(scales_array) // 2])
    return best_scale


# ---------- analyze one sequence, save plot, return runs and sample info ----------
def analyze_one_sequence(
    seq_id,
    seq_str,
    k=20,
    scales=(200, 401),
    peak_keep_frac=0.5,
    prominence_frac=0.10,
    width_threshold=5,
    out_dir=".",
    save_plot=True,
):
    """
    Run k-mer-based CWT analysis on a single sequence.

    Returns:
        optimal_scale (float or None),
        result_events (list of dict),
        fig_path (str or None)
    """
    result_events = []

    # build signal
    signal = build_kmer_signal(seq_str, k=k)
    if signal.size == 0:
        return None, result_events, None

    # CWT
    sc_arr = np.arange(scales[0], scales[1], dtype=int)
    coeffs, _ = pywt.cwt(signal, sc_arr, "morl")
    abs_coeffs = np.abs(coeffs)

    # optimal scale
    optimal_scale = find_optimal_scale(abs_coeffs, sc_arr)
    scale_idx = int(np.argmin(np.abs(sc_arr - optimal_scale)))
    primary = abs_coeffs[scale_idx, :]
    smin, smax = float(np.min(primary)), float(np.max(primary))
    srange = smax - smin
    if srange <= 0:
        return float(sc_arr[scale_idx]), result_events, None

    # thresholds
    prom = srange * prominence_frac
    if prom == 0:
        prom = float(np.mean(primary)) * 0.01
    peaks_all, _ = find_peaks(primary, prominence=prom, width=width_threshold)
    troughs_all, _ = find_peaks(-primary, prominence=prom, width=width_threshold)

    # keep strongest X% peaks
    if peaks_all.size > 0:
        vals = primary[peaks_all]
        n_keep = int(round(peak_keep_frac * len(peaks_all)))
        n_keep = max(n_keep, 0)
        if n_keep > 0:
            top_idx = np.argsort(vals)[::-1][:n_keep]
            peaks = np.sort(peaks_all[top_idx])
        else:
            peaks = np.array([], dtype=int)
    else:
        peaks = np.array([], dtype=int)

    troughs = np.sort(troughs_all) if troughs_all.size > 0 else np.array([], dtype=int)

    # group to runs (k-mer indices)
    gap_for_grouping = int(k * 2.5)
    peak_runs = group_contiguous_indices(peaks, min_gap=gap_for_grouping)
    trough_runs = group_contiguous_indices(troughs, min_gap=gap_for_grouping)

    # prepare plotting
    fig_path = None
    if save_plot:
        plt.figure(figsize=(12, 7))
        plt.imshow(
            np.abs(abs_coeffs),
            aspect="auto",
            interpolation="nearest",
            origin="lower",
            extent=[0, signal.shape[0], sc_arr[0], sc_arr[-1]],
        )
        plt.xlabel("Sequence Position (k-mer index)")
        plt.ylabel("Wavelet Scale")
        plt.title(f"Scalogram (CWT) of k-mer Repetitiveness — {seq_id}")
        cbar = plt.colorbar()
        cbar.set_label("|Coefficient| (Signal Strength)")
        fig_path = os.path.join(out_dir, f"{seq_id}_scalogram.png")
        plt.tight_layout()
        plt.savefig(fig_path, dpi=200)
        plt.close()

    # analyze both strong and weak runs and append to result_events
    window_size_for_atcg = k * 2

    def process_runs(runs, signal_type):
        for (start_idx, end_idx) in runs:
            # run k-mer indices
            kmer_indices = np.arange(start_idx, end_idx + 1, dtype=int)
            # compute mean AT/CG over run
            atcg_vals = atcg_ratio_calc(
                seq_str, kmer_indices, window_size=window_size_for_atcg
            )
            atcg_vals_clean = atcg_vals[
                ~np.isnan(atcg_vals) & ~np.isinf(atcg_vals)
            ]
            mean_atcg = (
                float(np.mean(atcg_vals_clean))
                if atcg_vals_clean.size > 0
                else np.nan
            )

            # representative index (use start of run)
            rep_kmer_idx = int(start_idx)

            # sample sequence window (base coordinates)
            sample_seq_start = max(0, rep_kmer_idx - window_size_for_atcg // 2)
            sample_seq_end = min(
                len(seq_str),
                rep_kmer_idx
                + window_size_for_atcg // 2
                + (window_size_for_atcg % 2),
            )
            sample_seq = seq_str[sample_seq_start:sample_seq_end]

            # convert k-mer indices to base coordinates (0-based)
            start_base = int(start_idx)
            end_base = int(end_idx + k - 1)  # inclusive

            result_events.append(
                {
                    "record_id": seq_id,
                    "signal_type": signal_type,
                    "start_kmer_idx": int(start_idx),
                    "end_kmer_idx": int(end_idx),
                    "start_base": start_base,
                    "end_base": end_base,
                    "run_kmer_len": int(end_idx - start_idx + 1),
                    "scale_used": float(sc_arr[scale_idx]),
                    "mean_atcg": mean_atcg,
                    "sample_seq_start": int(sample_seq_start),
                    "sample_seq_end": int(sample_seq_end - 1),  # inclusive
                    "sample_seq_len": len(sample_seq),
                    "sample_seq": sample_seq,
                }
            )

    process_runs(peak_runs, "strong")
    process_runs(trough_runs, "weak")

    return float(sc_arr[scale_idx]), result_events, fig_path


# ---------- main batch ----------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch CWT analysis on FASTA records and save scalogram plots and a run-level CSV "
            "containing sample sequences, AT/CG ratios, and coordinates."
        )
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file with records (e.g. chr1p, chr1q).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory to save plots and the summary CSV.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=20,
        help="k-mer length (default: 20)",
    )
    parser.add_argument(
        "--scale_min",
        type=int,
        default=200,
        help="Minimum wavelet scale (inclusive, default: 200)",
    )
    parser.add_argument(
        "--scale_max",
        type=int,
        default=401,
        help="Maximum wavelet scale (exclusive, default: 401)",
    )
    parser.add_argument(
        "--peak_keep_frac",
        type=float,
        default=0.5,
        help="Fraction of strongest peaks to keep (default: 0.5)",
    )
    parser.add_argument(
        "--prominence_frac",
        type=float,
        default=0.10,
        help="Prominence fraction of the dynamic range (default: 0.10)",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=5,
        help="Minimum peak width for find_peaks (default: 5)",
    )
    parser.add_argument(
        "--log_each",
        action="store_true",
        help="If set, write a small text log per record.",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    all_rows = []
    n_records = 0

    for rec in SeqIO.parse(args.fasta, "fasta"):
        n_records += 1
        seq_id = rec.id  # expected like 'chr1p'/'chr1q'
        seq_str = str(rec.seq).upper().replace("\n", "").replace(" ", "")

        scale_used, events, fig_path = analyze_one_sequence(
            seq_id=seq_id,
            seq_str=seq_str,
            k=args.k,
            scales=(args.scale_min, args.scale_max),
            peak_keep_frac=args.peak_keep_frac,
            prominence_frac=args.prominence_frac,
            width_threshold=args.width,
            out_dir=args.outdir,
            save_plot=True,
        )

        # extend master list
        if events:
            all_rows.extend(events)

        # optional per-record log
        if args.log_each:
            log_path = os.path.join(args.outdir, f"{seq_id}_log.txt")
            with open(log_path, "w") as f:
                f.write(f"[{seq_id}] scale_used: {scale_used}\n")
                f.write(f"events_count: {len(events)}\n")
                f.write(f"figure: {fig_path if fig_path else 'NA'}\n")

    # write summary CSV
    csv_path = os.path.join(args.outdir, "summary_runs.csv")
    if all_rows:
        df = pd.DataFrame(all_rows)
        cols = [
            "record_id",
            "signal_type",
            "start_kmer_idx",
            "end_kmer_idx",
            "start_base",
            "end_base",
            "run_kmer_len",
            "scale_used",
            "mean_atcg",
            "sample_seq_start",
            "sample_seq_end",
            "sample_seq_len",
            "sample_seq",
        ]
        df = df[cols]
        df.to_csv(csv_path, index=False)
    else:
        pd.DataFrame(
            columns=[
                "record_id",
                "signal_type",
                "start_kmer_idx",
                "end_kmer_idx",
                "start_base",
                "end_base",
                "run_kmer_len",
                "scale_used",
                "mean_atcg",
                "sample_seq_start",
                "sample_seq_end",
                "sample_seq_len",
                "sample_seq",
            ]
        ).to_csv(csv_path, index=False)

    print(f"\nDone. Records processed: {n_records}")
    print(f"Plots & CSV saved to: {args.outdir}")
    print(f"Summary CSV: {csv_path}")


if __name__ == "__main__":
    main()
