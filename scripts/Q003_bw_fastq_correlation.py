import numpy as np
import pywt
import pyBigWig
import matplotlib.pyplot as plt
from Bio import SeqIO
import pandas as pd
from scipy.signal import find_peaks
from collections import Counter
from scipy.stats import chi2_contingency, fisher_exact
import re
import os

# ===== Helper Functions =====

def extract_arm_from_label(label_str):
    """
    Extracts a standard chromosome arm label (e.g., 'chr1p', 'chrxq') from a string using regex.
    """
    # Regex to find patterns like chr1p, chr22q, chrxp, etc.
    match = re.search(r'(chr([1-9]|1[0-9]|2[0-2]|X|Y)[pq])', label_str.lower())
    if match:
        return match.group(1)
    return None

def atcg_ratio_calc(sequence_str, indices, window_size=25):
    """Calculates AT/CG ratio for given indices within a sequence using a sliding window."""
    ratios = []
    L = len(sequence_str)
    for idx in indices:
        seq_idx = idx
        start = max(0, seq_idx - window_size // 2)
        end = min(L, seq_idx + window_size // 2 + (window_size % 2))
        window_seq = sequence_str[start:end]

        if not window_seq:
            ratios.append(np.nan)
            continue

        at_count = window_seq.count('A') + window_seq.count('T')
        cg_count = window_seq.count('C') + window_seq.count('G')

        if cg_count == 0:
            ratios.append(np.inf if at_count > 0 else np.nan)
        else:
            ratios.append(at_count / cg_count)
    return np.array(ratios)

def group_contiguous_indices(indices, min_gap=10):
    """Groups contiguous indices (peaks or troughs) into runs."""
    if not indices.size:
        return []
    if indices.size == 1:
        return [(indices[0], indices[0])]
    sorted_indices = np.sort(indices)
    groups = []
    current_start = sorted_indices[0]
    current_end = sorted_indices[0]
    for i in range(1, len(sorted_indices)):
        if sorted_indices[i] - current_end <= min_gap:
            current_end = sorted_indices[i]
        else:
            groups.append((current_start, current_end))
            current_start = sorted_indices[i]
            current_end = sorted_indices[i]
    groups.append((current_start, current_end))
    return groups

def find_optimal_scale_from_cwt(abs_coeffs_matrix, scales_array):
    """Finds the scale with the highest average peak strength from CWT coefficients."""
    max_avg_peak_strength = -1
    optimal_scale_value = -1
    for i, scale_val in enumerate(scales_array):
        current_signal_strength = abs_coeffs_matrix[i, :]
        s_range = np.max(current_signal_strength) - np.min(current_signal_strength)
        if s_range == 0:
            continue
        prominence_for_scan = s_range * 0.05
        peaks, _ = find_peaks(current_signal_strength, prominence=prominence_for_scan, width=3)
        if peaks.size > 0:
            avg_peak_strength = np.mean(current_signal_strength[peaks])
            if avg_peak_strength > max_avg_peak_strength:
                max_avg_peak_strength = avg_peak_strength
                optimal_scale_value = scale_val
    if optimal_scale_value == -1 and len(scales_array) > 0:
        return scales_array[len(scales_array) // 2] # Return mid-scale as default
    return optimal_scale_value

# ===== Configuration =====
# Define chromosome arms to be processed
arms_to_process = [f"chr{i}{p}" for i in range(1, 23) for p in ("p", "q")]
arms_to_process.extend(["chrxp", "chrxq"])

# File Paths (Please verify these paths)
bed_path  = "/Users/daecheolmoon/Downloads/6p_data/Original/mycode/tar_fasta_masked/chm13/tar1_blocks.bed"
bw_path   = "/Users/daecheolmoon/Downloads/6p_data/Original/mycode/github/mycode/fastq_/WT_H3K9me3.bw"
fasta_path= "/Users/daecheolmoon/Downloads/6p_data/Original/mycode/tar_fasta_masked/chm13/tar1_blocks.fa"

# CWT and Analysis Parameters
common_cwt_scales = np.arange(100, 401)
k_mer_size = 20
smooth_win = 20
window_size_for_atcg_stat = k_mer_size * 2

# Output directory for plots
output_plot_dir = "cwt_analysis_plots"
os.makedirs(output_plot_dir, exist_ok=True)


# ===== Data Loading and Preprocessing =====

# Load and parse BED file into a dictionary for easy lookup
print("Loading and parsing BED file...")
bed_blocks = {}
with open(bed_path) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            chrom, start, end, label = parts[0], int(parts[1]), int(parts[2]), parts[3]
            arm_label = extract_arm_from_label(label)
            if arm_label:
                bed_blocks[arm_label] = (chrom, start, end)
print(f"Loaded {len(bed_blocks)} blocks from BED file.")

# Load and parse FASTA file into a dictionary
print("Loading and parsing FASTA file...")
fasta_seqs = {}
for rec in SeqIO.parse(fasta_path, "fasta"):
    arm_label = extract_arm_from_label(rec.description)
    if arm_label:
        if arm_label not in fasta_seqs:
            fasta_seqs[arm_label] = []
        fasta_seqs[arm_label].append(str(rec.seq).upper())

# Join sequences for arms that were split into multiple records
for arm_label, seq_list in fasta_seqs.items():
    fasta_seqs[arm_label] = "".join(seq_list)
print(f"Loaded and consolidated sequences for {len(fasta_seqs)} arms from FASTA file.")

# Open bigWig file
try:
    bw = pyBigWig.open(bw_path)
    if not bw:
        raise RuntimeError(f"Could not open bigWig file at {bw_path}")
except Exception as e:
    print(f"Error opening bigWig file: {e}")
    exit()

# ===== Main Analysis Loop =====
results_list = []
print("\nStarting analysis for each chromosome arm...")
for arm_label in arms_to_process:
    print(f"\n===== Processing arm: {arm_label} =====")
    arm_results = {'Arm': arm_label}

    # Retrieve data from pre-processed dictionaries
    current_bed_block = bed_blocks.get(arm_label)
    fasta_seq_for_arm = fasta_seqs.get(arm_label)

    if not current_bed_block or not fasta_seq_for_arm:
        print(f"Data missing for arm '{arm_label}'. Skipping.")
        results_list.append(arm_results)
        continue

    chrom, bed_start, bed_end = current_bed_block
    print(f"Found BED block: {chrom}:{bed_start}-{bed_end} | FASTA length: {len(fasta_seq_for_arm)}")

    # Initialize variables for this arm
    optimal_scale_fasta = np.nan
    primary_signal_strength_kmer = None

    # --- FASTA k-mer CWT Analysis ---
    if len(fasta_seq_for_arm) < max(common_cwt_scales) + k_mer_size:
        print(f"FASTA sequence for {arm_label} is too short for analysis.")
    else:
        # Build k-mer signal
        kmers = [fasta_seq_for_arm[i:i+k_mer_size] for i in range(len(fasta_seq_for_arm) - k_mer_size + 1)]
        freq = Counter(kmers)
        kmer_signal = np.array([freq[k_m] for k_m in kmers])

        if len(kmer_signal) < max(common_cwt_scales):
            print(f"K-mer signal for {arm_label} is too short for CWT.")
        else:
            # Perform CWT and find optimal scale
            coeffs_kmer, _ = pywt.cwt(kmer_signal, common_cwt_scales, 'morl')
            abs_coeffs_kmer = np.abs(coeffs_kmer)
            optimal_scale_fasta = find_optimal_scale_from_cwt(abs_coeffs_kmer, common_cwt_scales)

            if optimal_scale_fasta is not None:
                arm_results['Optimal Scale (FASTA k-mer)'] = optimal_scale_fasta
                print(f"Optimal FASTA k-mer CWT scale: {optimal_scale_fasta:.1f}")

                # --- Statistical Analysis (based on Representative Sample Sequences) ---
                optimal_scale_idx = np.argmin(np.abs(common_cwt_scales - optimal_scale_fasta))
                primary_signal_strength_kmer = abs_coeffs_kmer[optimal_scale_idx, :]

                # Find and filter peaks/troughs
                signal_range = np.max(primary_signal_strength_kmer) - np.min(primary_signal_strength_kmer)
                prominence = signal_range * 0.10 if signal_range > 0 else 1e-5

                peak_indices, _ = find_peaks(primary_signal_strength_kmer, prominence=prominence, width=5)
                trough_indices, _ = find_peaks(-primary_signal_strength_kmer, prominence=prominence, width=5)

                # ----- START: Modified Period Calculation Code -----
                # 모든 유의미한 피크들 사이의 '연속적인' 거리를 계산하여 주기를 추정합니다.
                # 이렇게 하면 전체 신호에 걸친 반복 패턴을 더 안정적으로 나타낼 수 있습니다.
                if peak_indices.size >= 2:
                    # find_peaks는 위치 순서대로 인덱스를 반환하므로, np.diff로 연속된 피크 간의 거리를 계산합니다.
                    inter_peak_distances = np.diff(peak_indices)
                    # 계산된 거리들의 중앙값(median)을 대표적인 주기로 사용합니다.
                    median_period = np.median(inter_peak_distances)
                    arm_results['Estimated Period (bp)'] = median_period
                    print(f"Estimated median period from {len(inter_peak_distances)} peak intervals: {median_period:.1f} bp")
                else:
                    arm_results['Estimated Period (bp)'] = np.nan
                    print("Not enough peaks found to estimate period.")
                # ----- END: Modified Period Calculation Code -----


                if peak_indices.size > 0:
                    num_strong = int(round(0.5 * len(peak_indices)))
                    # 통계 분석을 위한 강한 피크 선택 로직은 신호 강도 기준으로 유지합니다.
                    strongest_peaks = peak_indices[np.argsort(primary_signal_strength_kmer[peak_indices])][::-1][:num_strong]
                else: strongest_peaks = np.array([])

                if trough_indices.size > 0:
                    num_weak = int(round(0.9 * len(trough_indices)))
                    weakest_troughs = trough_indices[np.argsort(primary_signal_strength_kmer[trough_indices])][:num_weak]
                else: weakest_troughs = np.array([])

                # Group indices and extract representative sample sequences
                peak_runs = group_contiguous_indices(strongest_peaks, min_gap=k_mer_size * 2.5)
                trough_runs = group_contiguous_indices(weakest_troughs, min_gap=k_mer_size * 2.5)

                strong_seqs = []
                for start_idx, _ in peak_runs:
                    s_start = max(0, start_idx - window_size_for_atcg_stat // 2)
                    s_end = min(len(fasta_seq_for_arm), s_start + window_size_for_atcg_stat)
                    strong_seqs.append(fasta_seq_for_arm[s_start:s_end])

                weak_seqs = []
                for start_idx, _ in trough_runs:
                    s_start = max(0, start_idx - window_size_for_atcg_stat // 2)
                    s_end = min(len(fasta_seq_for_arm), s_start + window_size_for_atcg_stat)
                    weak_seqs.append(fasta_seq_for_arm[s_start:s_end])

                # Perform statistical tests
                if strong_seqs and weak_seqs:
                    s_concat, w_concat = "".join(strong_seqs), "".join(weak_seqs)
                    at_s, cg_s = (s_concat.count('A') + s_concat.count('T')), (s_concat.count('C') + s_concat.count('G'))
                    at_w, cg_w = (w_concat.count('A') + w_concat.count('T')), (w_concat.count('C') + w_concat.count('G'))

                    arm_results.update({
                        'Strong Fasta AT': at_s, 'Strong Fasta CG': cg_s,
                        'Strong Fasta AT/CG Ratio': at_s / cg_s if cg_s > 0 else np.inf,
                        'Strong Fasta Length': len(s_concat),
                        'Weak Fasta AT': at_w, 'Weak Fasta CG': cg_w,
                        'Weak Fasta AT/CG Ratio': at_w / cg_w if cg_w > 0 else np.inf,
                        'Weak Fasta Length': len(w_concat)
                    })

                    if (at_s + cg_s > 0) and (at_w + cg_w > 0):
                        table = [[at_s, cg_s], [at_w, cg_w]]
                        _, p_chi2, _, _ = chi2_contingency(table)
                        _, p_fisher = fisher_exact(table)
                        arm_results['Chi2 P-value (FASTA)'] = p_chi2
                        arm_results['Fisher P-value (FASTA)'] = p_fisher

    # --- BigWig CWT Analysis and Correlation ---
    bw_vals = bw.values(chrom, bed_start, bed_end, numpy=True)
    if bw_vals is not None and len(bw_vals) > max(common_cwt_scales) and not np.isnan(optimal_scale_fasta):
        bw_vals_clean = np.nan_to_num(bw_vals, nan=0.0)
        coeffs_bw, _ = pywt.cwt(bw_vals_clean, common_cwt_scales, 'morl')
        abs_coeffs_bw = np.abs(coeffs_bw)

        # Get CWT signals at the optimal scale found from FASTA
        optimal_scale_idx = np.argmin(np.abs(common_cwt_scales - optimal_scale_fasta))
        bw_signal_at_optimal = abs_coeffs_bw[optimal_scale_idx, :]

        if primary_signal_strength_kmer is not None:
            # Smooth, match length, and normalize
            fasta_smooth = np.convolve(primary_signal_strength_kmer, np.ones(smooth_win)/smooth_win, mode='same')
            bw_smooth = np.convolve(bw_signal_at_optimal, np.ones(smooth_win)/smooth_win, mode='same')

            L = min(len(fasta_smooth), len(bw_smooth))
            if L > smooth_win:
                fasta_norm = (fasta_smooth[:L] - np.mean(fasta_smooth[:L])) / (np.std(fasta_smooth[:L]) + 1e-9)
                bw_norm = (bw_smooth[:L] - np.mean(bw_smooth[:L])) / (np.std(bw_smooth[:L]) + 1e-9)

                correlation = np.corrcoef(fasta_norm, bw_norm)[0, 1]
                arm_results['FASTA-BW CWT Corr (Optimal Scale)'] = correlation
                print(f"FASTA-BW CWT Correlation: {correlation:.3f}")

                # --- Visualization ---
                plt.figure(figsize=(12, 5))
                plt.plot(fasta_norm, label=f'FASTA k-mer CWT (norm)', color='dodgerblue')
                plt.plot(bw_norm, label=f'H3K9me3 bigWig CWT (norm)', color='orangered', alpha=0.7)

                # Plot peaks found in the FASTA k-mer signal
                if 'peak_indices' in locals() and peak_indices.size > 0:
                    peaks_to_plot = peak_indices[peak_indices < L]
                    plt.scatter(peaks_to_plot, fasta_norm[peaks_to_plot], color='black', s=30, marker='o', label='FASTA k-mer Peaks', zorder=5)

                title_text = f"Arm: {arm_label} - CWT Analysis at Optimal FASTA Scale ({optimal_scale_fasta:.0f})"
                if 'Estimated Period (bp)' in arm_results and not np.isnan(arm_results['Estimated Period (bp)']):
                    title_text += f"\nEstimated Median Period: {arm_results['Estimated Period (bp)']:.1f} bp"
                plt.title(title_text)

                plt.xlabel(f"Position (bp, relative to start of analysis region)")
                plt.ylabel("Normalized CWT Amplitude")
                plt.legend()
                plt.grid(True, linestyle='--', alpha=0.6)
                plt.tight_layout()
                plt.savefig(os.path.join(output_plot_dir, f"{arm_label}_cwt_analysis.png"))
                plt.close()
                print(f"Plot saved to {os.path.join(output_plot_dir, f'{arm_label}_cwt_analysis.png')}")

    results_list.append(arm_results)

# ===== Finalization =====
if bw:
    bw.close()

# Create DataFrame and save to CSV
df_results = pd.DataFrame(results_list)
# Define column order for better readability
column_order = [
    'Arm', 'Optimal Scale (FASTA k-mer)', 'Estimated Period (bp)', 'FASTA-BW CWT Corr (Optimal Scale)',
    'Strong Fasta AT', 'Strong Fasta CG', 'Strong Fasta AT/CG Ratio', 'Strong Fasta Length',
    'Weak Fasta AT', 'Weak Fasta CG', 'Weak Fasta AT/CG Ratio', 'Weak Fasta Length',
    'Chi2 P-value (FASTA)', 'Fisher P-value (FASTA)'
]
# Reorder dataframe, adding columns that might be missing in some runs as NaN
df_results = df_results.reindex(columns=column_order)
csv_output_path = "chromosome_arm_cwt_analysis_summary.csv"
df_results.to_csv(csv_output_path, index=False, na_rep='NaN')

print(f"\nAnalysis complete. Results saved to {csv_output_path}")