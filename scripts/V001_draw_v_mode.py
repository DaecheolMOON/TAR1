# -*- coding: utf-8 -*-
"""
TAR1 periodic mode analysis:

1) Reference (CHM13) & HG002m:
   - SVD (U, V) analysis and plots per chromosome arm.
2) HG002_karimian:
   - SVD per read (V modes only).
3) V-mode comparison between HG002m and HG002_karimian
   using cosine similarity.

Example:
    python V001_draw_v_mode.py \
        --ref <reference FASTA> \
        --hg002 <HG002m FASTA> \
        --kar <Karimian FASTA> \
        --output_dir <output directory>
"""

import re
import argparse
from pathlib import Path
import numpy as np
import pywt
from scipy.linalg import svd
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine

# --- 0. Parameters ---
KMER_K      = 20
MIN_LEN     = 500
CWT_SCALES  = np.arange(200, 401)
CWT_WAVELET = "morl"
SVD_N_MODES = 3


# --- FASTA parsing (by chromosome arm) ---
def parse_fasta_to_arms(fp: Path):
    """
    Parse a FASTA file and merge sequences per chromosome arm.

    Arm names are extracted from headers using the pattern:
        chr([0-9XY]+[pq])
    Multiple sequences with the same arm are concatenated.
    Only arms with length >= MIN_LEN + KMER_K - 1 are kept.
    """
    arms = defaultdict(list)
    hdr, buf = None, []
    for ln in fp.read_text().splitlines():
        if ln.startswith(">"):
            if hdr and buf:
                seq = "".join(buf).upper()
                m = re.search(r"chr([0-9XY]+[pq])", hdr)
                if m:
                    arms[m.group(1)].append(seq)
            hdr, buf = ln, []
        else:
            buf.append(ln.strip())
    if hdr and buf:
        seq = "".join(buf).upper()
        m = re.search(r"chr([0-9XY]+[pq])", hdr)
        if m:
            arms[m.group(1)].append(seq)

    merged = {arm: "".join(seqs) for arm, seqs in arms.items()}
    return {
        arm: seq
        for arm, seq in merged.items()
        if len(seq) >= MIN_LEN + KMER_K - 1
    }


# --- FASTA parsing (by read) ---
def parse_fasta_to_reads(fp: Path):
    """
    Parse a FASTA file into a dictionary of read_id -> sequence.

    Read IDs are extracted from headers using extract_read_id().
    Only reads with length >= MIN_LEN + KMER_K - 1 are kept.
    """
    reads, hdr, buf = {}, None, []
    for ln in fp.read_text().splitlines():
        if ln.startswith(">"):
            if hdr and buf:
                seq = "".join(buf).upper()
                rid = extract_read_id(hdr)
                reads[rid] = seq
            hdr, buf = ln, []
        else:
            buf.append(ln.strip())
    if hdr and buf:
        seq = "".join(buf).upper()
        rid = extract_read_id(hdr)
        reads[rid] = seq
    return {
        rid: seq
        for rid, seq in reads.items()
        if len(seq) >= MIN_LEN + KMER_K - 1
    }


def extract_read_id(header: str) -> str:
    """
    Extract a read ID (up to 20 characters) from a FASTA header.

    If '::' is present, use the substring after the first '::' up to
    the first whitespace. Otherwise, use the substring after '>'
    up to the first whitespace.
    """
    if "::" in header:
        return header.split("::", 1)[1].split()[0][:20]
    return header[1:].split()[0][:20]


# --- k-mer signal & CWT energy ---
def kmer_signal(seq: str) -> np.ndarray:
    """
    Build a k-mer repetitiveness signal from a sequence.

    signal[i] = count of k-mer seq[i:i+KMER_K] in the entire sequence.
    """
    kmers = [seq[i : i + KMER_K] for i in range(len(seq) - KMER_K + 1)]
    cnt = Counter(kmers)
    return np.array([cnt[k] for k in kmers], dtype=float)


def cwt_energy(sig: np.ndarray) -> np.ndarray:
    """
    Compute CWT energy per scale for a given signal.

    Returns a vector of length len(CWT_SCALES), where each element is
    the sum of absolute CWT coefficients across positions at that scale.
    """
    if len(sig) < max(CWT_SCALES):
        return np.zeros(len(CWT_SCALES), dtype=float)
    coeffs, _ = pywt.cwt(sig, CWT_SCALES, CWT_WAVELET, method="fft")
    return np.sum(np.abs(coeffs), axis=1)


# --- Build M matrix (missing keys -> zero vectors) ---
def build_M(seqs: dict, keys: list) -> np.ndarray:
    """
    Build an M matrix of shape (len(keys), len(CWT_SCALES)).

    Each row corresponds to one key (arm or read ID),
    and contains the CWT energy per scale.
    Keys missing in 'seqs' are assigned a zero vector.
    """
    n, s = len(keys), len(CWT_SCALES)
    M = np.zeros((n, s), dtype=float)
    for i, k in enumerate(keys):
        if k in seqs:
            sig = kmer_signal(seqs[k])
            M[i] = cwt_energy(sig)
    return M


# --- SVD and plotting ---
def svd_and_plot(M, keys, title, output_dir: Path, plot_U=True, plot_V=True):
    """
    Perform SVD on M and plot U and/or V modes.

    U modes:
      - bar plots across keys (arms or reads).
    V modes:
      - line plots across scales.

    The explained variance ratio (S^2 / sum S^2) is shown in the legend.
    """
    U, S, Vt = svd(M, full_matrices=False)
    V = Vt.T
    var_ratio = S**2 / np.sum(S**2)

    # --- U modes: bar plots ---
    if plot_U:
        n_modes = min(SVD_N_MODES, U.shape[1])
        indices = np.arange(len(keys))
        bar_width = 0.8 / n_modes

        plt.figure(figsize=(12, 6))
        for i in range(n_modes):
            offset = indices + i * bar_width
            plt.bar(
                offset,
                U[:, i],
                width=bar_width,
                label=f"U{i + 1} ({var_ratio[i] * 100:.1f}%)",
            )
        plt.title(f"{title} U modes")
        plt.xticks(
            indices + bar_width * (n_modes - 1) / 2,
            keys,
            rotation=90,
            fontsize=6,
        )
        plt.ylabel("Contribution")
        plt.legend()
        plt.tight_layout()

        safe_title = title.replace(" ", "_").replace("/", "_")
        plt.savefig(output_dir / f"{safe_title}_U_modes.png", dpi=300)
        plt.close()

    # --- V modes: line plots ---
    if plot_V:
        plt.figure(figsize=(8, 4))
        for i in range(min(SVD_N_MODES, V.shape[1])):
            plt.plot(
                CWT_SCALES,
                V[:, i],
                marker="o",
                label=f"V{i + 1} ({var_ratio[i] * 100:.1f}%)",
            )
        plt.title(f"{title} V modes")
        plt.xlabel("Scale (bp)")
        plt.ylabel("Contribution")
        plt.legend()
        plt.tight_layout()

        safe_title = title.replace(" ", "_").replace("/", "_")
        plt.savefig(output_dir / f"{safe_title}_V_modes.png", dpi=300)
        plt.close()

    return V


# --- Main ---
def main():
    parser = argparse.ArgumentParser(
        description="TAR1 periodic mode analysis based on CWT energy and SVD."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Path to reference genome (CHM13) FASTA file.",
    )
    parser.add_argument(
        "--hg002",
        type=str,
        required=True,
        help="Path to HG002m FASTA file (arm-level sequences).",
    )
    parser.add_argument(
        "--kar",
        type=str,
        required=True,
        help="Path to HG002_karimian FASTA file (read-level sequences).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the result plots.",
    )
    args = parser.parse_args()

    REF_FASTA = Path(args.ref)
    HG002_FASTA = Path(args.hg002)
    KAR_FASTA = Path(args.kar)
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Reference & HG002m (by arm)
    ref_seqs = parse_fasta_to_arms(REF_FASTA)
    hg2_seqs = parse_fasta_to_arms(HG002_FASTA)
    arms = sorted(set(ref_seqs) | set(hg2_seqs))

    print("Arms:", arms)
    M_ref = build_M(ref_seqs, arms)
    M_hg2 = build_M(hg2_seqs, arms)

    V_ref = svd_and_plot(M_ref, arms, "CHM13", OUTPUT_DIR, plot_U=True, plot_V=True)
    V_hg2 = svd_and_plot(M_hg2, arms, "HG002m", OUTPUT_DIR, plot_U=True, plot_V=True)

    # 2) HG002_karimian (by read)
    kar_seqs = parse_fasta_to_reads(KAR_FASTA)
    reads = sorted(kar_seqs)
    print("Reads count:", len(reads))
    M_kar = build_M(kar_seqs, reads)
    V_kar = svd_and_plot(
        M_kar, reads, "HG002_karimian", OUTPUT_DIR, plot_U=False, plot_V=True
    )

    # 3) V-mode comparison (HG002m vs KAR) using cosine similarity
    print("\nV-mode cosine similarities (HG002m vs KAR):")
    for i in range(min(V_hg2.shape[1], V_kar.shape[1], SVD_N_MODES)):
        sim = 1 - cosine(V_hg2[:, i], V_kar[:, i])
        print(f" Mode{i + 1}: {sim:.4f}")


if __name__ == "__main__":
    main()
