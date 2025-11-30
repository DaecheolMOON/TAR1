# -*- coding: utf-8 -*-
"""
Primate TAR1 periodic mode analysis (V001_draw_v_mode_primate.py):

For each primate FASTA file containing TAR1 sequence blocks,
this script:
  - builds a k-mer-based CWT energy matrix (blocks × scales),
  - performs SVD on this matrix,
  - and generates U-mode (block contributions) and V-mode (scale patterns) plots.

Usage:
    python V001_draw_v_mode_primate.py \
        --borang <Bornean orangutan FASTA> \
        --bonobo <bonobo FASTA> \
        --chimp <chimpanzee FASTA> \
        --output_dir <output directory>
"""

import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pywt
from scipy.linalg import svd
import matplotlib.pyplot as plt

# --- Parameters (same as human-mode script) ---
KMER_K = 20
MIN_LEN = 500  # minimum sequence length
CWT_SCALES = np.arange(200, 401)
CWT_WAVELET = "morl"
SVD_N_MODES = 3


# --- FASTA parsing (per sequence block) ---
def parse_fasta_to_sequences(fp: Path):
    """
    Read a FASTA file and return a dict {block_id: sequence}.

    Only sequences with length >= MIN_LEN + KMER_K - 1 are retained.
    Block IDs are derived from headers using extract_id_from_header().
    """
    sequences, hdr, buf = {}, None, []
    try:
        lines = fp.read_text().splitlines()
    except Exception as e:
        print(f"Error reading file {fp}: {e}")
        return {}

    for ln in lines:
        if ln.startswith(">"):
            if hdr and buf:
                seq = "".join(buf).upper()
                if len(seq) >= MIN_LEN + KMER_K - 1:
                    rid = extract_id_from_header(hdr)
                    sequences[rid] = seq
            hdr, buf = ln, []
        else:
            buf.append(ln.strip())
    if hdr and buf:
        seq = "".join(buf).upper()
        if len(seq) >= MIN_LEN + KMER_K - 1:
            rid = extract_id_from_header(hdr)
            sequences[rid] = seq
    return sequences


def extract_id_from_header(header: str) -> str:
    """
    Extract a block ID from a FASTA header.

    If '::' is present, use the substring after '::' (e.g. '::chr8:13587-15321'
    -> 'chr8:13587-15321'). Otherwise, use the first token after '>'.
    The ID is truncated to at most 30 characters.
    """
    if "::" in header:
        try:
            return header.split("::", 1)[1].strip()
        except IndexError:
            return header[1:].split()[0][:30]  # fallback
    return header[1:].split()[0][:30]


# --- k-mer signal & CWT energy ---
def kmer_signal(seq: str) -> np.ndarray:
    """
    Build a k-mer repetitiveness signal.

    signal[i] = count of k-mer seq[i:i+KMER_K] in the entire sequence.
    """
    kmers = [seq[i : i + KMER_K] for i in range(len(seq) - KMER_K + 1)]
    cnt = Counter(kmers)
    return np.array([cnt[k] for k in kmers], dtype=float)


def cwt_energy(sig: np.ndarray) -> np.ndarray:
    """
    Compute CWT energy per scale for a given signal.

    Returns a vector of length len(CWT_SCALES).
    """
    if len(sig) < max(CWT_SCALES):
        return np.zeros(len(CWT_SCALES), dtype=float)
    coeffs, _ = pywt.cwt(sig, CWT_SCALES, CWT_WAVELET, method="fft")
    return np.sum(np.abs(coeffs), axis=1)


# --- Build M matrix ---
def build_M(seqs: dict, keys: list) -> np.ndarray:
    """
    Build an M matrix of shape (len(keys), len(CWT_SCALES)).

    Each row contains the CWT energy per scale for one sequence block.
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

    - U modes: bar plots showing contributions of each sequence block.
    - V modes: line plots showing scale-dependent periodicity patterns.
    """
    if M.size == 0 or np.all(M == 0):
        print(f"  Warning: no valid data for '{title}', skipping SVD and plots.")
        return None

    U, S, Vt = svd(M, full_matrices=False)
    V = Vt.T
    var_ratio = S**2 / np.sum(S**2)
    safe_title = title.replace(" ", "_").replace("/", "_")

    # --- U modes: bar plots ---
    if plot_U:
        n_modes = min(SVD_N_MODES, U.shape[1])
        indices = np.arange(len(keys))
        bar_width = 0.8 / n_modes

        plt.figure(figsize=(16, 8))
        for i in range(n_modes):
            offset = indices + i * bar_width
            plt.bar(
                offset,
                U[:, i],
                width=bar_width,
                label=f"U{i + 1} ({var_ratio[i] * 100:.1f}%)",
            )

        plt.title(
            f"{title} U modes (Contribution of each sequence block)", fontsize=16
        )
        plt.xticks(
            indices + bar_width * (n_modes - 1) / 2,
            keys,
            rotation=90,
            fontsize=6,
        )
        plt.ylabel("Contribution", fontsize=12)
        plt.legend()
        plt.grid(axis="y", linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.savefig(output_dir / f"{safe_title}_U_modes.png", dpi=300)
        plt.close()

    # --- V modes: line plots ---
    if plot_V:
        plt.figure(figsize=(10, 6))
        for i in range(min(SVD_N_MODES, V.shape[1])):
            plt.plot(
                CWT_SCALES,
                V[:, i],
                marker="o",
                markersize=4,
                linestyle="-",
                label=f"V{i + 1} ({var_ratio[i] * 100:.1f}%)",
            )

        plt.title(f"{title} V modes (Periodicity patterns)", fontsize=16)
        plt.xlabel("Scale (Periodicity in bp)", fontsize=12)
        plt.ylabel("Contribution", fontsize=12)
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.savefig(output_dir / f"{safe_title}_V_modes.png", dpi=300)
        plt.close()

    return V


# --- Main ---
def main():
    parser = argparse.ArgumentParser(
        description="Primate TAR1 periodic mode analysis (CWT energy + SVD)."
    )
    parser.add_argument(
        "--borang",
        type=str,
        required=True,
        help="Path to Bornean orangutan TAR1 FASTA file.",
    )
    parser.add_argument(
        "--bonobo",
        type=str,
        required=True,
        help="Path to bonobo TAR1 FASTA file.",
    )
    parser.add_argument(
        "--chimp",
        type=str,
        required=True,
        help="Path to chimpanzee TAR1 FASTA file.",
    )
    # Additional primates can be added similarly:
    # parser.add_argument("--gorilla", type=str, help="Path to gorilla TAR1 FASTA file.")
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the result plots.",
    )
    args = parser.parse_args()

    species_files = {
        "borang": Path(args.borang),
        "bonobo": Path(args.bonobo),
        "chimp": Path(args.chimp),
    }

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("### Primate TAR1 periodicity analysis started ###")

    for species, fasta_path in species_files.items():
        print(f"\n--- Analyzing species: {species.capitalize()} ---")

        if not fasta_path.is_file():
            print(f"  Warning: file not found, skipping: {fasta_path}")
            continue

        # 1. Parse FASTA
        sequences = parse_fasta_to_sequences(fasta_path)
        keys = sorted(sequences.keys())

        if not keys:
            print(
                f"  Warning: no valid sequence blocks (>= {MIN_LEN} bp) found for {species}."
            )
            continue

        print(f"  > Found {len(keys)} TAR1 sequence blocks.")

        # 2. Build M matrix
        M = build_M(sequences, keys)

        # 3. SVD and plot
        print("  > Running SVD and generating plots...")
        svd_and_plot(
            M, keys, f"Primate_{species}", output_dir, plot_U=True, plot_V=True
        )
        print(f"  > Done. Results saved to '{output_dir}'.")

    print("\n### All primate analyses completed. ###")


if __name__ == "__main__":
    main()
