import numpy as np
import pywt
import matplotlib.pyplot as plt
from collections import Counter
from numpy.linalg import svd
import os
import argparse
from Bio import SeqIO

def analyze_sequence(seq, read_id, output_dir):
    """
    Analyze a single FASTA read: generate k-mer signal, compute CWT scalogram, 
    perform SVD on raw coefficients, and save plots.
    Usage example:
      python fasta_to_scalogram_svd.py \
        --input_fasta path/to/reads.fasta \
        --output_dir path/to/output
    """
    print(f"--- Processing read: {read_id} ---")

    # 1. k-mer signal generation (k=20)
    k = 20
    if len(seq) < k:
        print(f"Warning: sequence '{read_id}' shorter than k={k}. Skipping.")
        return
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    freq = Counter(kmers)
    signal = np.array([freq[kmer] for kmer in kmers])

    # 2. Continuous Wavelet Transform (scales 100-400)
    scales = np.arange(100, 401)
    coeffs, freqs = pywt.cwt(signal, scales, 'morl')

    # 3. Plot and save scalogram
    plt.figure(figsize=(10, 6))
    plt.imshow(np.abs(coeffs), aspect='auto', interpolation='nearest', origin='lower',
               extent=[0, len(signal), scales[0], scales[-1]])
    plt.xlabel('Sequence Position')
    plt.ylabel('Wavelet Scale')
    plt.title(f'Scalogram of Repeat Intensity for {read_id}')
    scalogram_path = os.path.join(output_dir, f"{read_id}_scalogram.png")
    plt.savefig(scalogram_path)
    plt.close()
    print(f"Saved scalogram to: {scalogram_path}")

    # 4. Prepare data matrix for SVD (use magnitude of coefficients)
    #    F matrix of shape (num_scales, signal_length)
    F = np.abs(coeffs)

    # 5. Perform SVD: F = U * S * Vt
    U, S, Vt = svd(F, full_matrices=False)

    # 6. Plot top-3 basis functions (scale modes)
    plt.figure(figsize=(8, 4))
    for i in range(min(3, len(S))):
        plt.plot(scales, U[:, i], label=f'Basis {i+1} (σ={S[i]:.1f})')
    plt.xlabel('Wavelet Scale (100–400)')
    plt.ylabel('Basis Value')
    plt.title(f'Top-3 Basis Functions from SVD for {read_id}')
    plt.legend()
    basis_path = os.path.join(output_dir, f"{read_id}_svd_basis.png")
    plt.savefig(basis_path)
    plt.close()
    print(f"Saved SVD basis functions to: {basis_path}")

    # 7. Plot top-3 components (position modes)
    x = np.arange(len(signal))
    plt.figure(figsize=(8, 4))
    for i in range(min(3, len(S))):
        plt.plot(x, Vt[i, :], label=f'Component {i+1} (σ={S[i]:.1f})')
    plt.xlabel('Sequence Position')
    plt.ylabel('Component Value')
    plt.title(f'Top-3 Components from SVD for {read_id}')
    plt.legend()
    components_path = os.path.join(output_dir, f"{read_id}_svd_components.png")
    plt.savefig(components_path)
    plt.close()
    print(f"Saved SVD components to: {components_path}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze repeats in FASTA sequences using CWT and SVD, then save plots."
    )
    parser.add_argument(
        "--input_fasta",
        type=str,
        required=True,
        help="Path to input FASTA file containing one or more reads."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where output plots will be saved."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory ensured at: {args.output_dir}")

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        read_id = record.id
        seq = str(record.seq).upper().replace("\n", "").replace(" ", "")
        analyze_sequence(seq, read_id, args.output_dir)

if __name__ == "__main__":
    main()
