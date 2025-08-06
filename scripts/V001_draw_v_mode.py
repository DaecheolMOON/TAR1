# -*- coding: utf-8 -*-
"""
TAR1 주기 모드 분석:
1) Reference (CHM13) & HG002m: 염색체 팔별 SVD (U, V) 분석 및 플롯
2) HG002_karimian: read 단위 SVD (V 모드만)
3) HG002m vs HG002_karimian V 모드 비교

[실행 방법]
python analyze_tar.py \
    --ref [참조유전체 FASTA 경로] \
    --hg002 [HG002m FASTA 경로] \
    --kar [Karimian FASTA 경로] \
    --output_dir [결과 저장 폴더]
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

# --- 0. 파라미터 설정 ---
KMER_K       = 20
MIN_LEN      = 500
CWT_SCALES   = np.arange(200, 401)
CWT_WAVELET  = 'morl'
SVD_N_MODES  = 3

# --- FASTA 파싱 (arm별) ---
def parse_fasta_to_arms(fp: Path):
    arms = defaultdict(list)
    hdr, buf = None, []
    for ln in fp.read_text().splitlines():
        if ln.startswith('>'):
            if hdr and buf:
                seq = ''.join(buf).upper()
                m = re.search(r'chr([0-9XY]+[pq])', hdr)
                if m:
                    arms[m.group(1)].append(seq)
            hdr, buf = ln, []
        else:
            buf.append(ln.strip())
    if hdr and buf:
        seq = ''.join(buf).upper()
        m = re.search(r'chr([0-9XY]+[pq])', hdr)
        if m:
            arms[m.group(1)].append(seq)
    merged = {arm: ''.join(seqs) for arm, seqs in arms.items()}
    return {
        arm: seq
        for arm, seq in merged.items()
        if len(seq) >= MIN_LEN + KMER_K - 1
    }

# --- FASTA 파싱 (read별) ---
def parse_fasta_to_reads(fp: Path):
    reads, hdr, buf = {}, None, []
    for ln in fp.read_text().splitlines():
        if ln.startswith('>'):
            if hdr and buf:
                seq = ''.join(buf).upper()
                rid = extract_read_id(hdr)
                reads[rid] = seq
            hdr, buf = ln, []
        else:
            buf.append(ln.strip())
    if hdr and buf:
        seq = ''.join(buf).upper()
        rid = extract_read_id(hdr)
        reads[rid] = seq
    return {
        rid: seq
        for rid, seq in reads.items()
        if len(seq) >= MIN_LEN + KMER_K - 1
    }

def extract_read_id(header: str) -> str:
    """헤더에서 read ID 추출 (최대 20자). '::' 있으면 그 뒤, 없으면 '>' 이후 첫 공백 전까지."""
    if "::" in header:
        return header.split("::", 1)[1].split()[0][:20]
    return header[1:].split()[0][:20]

# --- k-mer 신호 & CWT 에너지 ---
def kmer_signal(seq):
    kmers = [seq[i:i + KMER_K] for i in range(len(seq) - KMER_K + 1)]
    cnt = Counter(kmers)
    return np.array([cnt[k] for k in kmers])

def cwt_energy(sig):
    if len(sig) < max(CWT_SCALES):
        return np.zeros(len(CWT_SCALES))
    coeffs, _ = pywt.cwt(sig, CWT_SCALES, CWT_WAVELET, method='fft')
    return np.sum(np.abs(coeffs), axis=1)

# --- M 행렬 구축 (누락 key → 0 벡터) ---
def build_M(seqs: dict, keys: list):
    n, s = len(keys), len(CWT_SCALES)
    M = np.zeros((n, s))
    for i, k in enumerate(keys):
        if k in seqs:
            sig = kmer_signal(seqs[k])
            M[i] = cwt_energy(sig)
    return M

# --- SVD 및 플롯 ---
def svd_and_plot(M, keys, title, output_dir: Path, plot_U=True, plot_V=True):
    U, S, Vt = svd(M, full_matrices=False)
    V = Vt.T
    var_ratio = S**2 / np.sum(S**2)

    # --- U 모드: 막대 그래프 ---
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
                label=f'U{i + 1} ({var_ratio[i] * 100:.1f}%)'
            )
        plt.title(f'{title} U modes')
        plt.xticks(indices + bar_width * (n_modes - 1) / 2, keys, rotation=90, fontsize=6)
        plt.ylabel('Contribution')
        plt.legend()
        plt.tight_layout()

        safe_title = title.replace(" ", "_").replace("/", "_")
        plt.savefig(output_dir / f'{safe_title}_U_modes.png', dpi=300)
        plt.close()

    # --- V 모드: 선 그래프 ---
    if plot_V:
        plt.figure(figsize=(8, 4))
        for i in range(min(SVD_N_MODES, V.shape[1])):
            plt.plot(CWT_SCALES, V[:, i], marker='o',
                     label=f'V{i + 1} ({var_ratio[i] * 100:.1f}%)')
        plt.title(f'{title} V modes')
        plt.xlabel('Scale (bp)')
        plt.ylabel('Contribution')
        plt.legend()
        plt.tight_layout()

        safe_title = title.replace(" ", "_").replace("/", "_")
        plt.savefig(output_dir / f'{safe_title}_V_modes.png', dpi=300)
        plt.close()

    return V

# --- 메인 실행 ---
def main():
    parser = argparse.ArgumentParser(description="TAR1 주기 모드 분석 스크립트")
    parser.add_argument("--ref", type=str, required=True, help="참조 유전체(CHM13) FASTA 파일 경로")
    parser.add_argument("--hg002", type=str, required=True, help="HG002m FASTA 파일 경로")
    parser.add_argument("--kar", type=str, required=True, help="HG002_karimian FASTA 파일 경로")
    parser.add_argument("--output_dir", type=str, required=True, help="결과 그래프를 저장할 폴더 경로")
    args = parser.parse_args()

    REF_FASTA  = Path(args.ref)
    HG002_FASTA = Path(args.hg002)
    KAR_FASTA  = Path(args.kar)
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Reference & HG002m (arm별)
    ref_seqs = parse_fasta_to_arms(REF_FASTA)
    hg2_seqs = parse_fasta_to_arms(HG002_FASTA)
    arms = sorted(set(ref_seqs) | set(hg2_seqs))

    print("Arms:", arms)
    M_ref = build_M(ref_seqs, arms)
    M_hg2 = build_M(hg2_seqs, arms)

    V_ref = svd_and_plot(M_ref, arms, "CHM13",  OUTPUT_DIR, plot_U=True, plot_V=True)
    V_hg2 = svd_and_plot(M_hg2, arms, "HG002m", OUTPUT_DIR, plot_U=True, plot_V=True)

    # 2) HG002_karimian (read별)
    kar_seqs = parse_fasta_to_reads(KAR_FASTA)
    reads = sorted(kar_seqs)
    print("Reads count:", len(reads))
    M_kar = build_M(kar_seqs, reads)
    V_kar = svd_and_plot(M_kar, reads, "HG002_karimian", OUTPUT_DIR,
                         plot_U=False, plot_V=True)

    # 3) V 모드 비교 (HG002m vs KAR)
    print("\nV-mode cosine similarities (HG002m vs KAR):")
    for i in range(min(V_hg2.shape[1], V_kar.shape[1], SVD_N_MODES)):
        sim = 1 - cosine(V_hg2[:, i], V_kar[:, i])
        print(f" Mode{i + 1}: {sim:.4f}")

if __name__ == "__main__":
    main()
