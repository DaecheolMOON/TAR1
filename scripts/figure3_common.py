from __future__ import annotations

import os
import re
from pathlib import Path
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pywt
from scipy.signal import find_peaks
from Bio import SeqIO


REFERENCE_FASTA_MAP: Dict[str, str] = {
    "human": os.path.join("reference", "chm13v2.0.fa"),
    "chimp": os.path.join("reference", "GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"),
    "bonobo": os.path.join("reference", "GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna"),
    "gorilla": os.path.join("reference", "GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna"),
    "borang": os.path.join("reference", "GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna"),
    "sorang": os.path.join("reference", "GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna"),
}

TAR1_FASTA_MAP: Dict[str, str] = {
    "human": os.path.join("tar1_blocks", "human_CHM13", "tar1_blocks.fa"),
    "chimp": os.path.join("tar1_blocks", "chimp", "tar1_blocks.fa"),
    "bonobo": os.path.join("tar1_blocks", "bonobo", "tar1_blocks.fa"),
    "gorilla": os.path.join("tar1_blocks", "gorilla", "tar1_blocks.fa"),
    "borang": os.path.join("tar1_blocks", "borang", "tar1_blocks.fa"),
    "sorang": os.path.join("tar1_blocks", "sorang", "tar1_blocks.fa"),
}

CUTOFFS_DIST = [
    10_000, 20_000, 30_000, 50_000, 70_000, 100_000, 200_000, 500_000,
    1_000_000, 2_000_000, 3_000_000, 5_000_000, 7_000_000, 8_000_000,
    9_000_000, 10_000_000, 20_000_000, 30_000_000, 40_000_000, 50_000_000,
    60_000_000, 70_000_000, 80_000_000, 90_000_000, 100_000_000
]

CHIMP_CAP_PLUS = {
    "2p", "2q", "7p", "8p", "10p", "11p", "12p", "13p", "15p", "16p", "16q", "17p", "17q",
    "18p", "19p", "19q", "20p", "20q", "21q", "22p", "22q", "Xq", "14p", "6p"
}

TOTAL_ARMS = {
    "human": 46,
    "chimp": 48,
    "bonobo": 48,
    "gorilla": 48,
    "borang": 48,
    "sorang": 48,
}

HDR_RE = re.compile(r"^(chr[\w]+[pq]):(\d+)-(\d+)")
CHR_ARM_RE = re.compile(r"^(chr(?:[0-9]+|X|Y))([pq])$", flags=re.IGNORECASE)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sanitize_id_for_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def load_manifest(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def infer_chr_only(arm: str) -> Optional[str]:
    m = CHR_ARM_RE.match(str(arm).split()[0])
    if m:
        return m.group(1)
    return None


def is_q_arm(arm: str) -> bool:
    return str(arm).strip().lower().endswith("q")


def normalize_chr_token(token: str) -> str:
    t = str(token).strip().lower()
    t = t.replace("chromosome", "").strip()
    t = t.replace("chr", "").strip()
    return t


def chromosome_matches_record(rec, chr_id: str) -> bool:
    target = normalize_chr_token(chr_id)
    texts = [
        str(rec.id),
        str(getattr(rec, "name", "")),
        str(getattr(rec, "description", "")),
    ]

    normalized = set()
    for text in texts:
        if not text:
            continue

        normalized.add(normalize_chr_token(text))

        m1 = re.search(r"\bchromosome\s+([0-9XY]+)\b", text, flags=re.IGNORECASE)
        if m1:
            normalized.add(normalize_chr_token(m1.group(1)))

        m2 = re.search(r"\bchr([0-9XY]+)\b", text, flags=re.IGNORECASE)
        if m2:
            normalized.add(normalize_chr_token(m2.group(1)))

        if text.strip().lower() in {f"chr{target}", target}:
            normalized.add(target)

    return target in normalized


def resolve_data_root(data_root: str | Path | None = None) -> Path:
    if data_root is None:
        data_root = os.environ.get("TAR1_DATA_ROOT", "./data")
    return Path(data_root)


def get_reference_fasta(species: str, data_root: str | Path | None = None) -> Path:
    if species not in REFERENCE_FASTA_MAP:
        raise KeyError(f"Unknown species: {species}")
    root = resolve_data_root(data_root)
    p = root / REFERENCE_FASTA_MAP[species]
    if not p.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {p}")
    return p


def get_tar1_fasta(species: str, data_root: str | Path | None = None) -> Path:
    if species not in TAR1_FASTA_MAP:
        raise KeyError(f"Unknown species: {species}")
    root = resolve_data_root(data_root)
    p = root / TAR1_FASTA_MAP[species]
    if not p.exists():
        raise FileNotFoundError(f"TAR1 FASTA not found: {p}")
    return p


def get_chromosome_sequence(
    species: str,
    chrom: str,
    data_root: str | Path | None = None,
) -> Tuple[str, str]:
    ref_path = get_reference_fasta(species, data_root=data_root)

    matched_desc = None
    matched_seq = None
    for rec in SeqIO.parse(str(ref_path), "fasta"):
        if chromosome_matches_record(rec, chrom):
            matched_desc = str(rec.description)
            matched_seq = str(rec.seq).upper()
            break

    if matched_seq is None:
        previews = []
        for i, rec in enumerate(SeqIO.parse(str(ref_path), "fasta")):
            if i >= 5:
                break
            previews.append(f"id={rec.id} | desc={rec.description}")
        raise ValueError(
            f"Chromosome {chrom} not found in {ref_path}. "
            f"First records: {' || '.join(previews)}"
        )

    return matched_desc, matched_seq


def build_telomere_window_sequence(
    species: str,
    arm: str,
    window_bp: int,
    data_root: str | Path | None = None,
    orient_telomere_left: bool = False,
) -> str:
    chrom = infer_chr_only(arm)
    if chrom is None:
        raise ValueError(f"Could not infer chromosome from arm={arm}")

    _, full_seq = get_chromosome_sequence(species, chrom, data_root=data_root)
    L = len(full_seq)
    w = min(int(window_bp), L)

    if is_q_arm(arm):
        seq = full_seq[L - w:]
        if orient_telomere_left:
            seq = seq[::-1]
    else:
        seq = full_seq[:w]

    return seq


def parse_tar1_fasta(fasta_path: str | Path) -> List[Tuple[str, int, int]]:
    items: List[Tuple[str, int, int]] = []
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        return items

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip().split()[0]
            m = HDR_RE.match(header)
            if m:
                arm = m.group(1)
                s = int(m.group(2))
                e = int(m.group(3))
                if s > e:
                    s, e = e, s
                items.append((arm, s, e))
    return items


def tar1_items_by_arm(fasta_path: str | Path) -> Dict[str, List[Tuple[int, int]]]:
    d: Dict[str, List[Tuple[int, int]]] = {}
    for arm, s, e in parse_tar1_fasta(fasta_path):
        d.setdefault(arm, []).append((s, e))
    for arm in d:
        d[arm] = sorted(d[arm], key=lambda x: (x[0], x[1]))
    return d


def tar1_intervals_within_window(
    species: str,
    arm: str,
    window_bp: int,
    data_root: str | Path | None = None,
) -> List[Tuple[int, int]]:
    tar1_path = get_tar1_fasta(species, data_root=data_root)
    d = tar1_items_by_arm(tar1_path)
    intervals = d.get(arm, [])
    out = []
    for s, e in intervals:
        if s <= window_bp and e <= window_bp:
            out.append((s, e))
    return out


def kmer_repetitiveness_signal(seq: str, k: int = 20) -> np.ndarray:
    seq = seq.upper()
    if len(seq) < k:
        return np.array([], dtype=float)
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    freq = Counter(kmers)
    return np.array([freq[km] for km in kmers], dtype=float)


def group_contiguous_indices(indices, min_gap=10):
    indices = np.array(indices, dtype=int)
    if indices.size == 0:
        return []
    sidx = np.sort(indices)
    groups = []
    start = end = int(sidx[0])
    for x in sidx[1:]:
        x = int(x)
        if x - end <= min_gap:
            end = x
        else:
            groups.append((start, end))
            start = end = x
    groups.append((start, end))
    return groups


def atcg_ratio_calc(sequence: str, indices, window_size: int = 25) -> np.ndarray:
    seq = sequence.upper()
    L = len(seq)
    ratios = []
    for idx in indices:
        st = max(0, int(idx) - window_size // 2)
        en = min(L, int(idx) + window_size // 2 + (window_size % 2))
        wseq = seq[st:en]
        if not wseq:
            ratios.append(np.nan)
            continue
        at_count = wseq.count("A") + wseq.count("T")
        cg_count = wseq.count("C") + wseq.count("G")
        if cg_count == 0:
            ratios.append(np.inf if at_count > 0 else np.nan)
        else:
            ratios.append(at_count / cg_count)
    return np.array(ratios, dtype=float)


def find_optimal_scale(abs_coeffs_matrix, scales_array, peak_prom_frac=0.05, width_for_scan=3):
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


def analyze_one_record_cwt(
    seq_id: str,
    seq_str: str,
    k=20,
    scale_min=200,
    scale_max=401,
    peak_keep_frac=0.5,
    prominence_frac=0.10,
    width_threshold=5,
):
    events = []
    signal = kmer_repetitiveness_signal(seq_str, k=k)
    if signal.size == 0:
        return None, events, None, None

    scales = np.arange(scale_min, scale_max, dtype=int)
    coeffs, _ = pywt.cwt(signal, scales, "morl")
    abs_coeffs = np.abs(coeffs)

    opt_scale = find_optimal_scale(abs_coeffs, scales)
    scale_idx = int(np.argmin(np.abs(scales - opt_scale)))
    primary = abs_coeffs[scale_idx, :]

    rmin, rmax = float(np.min(primary)), float(np.max(primary))
    rrange = rmax - rmin
    if rrange <= 0:
        return float(scales[scale_idx]), events, signal, abs_coeffs

    prom = rrange * prominence_frac
    if prom == 0:
        prom = float(np.mean(primary)) * 0.01

    peaks_all, _ = find_peaks(primary, prominence=prom, width=width_threshold)
    troughs_all, _ = find_peaks(-primary, prominence=prom, width=width_threshold)

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

    gap_for_grouping = int(k * 2.5)
    peak_runs = group_contiguous_indices(peaks, min_gap=gap_for_grouping)
    trough_runs = group_contiguous_indices(troughs, min_gap=gap_for_grouping)

    window_size_for_atcg = k * 2

    def process_runs(runs, signal_type):
        for (start_idx, end_idx) in runs:
            kmer_indices = np.arange(start_idx, end_idx + 1, dtype=int)
            atcg_vals = atcg_ratio_calc(seq_str, kmer_indices, window_size=window_size_for_atcg)
            atcg_vals_clean = atcg_vals[~np.isnan(atcg_vals) & ~np.isinf(atcg_vals)]
            mean_atcg = float(np.mean(atcg_vals_clean)) if atcg_vals_clean.size > 0 else np.nan

            start_base = int(start_idx)
            end_base = int(end_idx + k - 1)

            events.append({
                "record_id": seq_id,
                "signal_type": signal_type,
                "start_kmer_idx": int(start_idx),
                "end_kmer_idx": int(end_idx),
                "start_base": start_base,
                "end_base": end_base,
                "run_kmer_len": int(end_idx - start_idx + 1),
                "scale_used": float(scales[scale_idx]),
                "mean_atcg": mean_atcg,
            })

    process_runs(peak_runs, "strong")
    process_runs(trough_runs, "weak")
    return float(scales[scale_idx]), events, signal, abs_coeffs


def run_counting(fasta_path: str | Path, cutoffs: List[int]) -> Tuple[List[int], int]:
    items = parse_tar1_fasta(fasta_path)
    counts = []
    for c in cutoffs:
        found_arms = set()
        for arm, s, e in items:
            if s <= c and e <= c:
                found_arms.add(arm)
        counts.append(len(found_arms))
    return counts, len(items)


def build_universe() -> set:
    s = {f"{i}{a}" for i in range(1, 23) for a in ("p", "q")}
    s |= {"Xp", "Xq"}
    return s


def normalize_arm_token(tok: str) -> str:
    s = str(tok).strip()
    s = re.sub(r"[,\.;:\s]+$", "", s)
    s = s.replace(" ", "").lower()
    if s.startswith("chr"):
        s = s[3:]
    m = re.match(r"^([0-9]+|x)([pq])$", s)
    if not m:
        raise ValueError(f"Cannot normalize arm token: {tok}")
    chrom = "X" if m.group(1) == "x" else str(int(m.group(1)))
    return f"{chrom}{m.group(2)}"