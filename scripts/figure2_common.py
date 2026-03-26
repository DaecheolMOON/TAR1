from __future__ import annotations

import os
import re
from pathlib import Path
from collections import Counter
from typing import Optional, Iterable

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.signal import find_peaks


# ============================================================
# Shared plotting / distribution constants
# ============================================================

CUTOFFS_DIST = [
    10_000, 20_000, 30_000, 50_000, 70_000, 100_000, 200_000, 500_000,
    1_000_000, 2_000_000, 3_000_000, 5_000_000, 7_000_000, 8_000_000,
    9_000_000, 10_000_000, 20_000_000, 30_000_000, 40_000_000, 50_000_000,
    60_000_000, 70_000_000, 80_000_000, 90_000_000, 100_000_000,
]

TOTAL_ARMS = {
    "human": 46,
    "chimp": 48,
    "bonobo": 48,
    "gorilla": 48,
    "borang": 48,
    "sorang": 48,
}

STYLE_MAP = {
    "human":   {"color": "black",       "marker": "X", "style": "-",  "label": "Human"},
    "chimp":   {"color": "forestgreen", "marker": "s", "style": ":",  "label": "Chimpanzee"},
    "bonobo":  {"color": "royalblue",   "marker": "P", "style": "--", "label": "Bonobo"},
    "gorilla": {"color": "crimson",     "marker": "D", "style": "-.", "label": "Gorilla"},
    "borang":  {"color": "gold",        "marker": "o", "style": "-",  "label": "Bornean Orangutan"},
    "sorang":  {"color": "darkorange",  "marker": "^", "style": "--", "label": "Sumatran Orangutan"},
}


# ============================================================
# Basic paths / repo helpers
# ============================================================

def get_repo_root() -> Path:
    here = Path(__file__).resolve()
    for p in [here.parent] + list(here.parents):
        if (p / "scripts").exists() and (p / "metadata").exists():
            return p
    return here.parent.parent


def get_data_root(data_root: Optional[str] = None) -> Path:
    if data_root is not None:
        return Path(data_root).resolve()

    env = os.environ.get("TAR1_DATA_ROOT", None)
    if env:
        return Path(env).resolve()

    return (get_repo_root() / "data").resolve()


def ensure_dir(path: Path | str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sanitize_id_for_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def write_single_fasta(header: str, sequence: str, out_path: Path | str, wrap: int = 60) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sequence = str(sequence)
    with open(out_path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), wrap):
            f.write(sequence[i:i + wrap] + "\n")


def read_single_fasta(path: Path | str) -> tuple[str, str]:
    path = Path(path)
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) == 0:
        raise ValueError(f"No FASTA records found: {path}")
    if len(records) > 1:
        raise ValueError(f"Expected a single FASTA record, found {len(records)} in {path}")
    rec = records[0]
    return rec.id, str(rec.seq).upper().replace("\n", "").replace(" ", "")


# ============================================================
# Manifest
# ============================================================

def load_manifest(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise ValueError(f"Manifest is empty: {path}")

    return df


# ============================================================
# Reference FASTA resolution
# ============================================================

REFERENCE_FILES = {
    "human": "chm13v2.0.fa",
    "human_chm13": "chm13v2.0.fa",
    "chimp": "GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna",
    "bonobo": "GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna",
    "gorilla": "GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna",
    "borang": "GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna",
    "sorang": "GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna",
}


def normalize_species_key(species: str) -> str:
    s = str(species).strip().lower()
    if s == "human_chm13":
        return "human"
    return s


def get_reference_fasta(species: str, data_root: Optional[str] = None) -> Path:
    species_key = normalize_species_key(species)
    if species_key not in REFERENCE_FILES:
        raise ValueError(f"Unknown species key: {species}")

    data_dir = get_data_root(data_root)
    p = data_dir / "reference" / REFERENCE_FILES[species_key]
    if not p.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {p}")
    return p


# ============================================================
# FASTA record utilities
# ============================================================

def parse_chrom_name(header: str) -> str:
    clean = str(header).lstrip(">").strip()
    tok = clean.split()[0]

    if tok.startswith("chr"):
        return tok

    m = re.search(r"chromosome\s+([0-9A-Za-z]+)", clean, re.IGNORECASE)
    if m:
        return f"chr{m.group(1)}"

    if tok in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        return f"chr{tok}"

    return tok


def get_chromosome_sequence(species: str, chrom: str, data_root: Optional[str] = None) -> tuple[str, str]:
    ref_path = get_reference_fasta(species, data_root=data_root)
    chrom = str(chrom)

    exact_match = None
    fallback_match = None

    for rec in SeqIO.parse(str(ref_path), "fasta"):
        rid = parse_chrom_name(rec.description)

        if rid == chrom:
            exact_match = (rid, str(rec.seq).upper())
            break

        if rid.replace("chr", "") == chrom.replace("chr", ""):
            fallback_match = (rid, str(rec.seq).upper())

    if exact_match is not None:
        return exact_match
    if fallback_match is not None:
        return fallback_match

    raise ValueError(f"Chromosome {chrom} not found in {ref_path}")


# ============================================================
# Arm/window helpers
# ============================================================

def is_q_arm(arm: str) -> bool:
    return str(arm).strip().lower().endswith("q")


def arm_to_chrom(arm: str) -> str:
    arm = str(arm).strip()
    if arm.endswith("p") or arm.endswith("q"):
        return arm[:-1]
    return arm


def build_telomere_window_sequence(
    species: str,
    arm: str,
    window_bp: int,
    data_root: Optional[str] = None,
    orient_telomere_left: bool = True,
) -> str:
    chrom = arm_to_chrom(arm)
    _, seq = get_chromosome_sequence(species, chrom, data_root=data_root)

    window_bp = int(window_bp)
    if window_bp <= 0:
        raise ValueError(f"window_bp must be positive: {window_bp}")

    if len(seq) <= window_bp:
        window_seq = seq
    elif is_q_arm(arm):
        window_seq = seq[-window_bp:]
    else:
        window_seq = seq[:window_bp]

    if not orient_telomere_left:
        return window_seq

    # Keep telomere-proximal side on the left in the extracted window.
    # p-arm already has telomere at left; q-arm extracted suffix has telomere at right.
    if is_q_arm(arm):
        return window_seq[::-1]

    return window_seq


# ============================================================
# Signal / downsampling helpers
# ============================================================

def sampled_kmer_repetitiveness_signal(
    seq: str,
    k: int = 20,
    stride_bp: int = 500,
) -> tuple[np.ndarray, np.ndarray]:
    seq = str(seq).upper()
    max_start = len(seq) - k
    if max_start < 0:
        return np.array([], dtype=np.int64), np.array([], dtype=float)

    positions = np.arange(0, max_start + 1, stride_bp, dtype=np.int64)
    kmers = [seq[pos:pos + k] for pos in positions]
    freq = Counter(kmers)
    signal = np.array([freq[km] for km in kmers], dtype=float)
    return positions, signal


def choose_downsample_factor(signal_len: int, target_len: int = 25000) -> int:
    if signal_len <= 0:
        return 1
    return max(1, int(np.ceil(signal_len / float(target_len))))


def downsample_signal_mean(signal: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return signal.copy()

    n = len(signal)
    m = n // factor
    if m == 0:
        return signal.copy()

    trimmed = signal[:m * factor]
    return trimmed.reshape(m, factor).mean(axis=1)


def downsample_positions_mean(positions: np.ndarray, factor: int) -> np.ndarray:
    if factor <= 1:
        return positions.copy()

    n = len(positions)
    m = n // factor
    if m == 0:
        return positions.copy()

    trimmed = positions[:m * factor]
    return trimmed.reshape(m, factor).mean(axis=1).astype(np.int64)


def build_log_scales(scale_min: int = 8, scale_max: int = 1200, n_scales: int = 96) -> np.ndarray:
    raw = np.geomspace(scale_min, scale_max, num=n_scales)
    scales = np.unique(np.round(raw).astype(int))
    scales = scales[scales >= 2]
    return scales


def find_optimal_scale(
    abs_coeffs_matrix: np.ndarray,
    scales_array: np.ndarray,
    peak_prom_frac: float = 0.05,
    width_for_scan: int = 3,
) -> int:
    max_avg_peak = -1.0
    best_scale = None

    for i, sc in enumerate(scales_array):
        row = abs_coeffs_matrix[i, :]
        if row.size == 0:
            continue

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
                best_scale = int(sc)

    if best_scale is None and len(scales_array) > 0:
        best_scale = int(scales_array[len(scales_array) // 2])

    return best_scale


# ============================================================
# Human annotation / RepeatMasker helpers
# ============================================================

def get_centromere_bed(data_root: Optional[str] = None) -> Path:
    data_dir = get_data_root(data_root)
    p = data_dir / "annotation" / "chm13v2.0_merged_centromeres.bed"
    if not p.exists():
        raise FileNotFoundError(f"Centromere BED not found: {p}")
    return p


def read_centromere_interval_from_bed(bed_path: str | Path, chrom: str) -> Optional[tuple[int, int]]:
    chrom = str(chrom)
    bed_path = Path(bed_path)

    with open(bed_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            if parts[0] != chrom:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                return start, end
            except Exception:
                continue

    return None


def get_repeatmasker_out(data_root: Optional[str] = None) -> Path:
    data_dir = get_data_root(data_root)
    p = data_dir / "repeatmasker_out" / "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out"
    if not p.exists():
        raise FileNotFoundError(f"RepeatMasker output not found: {p}")
    return p


def parse_repeatmasker_tar1(path: str | Path, chrom_filter: Optional[str] = None) -> pd.DataFrame:
    rows = []
    chrom_filter = None if chrom_filter is None else str(chrom_filter)

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("SW") or s.startswith("score") or s.startswith("perc"):
                continue

            parts = re.split(r"\s+", s)
            if len(parts) < 14:
                continue

            try:
                query_seq = parts[4]
                q_begin = int(parts[5])
                q_end = int(parts[6])
                strand = parts[8]
                rep_name = parts[9]
                rep_class_family = parts[10]
            except Exception:
                continue

            text = f"{rep_name} {rep_class_family}".upper()
            if "TAR1" not in text:
                continue

            chrom = parse_chrom_name(query_seq)
            if chrom_filter is not None and chrom != chrom_filter:
                continue

            start = min(q_begin, q_end)
            end = max(q_begin, q_end)

            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "rep_name": rep_name,
                "rep_class_family": rep_class_family,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    return df


# ============================================================
# TAR1 block helpers
# ============================================================

def _tar1_blocks_path_for_species(species: str, data_root: Optional[str] = None) -> Path:
    data_dir = get_data_root(data_root)

    mapping = {
        "human": "human_CHM13",
        "chimp": "chimp",
        "bonobo": "bonobo",
        "gorilla": "gorilla",
        "borang": "borang",
        "sorang": "sorang",
    }

    key = normalize_species_key(species)
    if key not in mapping:
        raise ValueError(f"Unknown species for TAR1 blocks: {species}")

    p = data_dir / "tar1_blocks" / mapping[key] / "tar1_blocks.fa"
    if not p.exists():
        raise FileNotFoundError(f"TAR1 blocks FASTA not found: {p}")
    return p


def get_tar1_fasta(species: str, data_root: Optional[str] = None) -> Path:
    return _tar1_blocks_path_for_species(species, data_root=data_root)


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted((int(s), int(e)) for s, e in intervals)
    if not intervals:
        return []

    merged = []
    cur_s, cur_e = intervals[0]

    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e

    merged.append((cur_s, cur_e))
    return merged


def parse_tar1_fasta_headers(fasta_path: Path) -> list[tuple[str, str, int, int]]:
    items = []
    header_re = re.compile(r"^(chr[0-9A-Za-z]+)([pq]):(\d+)-(\d+)$")

    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        rid = rec.id.split()[0]
        m = header_re.match(rid)
        if not m:
            continue

        chrom = m.group(1)
        arm = m.group(2)
        start = int(m.group(3))
        end = int(m.group(4))
        items.append((chrom, arm, min(start, end), max(start, end)))

    return items


def run_counting(fasta_path: Path | str, cutoffs: Iterable[int]) -> tuple[list[int], int]:
    fasta_path = Path(fasta_path)
    items = parse_tar1_fasta_headers(fasta_path)

    counts = []
    for cutoff in cutoffs:
        found_arms = set()
        cutoff = int(cutoff)

        for chrom, arm, s, e in items:
            if s <= cutoff and e <= cutoff:
                found_arms.add(f"{chrom}{arm}")

        counts.append(len(found_arms))

    return counts, len(items)


def tar1_intervals_within_window(
    species: str,
    arm: str,
    window_bp: int,
    data_root: Optional[str] = None,
) -> list[tuple[int, int]]:
    """
    Return TAR1 intervals inside the extracted telomere window, in the same
    coordinate system used by build_telomere_window_sequence(..., orient_telomere_left=True).
    So coordinate 0 is always telomere-proximal.
    """
    fa_path = _tar1_blocks_path_for_species(species, data_root=data_root)
    arm = str(arm)
    window_bp = int(window_bp)

    out = []
    header_re = re.compile(r"^(chr[0-9A-Za-z]+)([pq]):(\d+)-(\d+)$")

    for rec in SeqIO.parse(str(fa_path), "fasta"):
        rid = rec.id.split()[0]
        m = header_re.match(rid)
        if not m:
            continue

        rec_arm = f"{m.group(1)}{m.group(2)}"
        if rec_arm != arm:
            continue

        s = int(m.group(3))
        e = int(m.group(4))
        s, e = min(s, e), max(s, e)

        if s > window_bp:
            continue

        local_s = max(0, s)
        local_e = min(window_bp, e)

        if local_e < local_s:
            continue

        # For q arms, build_telomere_window_sequence(..., orient_telomere_left=True)
        # reverses the extracted suffix. TAR1 coordinates should match that orientation.
        if is_q_arm(arm):
            new_s = max(0, window_bp - local_e)
            new_e = max(0, window_bp - local_s)
            local_s, local_e = min(new_s, new_e), max(new_s, new_e)

        out.append((int(local_s), int(local_e)))

    out.sort(key=lambda x: (x[0], x[1]))
    return out


def whole_chromosome_tar1_intervals(
    species: str,
    chrom: str,
    chrom_len: Optional[int] = None,
    data_root: Optional[str] = None,
) -> list[tuple[int, int, str]]:
    fa_path = _tar1_blocks_path_for_species(species, data_root=data_root)
    chrom = str(chrom)

    p_intervals = []
    q_intervals = []

    header_re = re.compile(r"^(chr[0-9A-Za-z]+)([pq]):(\d+)-(\d+)$")

    for rec in SeqIO.parse(str(fa_path), "fasta"):
        rid = rec.id.split()[0]
        m = header_re.match(rid)
        if not m:
            continue

        rec_chrom = m.group(1)
        arm = m.group(2)
        start = int(m.group(3))
        end = int(m.group(4))

        if rec_chrom != chrom:
            continue

        s = min(start, end)
        e = max(start, end)

        if arm == "p":
            p_intervals.append((s, e))
        else:
            q_intervals.append((s, e))

    merged_p = merge_intervals(p_intervals)
    merged_q = merge_intervals(q_intervals)

    out = []
    for s, e in merged_p:
        out.append((s, e, "p"))

    if chrom_len is not None:
        for s, e in merged_q:
            q_start = max(0, chrom_len - e)
            q_end = max(0, chrom_len - s)
            if q_end < q_start:
                q_start, q_end = q_end, q_start
            out.append((q_start, q_end, "q"))
    else:
        for s, e in merged_q:
            out.append((s, e, "q"))

    out.sort(key=lambda x: (x[0], x[1], x[2]))
    return out


# ============================================================
# Distribution helpers for Figure 2C
# ============================================================

def cumulative_to_bin_counts(x: np.ndarray, cum: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cum = np.asarray(cum, dtype=float)
    x = np.asarray(x, dtype=float)

    if cum.size == 0:
        return x, cum

    n = np.empty_like(cum)
    n[0] = cum[0]
    n[1:] = np.diff(cum)
    n[n < 0] = 0.0
    return x, n


def unify_grid(species_bin_dict: dict[str, tuple[np.ndarray, np.ndarray]]) -> pd.DataFrame:
    all_x = sorted(set(np.concatenate([v[0] for v in species_bin_dict.values()])))
    grid = pd.DataFrame({"cutoff": all_x})

    for sp, (x, nbin) in species_bin_dict.items():
        df_sp = pd.DataFrame({"cutoff": x, sp: nbin})
        grid = grid.merge(df_sp, on="cutoff", how="left")

    return grid.fillna(0.0)


def to_log_coord(x: np.ndarray, x0_ref: float = 1.0) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.maximum(x, 1e-12)
    return np.log(x / x0_ref)