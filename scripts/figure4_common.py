from __future__ import annotations

import os
import re
import glob
from pathlib import Path
from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy import stats


# =========================================================
# Basic paths
# =========================================================

def get_repo_root() -> Path:
    here = Path(__file__).resolve()
    for p in [here.parent] + list(here.parents):
        if (p / "scripts").exists() and (p / "data").exists():
            return p
    return here.parent.parent


def get_data_root(data_root: Optional[str] = None) -> Path:
    if data_root is not None:
        return Path(data_root).resolve()
    env = os.environ.get("TAR1_DATA_ROOT")
    if env:
        return Path(env).resolve()
    return (get_repo_root() / "data").resolve()


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# =========================================================
# Canonical constants from original dump
# =========================================================

QUERY_SEQ = "GGCGCGCCGCGCCGGCGCAGGCGCAGAGA"
SCORE_THRESHOLD = 20

MIN_POS = 500
MAX_POS = 30000
HARD_CUTOFF = 2500

AGE_GROUP_ORDER = ["Age 0", "Age 1-59", "Age 60+"]
GROUP_ORDER_RATIO = ["Newborn (Age 0)", "Adult (Age 1+)", "Cancer"]
FINAL_GROUP_ORDER = [
    "Cancer (Fast Growth)",
    "Cancer (Slow Growth)",
    "Short Telomere Disease",
    "Age 0",
    "Age 1-59",
    "Age 60+",
]
SHAPE_ORDER = ["TERT", "ALT", "Normal", "Unknown"]

COLOR_DICT = {
    "Cancer (Fast Growth)": "#9ecae1",
    "Cancer (Slow Growth)": "#3182bd",
    "Short Telomere Disease": "green",
    "Age 0": "#cb181d",
    "Age 1-59": "#fb6a4a",
    "Age 60+": "#fcae91",
}

MARKER_DICT = {
    "TERT": "s",
    "ALT": "^",
    "Normal": "o",
    "Unknown": "X",
}

DT_CUTOFF = 85.0

CELL_LINE_TMM = {
    "U2-OS": "ALT",
    "Saos-2": "ALT",
    "G-292": "ALT",
    "SK-N-F1": "ALT",
    "SK-N-AS": "ALT",
    "HT-1080": "TERT",
    "HT-29": "TERT",
    "HOS": "TERT",
    "Calu-3": "TERT",
    "SK-LU-1": "TERT",
}


# =========================================================
# File resolution
# =========================================================

def get_tar1_blocks_root(data_root: Optional[str] = None) -> Path:
    return get_data_root(data_root) / "tar1_blocks"


def get_cluster_master_csv(data_root: Optional[str] = None) -> Path:
    p = get_tar1_blocks_root(data_root) / "TAR1_clustering_results_master_BLAST_method_29bp_2clusters.csv"
    return p


def get_position_stats_csv(data_root: Optional[str] = None) -> Path:
    p = get_tar1_blocks_root(data_root) / "TAR1_MultiCluster_Position_Stats.csv"
    return p


def get_srr_only_median_csv(data_root: Optional[str] = None) -> Path:
    p = get_tar1_blocks_root(data_root) / "TAR1_SRR_Only_Median_Positions.csv"
    return p


def get_donor_meta_csv(data_root: Optional[str] = None) -> Path:
    p = get_data_root(data_root) / "csv" / "table_s2_donor_telomere_lengths.csv"
    if not p.exists():
        raise FileNotFoundError(f"Donor metadata CSV not found: {p}")
    return p


def get_schmidt_csv(data_root: Optional[str] = None) -> Path:
    p = get_data_root(data_root) / "csv" / "schmidt_corres.csv"
    if not p.exists():
        raise FileNotFoundError(f"Schmidt correspondence CSV not found: {p}")
    return p


def get_delta_a_csv(data_root: Optional[str] = None) -> Path:
    p = get_data_root(data_root) / "csv" / "delta_a_summary_SRR_only.csv"
    if not p.exists():
        raise FileNotFoundError(f"delta_a CSV not found: {p}")
    return p


# =========================================================
# Shared helpers
# =========================================================

def normalize_key(name: str) -> str:
    if pd.isna(name):
        return ""
    s = str(name).strip().replace("_", ".")
    return re.sub(r"\.+", ".", s)


def donor_id_from_sample_id(sample_id: str) -> str:
    s = str(sample_id).strip()
    if "." in s:
        return s.split(".", 1)[0]
    return s


def assign_group(sample_id: str) -> str:
    s = str(sample_id)
    if s.startswith("SRR"):
        return "Cancer"
    elif s.startswith("JH") or "nanopore" in s:
        return "Normal"
    return "Unknown"


def parse_header_interval_from_string(header: str) -> tuple[Optional[int], Optional[int]]:
    if pd.isna(header):
        return None, None
    m = re.search(r"(\d+)-(\d+)\s*$", str(header).strip())
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def parse_header_interval_from_fasta_line(line: str) -> tuple[Optional[int], Optional[int]]:
    line = line.strip()
    if not line.startswith(">") or ":" not in line:
        return None, None
    _, rest = line[1:].split(":", 1)
    m = re.search(r"(\d+)-(\d+)", rest)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def calculate_weighted_median(x_positions: np.ndarray, y_counts: np.ndarray, cutoff: int = HARD_CUTOFF) -> tuple[float, int]:
    mask = x_positions >= cutoff
    x_filt = x_positions[mask]
    y_filt = y_counts[mask]

    total_count = int(y_filt.sum())
    if total_count == 0:
        return np.nan, 0

    cumsum = np.cumsum(y_filt)
    idx = int(np.searchsorted(cumsum, total_count / 2.0))
    if idx >= len(x_filt):
        idx = len(x_filt) - 1
    return float(x_filt[idx]), total_count


def get_star_from_pvalue(p: float) -> str:
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# =========================================================
# Step 1: exact 29bp local-alignment clustering
# =========================================================

def build_cluster_master_from_tar1_blocks(data_root: Optional[str] = None, save_csv: bool = True) -> pd.DataFrame:
    base_dir = get_tar1_blocks_root(data_root)
    target_files = sorted(glob.glob(str(base_dir / "*" / "tar1_blocks.fa")))

    if not target_files:
        raise FileNotFoundError(f"No tar1_blocks.fa files found under: {base_dir}")

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -1.0

    all_results = []

    for fasta_path in target_files:
        sample_id = os.path.basename(os.path.dirname(fasta_path))
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            continue

        for rec in records:
            seq_str = str(rec.seq).upper()
            score = float(aligner.score(seq_str, QUERY_SEQ))

            if score >= SCORE_THRESHOLD:
                cluster_code = 0
                cluster_name = "Cluster_0 (Contains 29bp)"
            else:
                cluster_code = 1
                cluster_name = "Cluster_1 (No 29bp)"

            original_read_id = rec.id.split(":")[0] if ":" in rec.id else rec.id

            all_results.append({
                "Sample_ID": sample_id,
                "Group": assign_group(sample_id),
                "Full_Header": rec.id,
                "Read_ID": original_read_id,
                "Sequence": seq_str,
                "Length": len(seq_str),
                "Alignment_Score": score,
                "Cluster": cluster_name,
                "Cluster_Code": cluster_code,
            })

    if not all_results:
        raise RuntimeError("No TAR1 reads processed while building cluster master.")

    df = pd.DataFrame(all_results)

    if save_csv:
        out_csv = get_cluster_master_csv(data_root)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)

    return df


# =========================================================
# Step 2: exact position stats from cluster master + fasta headers
# =========================================================

def build_position_stats_from_cluster_master(data_root: Optional[str] = None, save_csv: bool = True) -> pd.DataFrame:
    cluster_csv = get_cluster_master_csv(data_root)
    if not cluster_csv.exists():
        df_meta = build_cluster_master_from_tar1_blocks(data_root=data_root, save_csv=True)
    else:
        df_meta = pd.read_csv(cluster_csv)

    unique_clusters = sorted(df_meta["Cluster_Code"].dropna().astype(int).unique().tolist())
    cluster_map: dict[str, dict[str, int]] = {}

    for _, row in df_meta.iterrows():
        sid = row["Sample_ID"]
        header = row["Full_Header"]
        cluster = int(row["Cluster_Code"])
        if sid not in cluster_map:
            cluster_map[sid] = {}
        cluster_map[sid][header] = cluster

    base_dir = get_tar1_blocks_root(data_root)
    samples = sorted([d for d in os.listdir(base_dir) if (base_dir / d).is_dir()])

    results = []
    x_positions = np.arange(MIN_POS, MAX_POS + 1)

    for entry in samples:
        fasta_path = base_dir / entry / "tar1_blocks.fa"
        if not fasta_path.exists():
            continue
        if entry not in cluster_map:
            continue

        sample_header_map = cluster_map[entry]
        sample_key = normalize_key(entry)

        cluster_covs = {c: np.zeros(MAX_POS + 2, dtype=int) for c in unique_clusters}

        with open(fasta_path, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    continue

                full_header = line.strip()[1:]
                c_code = sample_header_map.get(full_header)
                if c_code is None or c_code not in cluster_covs:
                    continue

                s_raw, e_raw = parse_header_interval_from_fasta_line(line)
                if s_raw is None:
                    continue

                if sample_key.startswith("JH") and (s_raw > 15000 or e_raw > 15000):
                    continue

                s = max(min(s_raw, e_raw), MIN_POS)
                e = min(max(s_raw, e_raw), MAX_POS)
                if s <= e:
                    cluster_covs[c_code][s:e + 1] += 1

        row_data = {
            "Sample_ID": entry,
            "Group": "Normal" if (entry.startswith("JH") or "nanopore" in entry) else "Cancer"
        }

        for c_code in unique_clusters:
            cov_slice = cluster_covs[c_code][MIN_POS:MAX_POS + 1]
            med, count = calculate_weighted_median(x_positions, cov_slice, cutoff=HARD_CUTOFF)
            row_data[f"Cluster_{c_code}_Pos_Median"] = med
            row_data[f"Cluster_{c_code}_Valid_Bases"] = int(count)

        results.append(row_data)

    df_res = pd.DataFrame(results)

    if save_csv:
        out_csv = get_position_stats_csv(data_root)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df_res.to_csv(out_csv, index=False)

        srr_df = df_res[df_res["Sample_ID"].astype(str).str.startswith("SRR")].copy()
        if not srr_df.empty:
            median_cols = sorted([c for c in srr_df.columns if "Pos_Median" in c])
            srr_report = srr_df[["Sample_ID"] + median_cols].copy()
            rename_map = {c: c.replace("_Pos_Median", "") for c in median_cols}
            srr_report = srr_report.rename(columns=rename_map).sort_values("Sample_ID")
            srr_report.to_csv(get_srr_only_median_csv(data_root), index=False)

    return df_res


# =========================================================
# Metadata loaders
# =========================================================

def load_donor_meta(data_root: Optional[str] = None) -> pd.DataFrame:
    df = pd.read_csv(get_donor_meta_csv(data_root))
    df["donor_id"] = df["donor_id"].astype(str).str.strip()
    df["age_years"] = pd.to_numeric(df["age_years"], errors="coerce")
    if "ipf_or_short_telomere_disorder" in df.columns:
        tmp = df["ipf_or_short_telomere_disorder"]
        if tmp.dtype != bool:
            tmp = tmp.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y", "t"])
        df["ipf_or_short_telomere_disorder"] = tmp.astype(bool)
    else:
        df["ipf_or_short_telomere_disorder"] = False
    return df


def load_schmidt_meta(data_root: Optional[str] = None) -> pd.DataFrame:
    df = pd.read_csv(get_schmidt_csv(data_root))
    if "Sample_Name" in df.columns:
        df["TMM_Status"] = df["Sample_Name"].map(CELL_LINE_TMM).fillna("Unknown")
    else:
        df["TMM_Status"] = "Unknown"
    return df


# =========================================================
# Final integrated position table for B-F and K
# =========================================================

def build_final_position_table(data_root: Optional[str] = None, save_intermediate: bool = True) -> pd.DataFrame:
    pos_csv = get_position_stats_csv(data_root)
    if not pos_csv.exists():
        df_pos = build_position_stats_from_cluster_master(data_root=data_root, save_csv=True)
    else:
        df_pos = pd.read_csv(pos_csv)

    donor_meta = load_donor_meta(data_root)
    schmidt = load_schmidt_meta(data_root)

    donor_age_map = dict(zip(donor_meta["donor_id"], donor_meta["age_years"]))
    donor_disease_map = dict(zip(donor_meta["donor_id"], donor_meta["ipf_or_short_telomere_disorder"]))

    dt_df = pd.DataFrame({
        "Sample_Name": ["Calu-3", "G-292", "HOS", "HT-1080", "HT-29", "SK-LU-1", "SK-N-AS", "Saos-2", "U2-OS", "SK-N-F1"],
        "Doubling_Time": [105.83, 102.19, 41.06, 78.31, 28.35, 98.21, 46.57, 99.13, 71.88, 97.22],
    })

    srr_dt_map = {}
    srr_tmm_map = {}
    if {"Sample_Name", "Run_ID"}.issubset(schmidt.columns):
        merged = pd.merge(schmidt, dt_df, on="Sample_Name", how="inner")
        srr_dt_map = dict(zip(merged["Run_ID"], merged["Doubling_Time"]))
        srr_tmm_map = dict(zip(merged["Run_ID"], merged["TMM_Status"]))

    df = df_pos.copy()
    df["Sample_ID"] = df["Sample_ID"].astype(str)
    df["Donor_ID"] = df["Sample_ID"].apply(donor_id_from_sample_id)
    df["Age"] = df["Donor_ID"].map(donor_age_map)
    df["Is_Disease"] = df["Donor_ID"].map(donor_disease_map).fillna(False).astype(bool)
    df["Doubling_Time"] = df["Sample_ID"].map(srr_dt_map)
    df["TMM_Status"] = df["Sample_ID"].map(srr_tmm_map).fillna("Unknown")

    df = df.dropna(subset=["Cluster_0_Pos_Median", "Cluster_1_Pos_Median"], how="all").copy()

    def assign_group_and_shape(row: pd.Series) -> pd.Series:
        sid = str(row["Sample_ID"])

        if sid.startswith("SRR"):
            dt = row["Doubling_Time"]
            shape_type = row["TMM_Status"] if row["TMM_Status"] in ["TERT", "ALT"] else "Unknown"

            if pd.isna(dt):
                return pd.Series({"Final_Group": "Exclude", "Shape_Type": shape_type})

            if float(dt) < DT_CUTOFF:
                fg = "Cancer (Fast Growth)"
            else:
                fg = "Cancer (Slow Growth)"
            return pd.Series({"Final_Group": fg, "Shape_Type": shape_type})

        if bool(row["Is_Disease"]):
            return pd.Series({"Final_Group": "Short Telomere Disease", "Shape_Type": "Normal"})

        age = row["Age"]
        if pd.isna(age):
            return pd.Series({"Final_Group": "Exclude", "Shape_Type": "Normal"})

        age = int(age)
        if age == 0:
            fg = "Age 0"
        elif 1 <= age <= 59:
            fg = "Age 1-59"
        elif age >= 60:
            fg = "Age 60+"
        else:
            fg = "Exclude"

        return pd.Series({"Final_Group": fg, "Shape_Type": "Normal"})

    info = df.apply(assign_group_and_shape, axis=1)
    df["Final_Group"] = info["Final_Group"]
    df["Shape_Type"] = info["Shape_Type"]
    df = df[df["Final_Group"] != "Exclude"].copy()

    def map_ratio_group(g: str) -> str:
        if g == "Age 0":
            return "Newborn (Age 0)"
        if g in ["Age 1-59", "Age 60+"]:
            return "Adult (Age 1+)"
        if "Cancer" in g:
            return "Cancer"
        return "Exclude"

    df["Ratio_Group"] = df["Final_Group"].apply(map_ratio_group)
    df["Pos_Ratio"] = df["Cluster_1_Pos_Median"] / df["Cluster_0_Pos_Median"]

    return df


# =========================================================
# Figure 4A-like overall TAR1 median vs telomere length
# =========================================================

def load_jh_ground_truth_map(data_root: Optional[str] = None) -> dict[str, float]:
    meta = load_donor_meta(data_root)
    if "fastq_file_name" not in meta.columns or "mean_telomere_length_bps" not in meta.columns:
        return {}
    meta["match_key"] = meta["fastq_file_name"].apply(normalize_key)
    return dict(zip(meta["match_key"], meta["mean_telomere_length_bps"]))


def load_srr_ground_truth_map(data_root: Optional[str] = None) -> dict[str, float]:
    schmidt = load_schmidt_meta(data_root)
    if {"Run_ID", "Mean_bp"}.issubset(schmidt.columns):
        valid = schmidt.dropna(subset=["Run_ID", "Mean_bp"])
        return dict(zip(valid["Run_ID"], valid["Mean_bp"]))
    return {}


def load_all_coverages(
    root_dir: str | Path,
    min_pos: int = MIN_POS,
    max_pos: int = MAX_POS,
    apply_jh_15kb_filter: bool = True,
) -> dict[str, np.ndarray]:
    root_dir = str(root_dir)
    coverage_data = {}
    samples = sorted([d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))])

    for entry in samples:
        sample_dir = os.path.join(root_dir, entry)
        fasta_path = os.path.join(sample_dir, "tar1_blocks.fa")
        if not os.path.exists(fasta_path):
            continue

        sample_key = normalize_key(entry)
        cov_array = np.zeros(max_pos + 2, dtype=int)

        with open(fasta_path, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                s_raw, e_raw = parse_header_interval_from_fasta_line(line)
                if s_raw is None:
                    continue

                if apply_jh_15kb_filter and sample_key.startswith("JH") and (s_raw > 15000 or e_raw > 15000):
                    continue

                s = max(min(s_raw, e_raw), min_pos)
                e = min(max(s_raw, e_raw), max_pos)
                if s <= e:
                    cov_array[s:e + 1] += 1

        coverage_data[sample_key] = cov_array[min_pos:max_pos + 1]

    return coverage_data


def calculate_tar1_stats(x_positions: np.ndarray, y_counts: np.ndarray, cutoff: int) -> tuple[float, float, float, int]:
    mask = x_positions >= cutoff
    x_filt = x_positions[mask]
    y_filt = y_counts[mask]
    total_count = int(y_filt.sum())

    if total_count == 0:
        return np.nan, np.nan, np.nan, 0

    cumsum = np.cumsum(y_filt)
    idx = int(np.searchsorted(cumsum, total_count / 2.0))
    if idx >= len(x_filt):
        idx = len(x_filt) - 1
    median_val = float(x_filt[idx])

    mean_val = float(np.average(x_filt, weights=y_filt))
    peak_val = float(x_filt[int(np.argmax(y_filt))])
    return median_val, mean_val, peak_val, total_count


def build_total_human_tar1_vs_telomere_table(data_root: Optional[str] = None) -> pd.DataFrame:
    tar1_root = get_tar1_blocks_root(data_root)
    x_positions = np.arange(MIN_POS, MAX_POS + 1)

    coverage_dict = load_all_coverages(
        root_dir=tar1_root,
        min_pos=MIN_POS,
        max_pos=MAX_POS,
        apply_jh_15kb_filter=True,
    )

    gt_jh = load_jh_ground_truth_map(data_root)
    gt_srr = load_srr_ground_truth_map(data_root)

    rows = []
    for sample_key, cov in coverage_dict.items():
        if sample_key.startswith("SRR"):
            gt = gt_srr.get(sample_key)
        else:
            gt = gt_jh.get(sample_key)

        med, mean_val, peak, count_valid = calculate_tar1_stats(x_positions, cov, cutoff=HARD_CUTOFF)
        rows.append({
            "sample_key": sample_key,
            "sample_group": "Cancer" if sample_key.startswith("SRR") else "Non-cancer",
            "ground_truth_bps": gt,
            "tar1_median_bps": med,
            "tar1_mean_bps": mean_val,
            "tar1_peak_bps": peak,
            "read_count_valid": count_valid,
        })

    return pd.DataFrame(rows)