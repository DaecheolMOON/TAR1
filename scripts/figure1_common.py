from __future__ import annotations

import os
import re
from pathlib import Path
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

DATA_REGISTRY: Dict[str, Dict[str, Optional[str]]] = {
    "human": {"tar1_dir": os.path.join("tar1_blocks", "human_CHM13")},
    "chimp": {"tar1_dir": os.path.join("tar1_blocks", "chimp")},
    "bonobo": {"tar1_dir": os.path.join("tar1_blocks", "bonobo")},
    "gorilla": {"tar1_dir": os.path.join("tar1_blocks", "gorilla")},
    "borang": {"tar1_dir": os.path.join("tar1_blocks", "borang")},
    "sorang": {"tar1_dir": os.path.join("tar1_blocks", "sorang")},
    "siamang": {"tar1_dir": None},
}

MONOMERS_P: Dict[str, Tuple[str, ...]] = {
    "61bp": ("gcacgcccgcctgctggcagctggggacactgctgggccctcttgctccaacagtagtggc",),
    "29bp": ("ggcgcgccgcgccggcgcaggcgcagaga",),
    "37bp": ("ccgctctgtgctgacgagaacgcaactcggccgtcgc",),
    "40bp": ("agcgcgggggtggcgcggtgcacgcgcagagacacacgtc",),
}

MONOMERS_Q: Dict[str, Tuple[str, ...]] = {
    "61bp": ("gccactactgttggagcaagagggcccagcagtgtccccagctgccagcaggcgggcgtgc",),
    "29bp": ("tctctgcgcctgcgccggcgcggcgcgcc",),
    "37bp": ("gcgacggccgagttgcgttctcgtcagcacagagcgg",),
    "40bp": ("gacgtgtgtctctgcgcgtgcaccgcgccacccccgcgct",),
}


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sanitize_id_for_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def is_q_arm(arm_name: str) -> bool:
    return str(arm_name).strip().lower().endswith("q")


def get_analysis_sequence(seq_full: str, seq_id: str, length: int) -> str:
    seq_full = str(seq_full).upper()
    if len(seq_full) <= length:
        return seq_full
    return seq_full[-length:] if is_q_arm(seq_id) else seq_full[:length]


def parse_fasta_simple(filename: str | Path) -> List[Tuple[str, str]]:
    sequences: List[Tuple[str, str]] = []
    current_header: Optional[str] = None
    current_seq: List[str] = []
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq).upper()))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        sequences.append((current_header, "".join(current_seq).upper()))
    return sequences


def write_single_fasta(header: str, sequence: str, out_path: str | Path, wrap: int = 60) -> None:
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), wrap):
            f.write(sequence[i:i + wrap] + "\n")


def read_single_fasta(fasta_path: str | Path) -> Tuple[str, str]:
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected exactly 1 FASTA record in {fasta_path}, found {len(records)}")
    rec = records[0]
    return rec.id, str(rec.seq).upper()


def resolve_tar1_blocks_fasta(species: str, data_root: str | Path) -> Path:
    if species not in DATA_REGISTRY:
        raise KeyError(f"Unknown species: {species}")
    tar1_dir = DATA_REGISTRY[species].get("tar1_dir")
    if not tar1_dir:
        raise FileNotFoundError(f"No tar1_dir configured for species={species}")
    fasta = Path(data_root) / tar1_dir / "tar1_blocks.fa"
    if not fasta.exists():
        raise FileNotFoundError(f"TAR1 blocks FASTA not found: {fasta}")
    return fasta


def load_manifest(manifest_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(manifest_path, sep="\t")
    required = ["panel", "species", "record_id", "match_mode", "use_for_cwt", "use_for_motif"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Manifest missing required columns: {missing}")

    if "min_length" not in df.columns:
        df["min_length"] = 0
    if "selection_mode" not in df.columns:
        df["selection_mode"] = "auto"
    if "merge_gap" not in df.columns:
        df["merge_gap"] = 500

    df["min_length"] = pd.to_numeric(df["min_length"], errors="coerce").fillna(0).astype(int)
    df["merge_gap"] = pd.to_numeric(df["merge_gap"], errors="coerce").fillna(500).astype(int)
    df["selection_mode"] = df["selection_mode"].fillna("auto").astype(str)
    return df


def match_header(header: str, record_id: str, match_mode: str = "exact") -> bool:
    header_core = header.split()[0]

    if match_mode == "exact":
        return header_core == record_id

    if match_mode == "contains":
        return (record_id in header) or (record_id in header_core)

    if match_mode == "startswith":
        return header_core.startswith(record_id)

    raise ValueError(f"Unsupported match_mode: {match_mode}")


def extract_arm_from_text(text: str) -> Optional[str]:
    m = re.search(r"(chr[0-9XYM]+[pq])", str(text), flags=re.IGNORECASE)
    if m:
        return m.group(1).lower()
    return None


def parse_coordinates_from_header(header: str) -> Tuple[Optional[int], Optional[int]]:
    text = str(header)

    patterns = [
        r":(\d+)-(\d+)",
        r"start=(\d+).*?end=(\d+)",
        r"from=(\d+).*?to=(\d+)",
        r"\b(\d+)-(\d+)\b",
    ]

    for pat in patterns:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            a = int(m.group(1))
            b = int(m.group(2))
            if a <= b:
                return a, b
            return b, a

    nums = [int(x) for x in re.findall(r"\d+", text)]
    for i in range(len(nums) - 1):
        a, b = nums[i], nums[i + 1]
        if a < b:
            return a, b

    return None, None


def build_record_info(header: str, seq: str) -> Dict[str, object]:
    header_core = header.split()[0]
    arm = extract_arm_from_text(header_core) or extract_arm_from_text(header)
    start, end = parse_coordinates_from_header(header)

    return {
        "header": header,
        "header_core": header_core,
        "seq": seq,
        "length": len(seq),
        "arm": arm,
        "start": start,
        "end": end,
    }


def collect_matching_records(
    fasta_path: str | Path,
    record_id: str,
    match_mode: str = "exact",
) -> List[Dict[str, object]]:
    matches: List[Dict[str, object]] = []
    for header, seq in parse_fasta_simple(fasta_path):
        if match_header(header, record_id=record_id, match_mode=match_mode):
            matches.append(build_record_info(header=header, seq=seq))
    return matches


def telomere_priority_tuple(
    record: Dict[str, object],
    requested_arm: Optional[str] = None,
) -> Tuple[float, float]:
    arm = requested_arm or record.get("arm")
    start = record.get("start")
    end = record.get("end")
    length = int(record["length"])

    if isinstance(arm, str) and arm.lower().endswith("p"):
        if start is not None:
            return (float(start), -float(length))
        return (float("inf"), -float(length))

    if isinstance(arm, str) and arm.lower().endswith("q"):
        if end is not None:
            return (-float(end), -float(length))
        return (float("inf"), -float(length))

    return (float("inf"), -float(length))


def filter_by_requested_arm(
    matches: Sequence[Dict[str, object]],
    requested_record_id: str,
) -> List[Dict[str, object]]:
    requested_arm = extract_arm_from_text(requested_record_id)
    if requested_arm is None:
        return list(matches)

    subset = [m for m in matches if str(m.get("arm") or "").lower() == requested_arm.lower()]
    return subset if subset else list(matches)


def merge_records_by_gap(
    records: Sequence[Dict[str, object]],
    requested_record_id: str,
    merge_gap: int = 500,
) -> List[Dict[str, object]]:
    if not records:
        return []

    requested_arm = extract_arm_from_text(requested_record_id)
    arm_filtered = filter_by_requested_arm(records, requested_record_id=requested_record_id)

    sortable = []
    unsortable = []
    for rec in arm_filtered:
        if rec.get("start") is None or rec.get("end") is None:
            unsortable.append(rec)
        else:
            sortable.append(rec)

    sortable = sorted(sortable, key=lambda r: (int(r["start"]), int(r["end"])))

    merged: List[Dict[str, object]] = []
    current_cluster: List[Dict[str, object]] = []

    for rec in sortable:
        if not current_cluster:
            current_cluster = [rec]
            continue

        prev = current_cluster[-1]
        gap = int(rec["start"]) - int(prev["end"])
        if gap <= int(merge_gap):
            current_cluster.append(rec)
        else:
            merged.append(collapse_record_cluster(current_cluster, requested_arm=requested_arm))
            current_cluster = [rec]

    if current_cluster:
        merged.append(collapse_record_cluster(current_cluster, requested_arm=requested_arm))

    for rec in unsortable:
        merged.append(collapse_record_cluster([rec], requested_arm=requested_arm))

    return merged


def collapse_record_cluster(
    cluster: Sequence[Dict[str, object]],
    requested_arm: Optional[str] = None,
) -> Dict[str, object]:
    if not cluster:
        raise ValueError("collapse_record_cluster() received an empty cluster.")

    ordered = list(cluster)
    ordered = sorted(
        ordered,
        key=lambda r: (
            float("inf") if r.get("start") is None else int(r["start"]),
            float("inf") if r.get("end") is None else int(r["end"]),
        ),
    )

    arm = requested_arm or ordered[0].get("arm")
    starts = [int(r["start"]) for r in ordered if r.get("start") is not None]
    ends = [int(r["end"]) for r in ordered if r.get("end") is not None]

    start = min(starts) if starts else None
    end = max(ends) if ends else None
    seq = "".join(str(r["seq"]) for r in ordered)
    total_length = len(seq)

    if start is not None and end is not None:
        header = f"{arm}:{start}-{end}|merged_n={len(ordered)}|seq_len={total_length}"
        header_core = f"{arm}:{start}-{end}"
    else:
        header = f"{arm or 'unknown'}|merged_n={len(ordered)}|seq_len={total_length}"
        header_core = f"{arm or 'unknown'}|merged"

    return {
        "header": header,
        "header_core": header_core,
        "seq": seq,
        "length": total_length,
        "arm": arm,
        "start": start,
        "end": end,
        "n_merged": len(ordered),
        "source_headers": [str(r["header"]) for r in ordered],
    }


def choose_best_record_or_cluster(
    matches: Sequence[Dict[str, object]],
    requested_record_id: str,
    min_length: int = 0,
    merge_gap: int = 500,
) -> Dict[str, object]:
    if not matches:
        raise ValueError("No matched records were provided.")

    requested_arm = extract_arm_from_text(requested_record_id)
    merged_candidates = merge_records_by_gap(
        records=matches,
        requested_record_id=requested_record_id,
        merge_gap=merge_gap,
    )

    passing = [m for m in merged_candidates if int(m["length"]) >= int(min_length)]
    candidate_pool = passing if passing else merged_candidates

    if not candidate_pool:
        raise ValueError("No candidates remained after merging/filtering.")

    ranked = sorted(
        candidate_pool,
        key=lambda rec: telomere_priority_tuple(rec, requested_arm=requested_arm),
    )
    return ranked[0]


def find_record_in_fasta(
    fasta_path: str | Path,
    record_id: str,
    match_mode: str = "exact",
    min_length: int = 0,
    selection_mode: str = "auto",
    merge_gap: int = 500,
) -> Tuple[str, str]:
    matches = collect_matching_records(
        fasta_path=fasta_path,
        record_id=record_id,
        match_mode=match_mode,
    )

    if not matches:
        raise ValueError(
            f"No record matched record_id={record_id!r} "
            f"with match_mode={match_mode!r} in {fasta_path}"
        )

    if len(matches) == 1 and int(matches[0]["length"]) >= int(min_length):
        return str(matches[0]["header"]), str(matches[0]["seq"])

    if selection_mode == "strict":
        raise ValueError(
            f"Multiple records matched record_id={record_id!r} with match_mode={match_mode!r} "
            f"in {fasta_path}. Please make the manifest more specific."
        )

    best = choose_best_record_or_cluster(
        matches=matches,
        requested_record_id=record_id,
        min_length=min_length,
        merge_gap=merge_gap,
    )
    return str(best["header"]), str(best["seq"])


def summarize_matches(
    matches: Sequence[Dict[str, object]],
    requested_record_id: str,
) -> pd.DataFrame:
    requested_arm = extract_arm_from_text(requested_record_id)
    rows = []
    for rec in matches:
        tel_key = telomere_priority_tuple(rec, requested_arm=requested_arm)
        rows.append({
            "header": rec["header"],
            "arm": rec["arm"],
            "start": rec["start"],
            "end": rec["end"],
            "length": rec["length"],
            "telomere_rank_key_1": tel_key[0],
            "telomere_rank_key_2": tel_key[1],
        })
    return pd.DataFrame(rows)


def summarize_merged_candidates(
    matches: Sequence[Dict[str, object]],
    requested_record_id: str,
    merge_gap: int = 500,
) -> pd.DataFrame:
    requested_arm = extract_arm_from_text(requested_record_id)
    merged = merge_records_by_gap(
        records=matches,
        requested_record_id=requested_record_id,
        merge_gap=merge_gap,
    )

    rows = []
    for rec in merged:
        tel_key = telomere_priority_tuple(rec, requested_arm=requested_arm)
        rows.append({
            "header": rec["header"],
            "arm": rec["arm"],
            "start": rec["start"],
            "end": rec["end"],
            "length": rec["length"],
            "n_merged": rec.get("n_merged", 1),
            "source_headers": " || ".join(rec.get("source_headers", [])),
            "telomere_rank_key_1": tel_key[0],
            "telomere_rank_key_2": tel_key[1],
        })
    return pd.DataFrame(rows)


def kmer_repetitiveness_signal(seq: str, k: int = 20) -> np.ndarray:
    seq = seq.upper()
    if len(seq) < k:
        return np.array([], dtype=float)
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    freq = Counter(kmers)
    return np.array([freq[km] for km in kmers], dtype=float)