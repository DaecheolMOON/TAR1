from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure1_common import (
    MONOMERS_P,
    MONOMERS_Q,
    ensure_dir,
    is_q_arm,
    load_manifest,
    read_single_fasta,
    sanitize_id_for_filename,
)


def edit_distance(s1: str, s2: str) -> int:
    s1 = str(s1)
    s2 = str(s2)

    if s1 == s2:
        return 0
    if len(s1) == 0:
        return len(s2)
    if len(s2) == 0:
        return len(s1)

    if len(s1) < len(s2):
        s1, s2 = s2, s1

    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1, start=1):
        current = [i]
        for j, c2 in enumerate(s2, start=1):
            insert_cost = current[j - 1] + 1
            delete_cost = previous[j] + 1
            sub_cost = previous[j - 1] + (0 if c1 == c2 else 1)
            current.append(min(insert_cost, delete_cost, sub_cost))
        previous = current
    return previous[-1]


def get_monomer_dict_from_header_or_id(header: str) -> Dict[str, Tuple[str, ...]]:
    return MONOMERS_Q if is_q_arm(header) else MONOMERS_P


def build_fasta_path(fasta_dir: str | Path, panel: str, species: str, record_id: str) -> Path:
    fname = f"{panel}_{species}_{sanitize_id_for_filename(record_id)}.fa"
    return Path(fasta_dir) / fname


def compute_match_count_profile(seq: str, motif_seq: str) -> np.ndarray:
    seq = seq.lower()
    motif_seq = motif_seq.lower()
    motif_len = len(motif_seq)

    n = len(seq) - motif_len + 1
    if n <= 0:
        return np.array([], dtype=int)

    match_counts = np.empty(n, dtype=int)
    for i in range(n):
        window = seq[i:i + motif_len]
        dist = edit_distance(window, motif_seq)
        match_counts[i] = motif_len - dist
    return match_counts


def detect_threshold_peaks(match_counts: np.ndarray, min_score_frac: float) -> Tuple[np.ndarray, int]:
    if match_counts.size == 0:
        return np.array([], dtype=int), 0

    global_peak = int(match_counts.max())
    threshold = int(global_peak * float(min_score_frac))

    peaks: List[int] = []
    if match_counts.size >= 3:
        for i in range(1, len(match_counts) - 1):
            if (
                match_counts[i] > match_counts[i - 1]
                and match_counts[i] > match_counts[i + 1]
                and match_counts[i] >= threshold
            ):
                peaks.append(i)

    return np.array(peaks, dtype=int), threshold


def summarize_profile(
    panel: str,
    species: str,
    record_id: str,
    header: str,
    motif_name: str,
    motif_len: int,
    match_counts: np.ndarray,
    peaks: np.ndarray,
    threshold: int,
) -> Dict[str, object]:
    return {
        "panel": panel,
        "species": species,
        "record_id": record_id,
        "header": header,
        "motif_name": motif_name,
        "motif_len": motif_len,
        "profile_length": int(len(match_counts)),
        "global_peak": int(match_counts.max()) if match_counts.size else 0,
        "threshold": int(threshold),
        "n_threshold_peaks": int(len(peaks)),
    }


def plot_match_count_profile(
    match_counts: np.ndarray,
    peaks: np.ndarray,
    motif_name: str,
    motif_len: int,
    out_png: str | Path,
    title_prefix: str | None = None,
) -> None:
    plt.figure(figsize=(8, 3))
    plt.plot(match_counts, label="match count")

    if len(peaks) > 0:
        plt.scatter(
            peaks,
            match_counts[peaks],
            s=20,
            label="≥threshold peaks",
        )

    if title_prefix:
        plt.title(f"{title_prefix} | {motif_name.replace('bp', ' bp')} (len={motif_len})", fontsize=9)
    else:
        plt.title(f"{motif_name.replace('bp', ' bp')} (len={motif_len})", fontsize=9)

    plt.xlabel("position")
    plt.ylabel(f"matches / {motif_len}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def save_profile_csv(
    out_csv: str | Path,
    match_counts: np.ndarray,
    peaks: np.ndarray,
    threshold: int,
) -> None:
    if match_counts.size == 0:
        df = pd.DataFrame(columns=["position", "match_count", "is_threshold_peak", "threshold"])
    else:
        is_peak = np.zeros(len(match_counts), dtype=int)
        if len(peaks) > 0:
            is_peak[peaks] = 1
        df = pd.DataFrame({
            "position": np.arange(len(match_counts)),
            "match_count": match_counts,
            "is_threshold_peak": is_peak,
            "threshold": threshold,
        })
    df.to_csv(out_csv, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Figure 1 motif match-count profiles")
    parser.add_argument("--manifest", required=True, help="Path to figure1_manifest.tsv")
    parser.add_argument("--fasta-dir", required=True, help="Directory with prepared Figure 1 FASTA files")
    parser.add_argument("--outdir", required=True, help="Output directory for motif results")
    parser.add_argument(
        "--min-score-frac",
        type=float,
        default=0.80,
        help="Peak threshold as fraction of the global peak, matching the old notebook logic",
    )
    parser.add_argument(
        "--motifs",
        nargs="*",
        default=["61bp", "40bp", "29bp", "37bp"],
        choices=["61bp", "40bp", "29bp", "37bp"],
        help="Motif classes to plot",
    )
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    plot_dir = ensure_dir(Path(outdir) / "profile_plots")
    csv_dir = ensure_dir(Path(outdir) / "profile_csv")
    summary_dir = ensure_dir(Path(outdir) / "summary")

    manifest = load_manifest(args.manifest)
    summary_rows: List[Dict[str, object]] = []

    for row in manifest.itertuples(index=False):
        if int(row.use_for_motif) != 1:
            continue

        fasta_path = build_fasta_path(args.fasta_dir, row.panel, row.species, row.record_id)
        if not fasta_path.exists():
            raise FileNotFoundError(f"Prepared FASTA not found: {fasta_path}")

        header, seq = read_single_fasta(fasta_path)
        monomers = get_monomer_dict_from_header_or_id(header)

        for motif_name in args.motifs:
            motif_seq = monomers[motif_name][0]
            motif_len = len(motif_seq)

            match_counts = compute_match_count_profile(seq=seq, motif_seq=motif_seq)
            peaks, threshold = detect_threshold_peaks(
                match_counts=match_counts,
                min_score_frac=args.min_score_frac,
            )

            base = f"{row.panel}_{row.species}_{sanitize_id_for_filename(row.record_id)}_{motif_name}"

            out_png = plot_dir / f"{base}.png"
            out_csv = csv_dir / f"{base}.csv"

            plot_match_count_profile(
                match_counts=match_counts,
                peaks=peaks,
                motif_name=motif_name,
                motif_len=motif_len,
                out_png=out_png,
                title_prefix=f"{row.panel}: {row.species} {row.record_id}",
            )
            save_profile_csv(
                out_csv=out_csv,
                match_counts=match_counts,
                peaks=peaks,
                threshold=threshold,
            )

            summary_rows.append(
                summarize_profile(
                    panel=row.panel,
                    species=row.species,
                    record_id=row.record_id,
                    header=header,
                    motif_name=motif_name,
                    motif_len=motif_len,
                    match_counts=match_counts,
                    peaks=peaks,
                    threshold=threshold,
                )
            )

            print(f"[OK] {row.panel}: {motif_name} profile -> {out_png}")
            print(f"[OK] {row.panel}: {motif_name} csv -> {out_csv}")

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = summary_dir / "figure1_motif_profile_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"[OK] summary -> {summary_csv}")


if __name__ == "__main__":
    main()