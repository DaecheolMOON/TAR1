#!/usr/bin/env python3
"""
Compute k-mer-based cosine similarity between sample_seq strings,
save full and top-N similarity matrices, and perform network and
cluster analysis based on cosine similarity.

Requirements:
    pip install pandas numpy matplotlib networkx
"""

import argparse
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx


ALLOWED_BASES = set("ACGT")


# ----------------------------- Utility functions ----------------------------- #

def clean_seq(s: str) -> str:
    """Uppercase and keep only A/C/G/T characters."""
    return "".join([ch for ch in str(s).upper() if ch in ALLOWED_BASES])


def build_kmer_vocab(k: int):
    """Build k-mer vocabulary and index dict."""
    vocab = ["".join(p) for p in product("ACGT", repeat=k)]
    idx = {mer: i for i, mer in enumerate(vocab)}
    return vocab, idx


def kmer_vec(s: str, k: int, idx: dict, d: int) -> np.ndarray:
    """Return L2-normalized k-mer count vector for sequence s."""
    v = np.zeros(d, dtype=float)
    if len(s) < k:
        return v
    for i in range(len(s) - k + 1):
        mer = s[i:i + k]
        j = idx.get(mer, None)
        if j is not None:
            v[j] += 1.0
    norm = np.linalg.norm(v)
    if norm > 0:
        v /= norm
    return v


def build_labels_from_seqs(seqs):
    """
    Build row/column labels from sequences:
    first 5 bases + running index, e.g. 'TATCT_01'.
    """
    labels = []
    counter = {}
    for s in seqs:
        tag = s[:5] if len(s) >= 5 else s
        cnt = counter.get(tag, 0) + 1
        counter[tag] = cnt
        labels.append(f"{tag}_{cnt:02d}")
    return labels


def ensure_input_path(path_str: str) -> Path:
    """Resolve input path; if no suffix and '.csv' exists, use that."""
    p = Path(path_str).expanduser().resolve()
    if p.exists():
        return p
    if p.suffix == "":
        alt = p.with_suffix(".csv")
        if alt.exists():
            return alt
    raise FileNotFoundError(f"Input file not found: {p}")


# ------------------------------- Core routines ------------------------------- #

def compute_full_similarity(df: pd.DataFrame, seq_col: str, k: int, out_base: Path,
                            preview_size: int, top_pairs: int):
    """
    Compute full cosine similarity matrix for all sequences in df[seq_col].
    Save CSV and print preview + most similar/dissimilar pairs.
    """
    seqs_raw = df[seq_col].astype(str).tolist()
    seqs = [clean_seq(s) for s in seqs_raw]
    n = len(seqs)

    if n == 0:
        raise ValueError("No sequences found in input file.")
    if n == 1:
        print("Only 1 sequence found. Pairwise similarity cannot be computed.")

    print(f"\n[Full] Number of sequences: {n}")
    print(f"[Full] k-mer size: k={k}")

    vocab, idx = build_kmer_vocab(k)
    d = len(vocab)
    print(f"[Full] k-mer vocabulary size: {d}")

    X = np.vstack([kmer_vec(s, k=k, idx=idx, d=d) for s in seqs])  # (n, d)
    S = X @ X.T  # cosine similarity matrix (n, n)

    labels = build_labels_from_seqs(seqs)
    sim_df = pd.DataFrame(S, index=labels, columns=labels)

    out_csv = out_base.with_name(out_base.stem + f"_cosine_k{k}.csv")
    sim_df.to_csv(out_csv, float_format="%.4f")
    print(f"[Full] Cosine similarity matrix saved to:\n  {out_csv}")

    # Preview
    p = min(preview_size, n)
    print(f"\n[Full] Preview (top-left {p}x{p}):")
    with pd.option_context("display.max_rows", p, "display.max_columns", p):
        print(sim_df.iloc[:p, :p])

    # Most similar / dissimilar pairs
    if n > 1:
        tri_idx = np.triu_indices(n, k=1)
        vals = S[tri_idx]
        order_desc = np.argsort(vals)[::-1]
        order_asc = np.argsort(vals)

        m = min(top_pairs, len(order_desc))

        print(f"\n[Full] Most similar pairs (top {m}, excluding self):")
        for r in range(m):
            i = tri_idx[0][order_desc[r]]
            j = tri_idx[1][order_desc[r]]
            print(f"  {labels[i]}  ~  {labels[j]} : cosine={S[i, j]:.4f}")

        print(f"\n[Full] Most dissimilar pairs (bottom {m}):")
        for r in range(m):
            i = tri_idx[0][order_asc[r]]
            j = tri_idx[1][order_asc[r]]
            print(f"  {labels[i]}  ~  {labels[j]} : cosine={S[i, j]:.4f}")

    return sim_df, labels, S, seqs


def compute_top_n_similarity(sim_df_full: pd.DataFrame,
                             labels_full,
                             S_full: np.ndarray,
                             df_full: pd.DataFrame,
                             seq_col: str,
                             top_n: int,
                             k: int,
                             out_base: Path):
    """
    Extract top-N sequences (by row order in df_full) and build their
    cosine similarity submatrix.
    Save CSV and print preview.
    """
    n_total = len(labels_full)
    top_n = min(top_n, n_total)

    if top_n < 1:
        print("\n[TopN] top-n < 1. Skipping top-N analysis.")
        return None, None, None, None

    if top_n == 1:
        print("\n[TopN] Only 1 sequence selected. Skipping pairwise analysis.")
        df_top = df_full.head(1).copy()
        labels_top = [labels_full[0]]
        S_top = S_full[:1, :1]
        sim_df_top = pd.DataFrame(S_top, index=labels_top, columns=labels_top)
        return sim_df_top, labels_top, df_top, top_n

    print(f"\n[TopN] Using top {top_n} sequences (by row order).")

    df_top = df_full.head(top_n).copy()
    labels_top = labels_full[:top_n]
    S_top = S_full[:top_n, :top_n]

    sim_df_top = pd.DataFrame(S_top, index=labels_top, columns=labels_top)
    out_csv_top = out_base.with_name(
        out_base.stem + f"_cosine_k{k}_top{top_n}.csv"
    )
    sim_df_top.to_csv(out_csv_top, float_format="%.4f")
    print(f"[TopN] Cosine similarity matrix (top {top_n}) saved to:\n  {out_csv_top}")

    p = min(top_n, 10)
    print(f"\n[TopN] Preview (top-left {p}x{p}):")
    with pd.option_context("display.max_rows", p, "display.max_columns", p):
        print(sim_df_top.iloc[:p, :p])

    # Simple most similar/dissimilar pairs inside top-N
    tri_idx = np.triu_indices(top_n, k=1)
    vals = S_top[tri_idx]
    if len(vals) > 0:
        order_desc = np.argsort(vals)[::-1]
        order_asc = np.argsort(vals)
        m = min(10, len(order_desc))

        print(f"\n[TopN] Most similar pairs (top {m}):")
        for r in range(m):
            i = tri_idx[0][order_desc[r]]
            j = tri_idx[1][order_desc[r]]
            print(f"  {labels_top[i]}  ~  {labels_top[j]} : cosine={S_top[i, j]:.4f}")

        print(f"\n[TopN] Most dissimilar pairs (bottom {m}):")
        for r in range(m):
            i = tri_idx[0][order_asc[r]]
            j = tri_idx[1][order_asc[r]]
            print(f"  {labels_top[i]}  ~  {labels_top[j]} : cosine={S_top[i, j]:.4f}")
    else:
        print("[TopN] Not enough pairs to compute similarity (top-n <= 1).")

    return sim_df_top, labels_top, df_top, top_n


def build_network_and_save(sim_df_top: pd.DataFrame,
                           labels_top,
                           network_threshold: float,
                           k: int,
                           top_n: int,
                           out_base: Path):
    """
    Build a similarity network from sim_df_top using the given threshold.
    Save the graph as a PNG image.
    """
    print("\n[Network] Building network graph "
          f"(threshold >= {network_threshold:.3f})")

    G = nx.Graph()
    labels = list(sim_df_top.index)
    G.add_nodes_from(labels)

    edges_to_add = []
    edge_weights_list = []
    edge_alphas_list = []
    edge_labels = {}

    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            label1 = labels[i]
            label2 = labels[j]
            weight = float(sim_df_top.loc[label1, label2])
            if weight >= network_threshold:
                edges_to_add.append((label1, label2))
                G.add_edge(label1, label2, weight=weight)
                edge_weights_list.append(weight * 5.0)
                edge_alphas_list.append(weight * 0.8 + 0.2)
                edge_labels[(label1, label2)] = f"{weight:.2f}"

    print(f"[Network] Nodes: {G.number_of_nodes()}, "
          f"Edges (>= {network_threshold:.3f}): {G.number_of_edges()}")

    if G.number_of_edges() == 0:
        print("[Network] No edges above threshold. Skipping drawing.")
        return

    edge_weights = np.array(edge_weights_list)
    edge_alphas_raw = np.array(edge_alphas_list)
    edge_alphas = np.clip(edge_alphas_raw, 0.0, 1.0)

    plt.figure(figsize=(14, 12))
    pos = nx.spring_layout(G, weight="weight", k=1.0, iterations=50, seed=42)

    nx.draw_networkx_nodes(
        G, pos,
        node_size=2500,
        node_color="skyblue",
        alpha=0.9,
    )

    nx.draw_networkx_edges(
        G, pos,
        edgelist=edges_to_add,
        width=edge_weights,
        alpha=edge_alphas,
        edge_color="gray",
    )

    nx.draw_networkx_labels(
        G, pos,
        font_size=10,
        font_weight="bold",
    )

    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_size=9,
        font_color="black",
    )

    thr_str = f"{network_threshold:.2f}".replace(".", "p")
    out_png = out_base.with_name(
        out_base.stem + f"_network_k{k}_top{top_n}_thr{thr_str}.png"
    )

    plt.title(
        f"Top {top_n} Sequence Similarity Network "
        f"(Cosine k={k}, thr >= {network_threshold:.2f})",
        fontsize=16,
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"[Network] Network figure saved to:\n  {out_png}")


def find_clusters(sim_df_top: pd.DataFrame, cluster_threshold: float):
    """
    Build a graph using cosine >= cluster_threshold and return connected components.
    """
    print("\n[Cluster] Finding clusters "
          f"(threshold >= {cluster_threshold:.3f})")

    G_cluster = nx.Graph()
    labels = list(sim_df_top.index)
    G_cluster.add_nodes_from(labels)

    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            label1 = labels[i]
            label2 = labels[j]
            weight = float(sim_df_top.loc[label1, label2])
            if weight >= cluster_threshold:
                G_cluster.add_edge(label1, label2, weight=weight)

    connected_components = list(nx.connected_components(G_cluster))
    clusters = sorted(connected_components, key=len, reverse=True)

    print(f"[Cluster] Found {len(clusters)} clusters "
          f"(families) at threshold >= {cluster_threshold:.3f}")
    return clusters


def cluster_analysis_basic(sim_df_top: pd.DataFrame,
                           clusters,
                           df_top: pd.DataFrame,
                           labels_top,
                           seq_col: str):
    """
    For each cluster, select a representative sequence and report its
    raw sequence and main chromosome arms (non-zero columns starting with 'chr').
    """
    print("\n[Cluster] Basic representative analysis")

    label_to_index = {label: i for i, label in enumerate(labels_top)}
    chr_cols = [c for c in df_top.columns if c.startswith("chr")]

    for i, cluster_nodes in enumerate(clusters):
        cluster_nodes_list = list(cluster_nodes)
        num_members = len(cluster_nodes_list)
        cluster_name = f"Cluster {i + 1}"

        print("-" * 30)
        print(f"{cluster_name} (Members: {num_members})")
        print("Members: " + ", ".join(cluster_nodes_list))

        # Find representative label
        if num_members == 1:
            rep_label = cluster_nodes_list[0]
            print(f"Representative: {rep_label} (single member cluster)")
        else:
            avg_similarities = {}
            for node_i in cluster_nodes_list:
                total_sim = 0.0
                for node_j in cluster_nodes_list:
                    if node_i == node_j:
                        continue
                    total_sim += float(sim_df_top.loc[node_i, node_j])
                avg_similarities[node_i] = total_sim / (num_members - 1)
            rep_label = max(avg_similarities, key=avg_similarities.get)
            print(
                "Representative: "
                f"{rep_label} (Avg. similarity: {avg_similarities[rep_label]:.4f})"
            )

        # Representative raw sequence and arms
        rep_index = label_to_index.get(rep_label, None)
        if rep_index is None:
            print(f"Error: could not map label to index: {rep_label}")
            print("-" * 30)
            continue

        rep_row = df_top.iloc[rep_index]
        rep_seq_raw = rep_row[seq_col]
        print(f"Rep. sequence: {rep_seq_raw}")

        if not chr_cols:
            print("Associated arms: no columns starting with 'chr' found.")
        else:
            arm_data = rep_row[chr_cols]
            active_arms = arm_data[arm_data > 0].sort_values(ascending=False)
            if active_arms.empty:
                print("Associated arms: none with count > 0.")
            else:
                print("Associated arms (count):")
                for arm, count in active_arms.items():
                    print(f"  - {arm}: {int(count)}")

        print("-" * 30 + "\n")


def cluster_analysis_full_arms(sim_df_top: pd.DataFrame,
                               clusters,
                               df_top: pd.DataFrame,
                               labels_top,
                               seq_col: str):
    """
    For each cluster, report representative sequence and the sum of
    chromosome-arm counts (for all members of the cluster).
    """
    print("\n[Cluster] Full arm summary analysis (cluster-wise sums)")

    label_to_index = {label: i for i, label in enumerate(labels_top)}
    chr_cols = [c for c in df_top.columns if c.startswith("chr")]

    for i, cluster_nodes in enumerate(clusters):
        cluster_nodes_list = list(cluster_nodes)
        num_members = len(cluster_nodes_list)
        cluster_name = f"Cluster {i + 1}"

        print("-" * 30)
        print(f"{cluster_name} (Members: {num_members})")
        print("Members: " + ", ".join(cluster_nodes_list))

        # Representative label (same logic as basic analysis)
        if num_members == 1:
            rep_label = cluster_nodes_list[0]
            print(f"Representative: {rep_label} (single member cluster)")
        else:
            avg_similarities = {}
            for node_i in cluster_nodes_list:
                total_sim = 0.0
                for node_j in cluster_nodes_list:
                    if node_i == node_j:
                        continue
                    total_sim += float(sim_df_top.loc[node_i, node_j])
                avg_similarities[node_i] = total_sim / (num_members - 1)
            rep_label = max(avg_similarities, key=avg_similarities.get)
            print(
                "Representative: "
                f"{rep_label} (Avg. similarity: {avg_similarities[rep_label]:.4f})"
            )

        # Representative sequence
        rep_index = label_to_index.get(rep_label, None)
        if rep_index is not None:
            rep_seq_raw = df_top.iloc[rep_index][seq_col]
            print(f"Rep. sequence: {rep_seq_raw}")
        else:
            print(f"Error: could not map label to index: {rep_label}")

        # Cluster-arm totals
        if not chr_cols:
            print("Total associated arms: no 'chr*' columns found.")
            print("-" * 30 + "\n")
            continue

        cluster_indices = []
        for label in cluster_nodes_list:
            idx = label_to_index.get(label, None)
            if idx is not None:
                cluster_indices.append(idx)

        if not cluster_indices:
            print("Could not find any indices for members in this cluster.")
            print("-" * 30 + "\n")
            continue

        cluster_rows = df_top.iloc[cluster_indices]
        cluster_arm_data = cluster_rows[chr_cols]
        cluster_arm_totals = cluster_arm_data.sum(axis=0)
        active_arms = cluster_arm_totals[cluster_arm_totals > 0].sort_values(
            ascending=False
        )

        if active_arms.empty:
            print("Total associated arms for cluster: none with count > 0.")
        else:
            print("Total associated arms for cluster (sum across members):")
            for arm, count in active_arms.items():
                print(f"  - {arm}: {int(count)}")

        print("-" * 30 + "\n")


# --------------------------------- CLI parsing -------------------------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute k-mer-based cosine similarity between sample_seq strings, "
            "build a similarity network, and perform cluster analysis.\n\n"
            "Typical input: summary_runs_topXXX_sample_seqs_by_arm.csv"
        )
    )
    parser.add_argument(
        "input_csv",
        help="Path to input CSV (pivot table with 'sample_seq' column).",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=3,
        help="k-mer size (default: 3).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=30,
        help="Number of top sequences (by row order) to use for network/cluster analysis (default: 30).",
    )
    parser.add_argument(
        "--network-threshold",
        type=float,
        default=0.5,
        help="Cosine similarity threshold for drawing network edges (default: 0.5).",
    )
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=0.7,
        help="Cosine similarity threshold for defining clusters (connected components) (default: 0.7).",
    )
    parser.add_argument(
        "--preview-size",
        type=int,
        default=10,
        help="Size of top-left block to preview for the full similarity matrix (default: 10).",
    )
    parser.add_argument(
        "--top-pairs",
        type=int,
        default=10,
        help="Number of most similar/dissimilar pairs to print (default: 10).",
    )
    return parser.parse_args()


# ------------------------------------ Main ------------------------------------ #

def main():
    args = parse_args()

    input_path = ensure_input_path(args.input_csv)
    print(f"Input file: {input_path}")

    df = pd.read_csv(input_path)

    if "sample_seq" in df.columns:
        seq_col = "sample_seq"
    else:
        # Fallback: use first column as sample_seq
        seq_col = df.columns[0]
        print(
            f"[Warning] 'sample_seq' column not found. "
            f"Using first column as sequence column: '{seq_col}'"
        )

    # Full similarity
    sim_df_full, labels_full, S_full, seqs_full = compute_full_similarity(
        df=df,
        seq_col=seq_col,
        k=args.k,
        out_base=input_path,
        preview_size=args.preview_size,
        top_pairs=args.top_pairs,
    )

    # Top-N similarity
    sim_df_top, labels_top, df_top, top_n_used = compute_top_n_similarity(
        sim_df_full=sim_df_full,
        labels_full=labels_full,
        S_full=S_full,
        df_full=df,
        seq_col=seq_col,
        top_n=args.top_n,
        k=args.k,
        out_base=input_path,
    )

    if sim_df_top is None or top_n_used is None or top_n_used < 2:
        print("\n[Main] Top-N analysis not performed (top_n < 2). Done.")
        return

    # Network
    build_network_and_save(
        sim_df_top=sim_df_top,
        labels_top=labels_top,
        network_threshold=args.network_threshold,
        k=args.k,
        top_n=top_n_used,
        out_base=input_path,
    )

    # Cluster analyses
    clusters = find_clusters(sim_df_top, cluster_threshold=args.cluster_threshold)

    cluster_analysis_basic(
        sim_df_top=sim_df_top,
        clusters=clusters,
        df_top=df_top,
        labels_top=labels_top,
        seq_col=seq_col,
    )

    cluster_analysis_full_arms(
        sim_df_top=sim_df_top,
        clusters=clusters,
        df_top=df_top,
        labels_top=labels_top,
        seq_col=seq_col,
    )

    print("\n[Main] All analyses complete.")


if __name__ == "__main__":
    main()
