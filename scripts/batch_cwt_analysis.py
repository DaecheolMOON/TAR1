#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import traceback

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# -------------------- Helper Function: Reverse Complement -------------------- #
def get_reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    if not isinstance(seq, str) or len(seq) == 0:
        return ""

    seq = seq.upper()
    complement_map = str.maketrans("ATGC", "TACG")
    reversed_seq = seq[::-1]
    return reversed_seq.translate(complement_map)


# ------------------------ Sample Sequence Processing ------------------------ #
def process_sample_seq_for_folder(
    folder: Path,
    target_filename: str,
    required_cols,
    top_k: int,
    top_m: int,
) -> None:
    """
    For a given folder:
      1. Load target CSV (summary_runs.csv by default)
      2. Apply reverse complement to 'q' arms
      3. Compute frequencies of sample_seq by record_id
      4. Save pivot CSV for top_k sample_seq
      5. Save rank–frequency line plot
      6. Save stacked bar plot for top_m sample_seq
    """
    in_path = folder / target_filename
    print(f"\n[SampleSeq] Processing file: {in_path}")

    if not in_path.exists():
        print(f"  [SKIPPING] File not found: {in_path}")
        return

    # 1. Load data
    try:
        df = pd.read_csv(in_path)

        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"  [SKIPPING] Missing required columns: {missing}")
            return

        df["record_id"] = df["record_id"].fillna("Unknown").astype(str)
        df["sample_seq"] = df["sample_seq"].fillna("").astype(str)

    except Exception as e:
        print(f"  [ERROR] During data loading: {e}")
        traceback.print_exc()
        return

    # 2. Apply reverse complement for q-arms
    try:
        is_q_arm = df["record_id"].str.endswith("q", na=False)
        q_arm_count = is_q_arm.sum()

        if q_arm_count > 0:
            df.loc[is_q_arm, "sample_seq"] = df.loc[is_q_arm, "sample_seq"].apply(
                get_reverse_complement
            )
            print(f"  Applied reverse complement to {q_arm_count} 'q-arm' entries.")
        else:
            print("  No 'q-arm' entries found. No reverse complement applied.")

    except Exception as e:
        print(f"  [ERROR] During reverse complement step: {e}")
        traceback.print_exc()
        return

    # 3. Analyze sample_seq frequencies and create plots
    try:
        # 3.1 Filter non-empty sample_seq
        df_seqs = df[df["sample_seq"] != ""].copy()

        if df_seqs.empty:
            print("  No non-empty sample_seq entries found. Skipping.")
            return

        # 3.2 Count by (sample_seq, record_id)
        grouped_counts = (
            df_seqs.groupby(["sample_seq", "record_id"])
            .size()
            .reset_index(name="count")
        )

        # 3.3 Total count per sample_seq
        seq_totals = (
            grouped_counts.groupby("sample_seq")["count"]
            .sum()
            .sort_values(ascending=False)
        )

        if seq_totals.empty:
            print("  No sequence counts found. Skipping.")
            return

        # top_k for pivot and rank–frequency plot
        top_totals = seq_totals.head(top_k)
        top_seqs_list = top_totals.index.tolist()

        # 3.4 Build pivot table for top_k
        top_grouped_data = grouped_counts[
            grouped_counts["sample_seq"].isin(top_seqs_list)
        ]

        if top_grouped_data.empty:
            print("  No grouped data for top sequences. Skipping.")
            return

        df_pivot = top_grouped_data.pivot_table(
            index="sample_seq", columns="record_id", values="count", fill_value=0
        )

        # sort pivot rows in the same order as seq_totals
        df_pivot = df_pivot.loc[top_seqs_list]

        # 3.5 Save pivot CSV
        top_out_path = in_path.with_name(
            in_path.stem + f"_top{top_k}_sample_seqs_by_arm.csv"
        )
        df_pivot.to_csv(top_out_path)
        print(
            f"  Top {top_k} sample_seqs (pivoted by arm) saved to: {top_out_path}"
        )

        # 3.6 Rank–frequency line plot
        fig1 = plt.figure(figsize=(8, 4))
        plt.plot(
            np.arange(1, len(top_totals) + 1),
            top_totals.values,
            marker="o",
            linewidth=1.5,
        )
        plt.title(
            f"Top {top_k} Sample Seqs: Rank–Frequency\n(Source: {folder.name})"
        )
        plt.xlabel("Rank (1 = most frequent)")
        plt.ylabel("Count")
        plt.grid(True, linestyle="--", alpha=0.4)
        plt.tight_layout()

        # 3.7 Stacked bar plot for top_m
        # If top_m > number of available sequences, head will just return all
        df_pivot_topm = df_pivot.head(top_m)

        if df_pivot_topm.empty:
            print("  No data available for stacked bar plot. Skipping plot.")
            fig1_path = Path(str(top_out_path).replace(".csv", "_rank_freq_lineplot.png"))
            fig1.savefig(fig1_path, dpi=150)
            plt.close(fig1)
            print(f"  Saved plot:\n  - {fig1_path}")
            return

        num_arms = len(df_pivot_topm.columns)
        cmap = plt.get_cmap("gist_ncar", num_arms)

        fig2, ax = plt.subplots(figsize=(15, 8))
        df_pivot_topm.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            colormap=cmap,
            width=0.8,
        )

        ax.set_title(
            f"Top {min(top_m, len(df_pivot))} Sample Seqs: "
            f"Stacked Counts by Chromosome Arm\n(Source: {folder.name})"
        )
        ax.set_xlabel("Sample Sequence", fontsize=12)
        ax.set_ylabel("Total Count", fontsize=12)
        ax.tick_params(axis="x", rotation=90)

        ax.legend(
            title="Record ID (Arm)",
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            fontsize="small",
            ncol=1 if num_arms <= 30 else 2,
        )

        fig2.tight_layout()

        # 3.8 Save figures
        fig1_path = Path(
            str(top_out_path).replace(".csv", "_rank_freq_lineplot.png")
        )
        fig2_path = Path(
            str(top_out_path).replace(".csv", f"_top{top_m}_stacked_bar.png")
        )

        fig1.savefig(fig1_path, dpi=150)
        fig2.savefig(fig2_path, dpi=150, bbox_inches="tight")

        print(f"  Saved plots:\n  - {fig1_path}\n  - {fig2_path}")

        plt.close(fig1)
        plt.close(fig2)

    except Exception as e:
        print(f"  [ERROR] During 'sample_seq' analysis/plotting: {e}")
        traceback.print_exc()
        return


# ------------------------ end_base Statistics (v3) ------------------------ #
def process_end_base_stats_for_folder(
    folder: Path,
    target_filename: str,
    target_col: str = "end_base",
) -> None:
    """
    For a given folder:
      1. Load target CSV
      2. Adjust q-arm means based on 10M or 30k limit inferred from folder name
      3. Aggregate statistics (all + 'strong' signals)
      4. Compute density columns
      5. Save statistics CSV and print summary
    """
    in_path = folder / target_filename
    print(f"\n[EndBase] Analyzing {target_col} stats for: {in_path.name}")
    print(f"  (Source folder: {folder})")

    if not in_path.exists():
        print(f"  [SKIPPING] File not found: {in_path}")
        return

    # Determine q-arm limit from folder name
    q_arm_limit = 0
    folder_str = str(folder)
    if "10000000" in folder_str:
        q_arm_limit = 10_000_000
        print("  (q-arm adjustment limit set to 10,000,000)")
    elif "30000" in folder_str:
        q_arm_limit = 30_000
        print("  (q-arm adjustment limit set to 30,000)")
    else:
        print(
            "  [Warning] Could not determine q-arm limit "
            "(expected '10000000' or '30000' in folder path)."
        )

    stats_required_cols = [
        "record_id",
        target_col,
        "signal_type",
    ]

    # 1. Load data
    try:
        df_stats = pd.read_csv(in_path)

        missing_stats = [c for c in stats_required_cols if c not in df_stats.columns]
        if missing_stats:
            print(f"  [SKIPPING] Missing required columns for stats: {missing_stats}")
            return

        # Drop rows with missing record_id or signal_type
        df_stats = df_stats.dropna(subset=["record_id", "signal_type"])

        # Convert target_col to numeric
        df_stats[target_col] = pd.to_numeric(
            df_stats[target_col], errors="coerce"
        )
        df_stats = df_stats.dropna(subset=[target_col])

        if df_stats.empty:
            print(f"  [SKIPPING] No valid '{target_col}' data found.")
            return

    except Exception as e:
        print(f"  [ERROR] During data loading for stats: {e}")
        traceback.print_exc()
        return

    # 2. Aggregate statistics
    try:
        agg_funcs = ["mean", "std", "count"]

        mean_col_name = f"{target_col}_mean"
        std_col_name = f"{target_col}_std"
        count_col_name = "Sample_n"

        # 2.1 All signals
        arm_stats = df_stats.groupby("record_id")[target_col].agg(agg_funcs)
        arm_stats = arm_stats.rename(
            columns={
                "mean": mean_col_name,
                "std": std_col_name,
                "count": count_col_name,
            }
        )

        # 2.2 Strong signals only
        strong_mean_col = f"{target_col}_mean_strong"
        strong_std_col = f"{target_col}_std_strong"
        strong_count_col = "sample_n_strong"

        df_strong = df_stats[df_stats["signal_type"] == "strong"]

        if df_strong.empty:
            print("  No 'strong' signals found for aggregation.")
            arm_stats_strong = pd.DataFrame(
                columns=[strong_mean_col, strong_std_col, strong_count_col]
            )
        else:
            arm_stats_strong = df_strong.groupby("record_id")[target_col].agg(
                agg_funcs
            )
            arm_stats_strong = arm_stats_strong.rename(
                columns={
                    "mean": strong_mean_col,
                    "std": strong_std_col,
                    "count": strong_count_col,
                }
            )

        # 2.3 Merge all vs strong stats
        arm_stats = arm_stats.merge(
            arm_stats_strong, left_index=True, right_index=True, how="left"
        )

        # 2.4 Adjust q-arm mean values
        if q_arm_limit > 0:

            def adjust_q_arm(row):
                is_q_arm = isinstance(row.name, str) and row.name.endswith("q")

                mean_val = row[mean_col_name]
                if is_q_arm and pd.notna(mean_val):
                    row[mean_col_name] = abs(q_arm_limit - mean_val)

                strong_mean_val = row[strong_mean_col]
                if is_q_arm and pd.notna(strong_mean_val):
                    row[strong_mean_col] = abs(q_arm_limit - strong_mean_val)

                return row

            arm_stats = arm_stats.apply(adjust_q_arm, axis=1)
            print("  Applied 'q' arm adjustments to mean values.")

        # 2.5 Add density columns
        arm_stats["density"] = np.where(
            arm_stats[count_col_name] > 0,
            (arm_stats[mean_col_name] / arm_stats[count_col_name]) * 1000.0,
            np.nan,
        )

        arm_stats["density_strong"] = np.where(
            arm_stats[strong_count_col] > 0,
            (arm_stats[strong_mean_col] / arm_stats[strong_count_col]) * 1000.0,
            np.nan,
        )

        # Sort by record_id
        arm_stats = arm_stats.sort_index()

        if arm_stats.empty:
            print("  No 'record_id' groups found to aggregate.")
            return

        # 3. Save and print
        stats_out_path = in_path.with_name(
            in_path.stem + f"_{target_col}_stats_by_arm.csv"
        )
        arm_stats.to_csv(stats_out_path, float_format="%.4f")
        print(
            f"  Statistics (all + strong, q-arm adj, density) saved to: {stats_out_path}"
        )

        print(
            f"\n  --- {target_col.upper()} Statistics by Chromosome Arm "
            f"(Source: {folder.name}) ---"
        )

        with pd.option_context(
            "display.max_rows",
            100,
            "display.max_columns",
            15,
            "display.width",
            180,
        ):
            formatters = {
                mean_col_name: "{:,.2f}".format,
                std_col_name: "{:,.2f}".format,
                count_col_name: "{:,.0f}".format,
                "density": "{:.4f}".format,
                strong_mean_col: "{:,.2f}".format,
                strong_std_col: "{:,.2f}".format,
                strong_count_col: "{:,.0f}".format,
                "density_strong": "{:.4f}".format,
            }

            print(arm_stats.to_string(formatters=formatters, na_rep="---"))
        print("  ------------------------------------------------" + "-" * len(folder.name))

    except Exception as e:
        print(f"  [ERROR] During '{target_col}' stats aggregation: {e}")
        traceback.print_exc()
        return


# ------------------------------ Main Entrypoint ------------------------------ #
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Batch process CWT summary CSV files:\n"
            "  1) sample_seq analysis (reverse complement, top-k, plots)\n"
            "  2) end_base statistics with q-arm adjustment (v3)"
        )
    )
    parser.add_argument(
        "folders",
        nargs="+",
        help="One or more folder paths containing the target CSV file.",
    )
    parser.add_argument(
        "-f",
        "--target-filename",
        default="summary_runs.csv",
        help="Name of the CSV file inside each folder (default: summary_runs.csv).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=100,
        help=(
            "Number of top sample sequences used for pivot table and rank–frequency "
            "plot (default: 100)."
        ),
    )
    parser.add_argument(
        "--top-m",
        type=int,
        default=30,
        help=(
            "Number of top sample sequences used for stacked bar plot "
            "(default: 30)."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # If top_m > top_k, clamp top_m to top_k
    if args.top_m > args.top_k:
        print(
            f"[Warning] top-m ({args.top_m}) is greater than top-k ({args.top_k}). "
            f"Using top-m = top-k."
        )
        args.top_m = args.top_k

    base_folders = [Path(p).expanduser().resolve() for p in args.folders]
    target_filename = args.target_filename
    top_k = args.top_k
    top_m = args.top_m

    required_cols = ["sample_seq", "record_id"]

    print(
        f"--- Starting batch processing for {len(base_folders)} folders "
        f"(target file: {target_filename}) ---"
    )

    # 1) sample_seq analysis
    for folder in base_folders:
        process_sample_seq_for_folder(
            folder=folder,
            target_filename=target_filename,
            required_cols=required_cols,
            top_k=top_k,
            top_m=top_m,
        )

    print("\n--- Sample sequence batch processing finished. ---")

    # 2) end_base statistics (v3)
    print("\n" + "=" * 50)
    print("--- Starting batch processing for 'end_base' statistics (v3) ---")
    print("=" * 50)

    for folder in base_folders:
        process_end_base_stats_for_folder(
            folder=folder,
            target_filename=target_filename,
            target_col="end_base",
        )

    print("\n--- 'end_base' statistics processing (v3) finished. ---")


if __name__ == "__main__":
    main()
