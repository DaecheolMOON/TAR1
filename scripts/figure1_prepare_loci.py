from __future__ import annotations

import argparse
import os
from pathlib import Path

from figure1_common import (
    collect_matching_records,
    ensure_dir,
    find_record_in_fasta,
    load_manifest,
    resolve_tar1_blocks_fasta,
    sanitize_id_for_filename,
    summarize_matches,
    write_single_fasta,
)


def build_output_name(panel: str, species: str, record_id: str) -> str:
    return f"{panel}_{species}_{sanitize_id_for_filename(record_id)}.fa"


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Figure 1 loci from species tar1_blocks.fa files")
    parser.add_argument("--manifest", required=True, help="Path to figure1_manifest.tsv")
    parser.add_argument("--data-root", default=os.environ.get("TAR1_DATA_ROOT", "./data"), help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output FASTA directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    manifest = load_manifest(args.manifest)

    for row in manifest.itertuples(index=False):
        tar1_fasta = resolve_tar1_blocks_fasta(row.species, args.data_root)

        matches = collect_matching_records(
            fasta_path=tar1_fasta,
            record_id=row.record_id,
            match_mode=row.match_mode,
        )

        if not matches:
            raise ValueError(
                f"No record matched record_id={row.record_id!r} "
                f"with match_mode={row.match_mode!r} in {tar1_fasta}"
            )

        if len(matches) > 1:
            summary_df = summarize_matches(matches, requested_record_id=row.record_id)
            print(
                f"[INFO] {row.panel}: multiple matches found for {row.species} {row.record_id} "
                f"(match_mode={row.match_mode}, min_length={row.min_length}, selection_mode={row.selection_mode})"
            )
            print(summary_df.to_string(index=False))

        header, seq = find_record_in_fasta(
            fasta_path=tar1_fasta,
            record_id=row.record_id,
            match_mode=row.match_mode,
            min_length=row.min_length,
            selection_mode=row.selection_mode,
        )

        out_name = build_output_name(row.panel, row.species, row.record_id)
        out_path = outdir / out_name
        write_single_fasta(header=header, sequence=seq, out_path=out_path)

        print(
            f"[OK] {row.panel}: {row.species} {row.record_id} -> {out_path} "
            f"(selected_header={header!r}, length={len(seq)})"
        )


if __name__ == "__main__":
    main()