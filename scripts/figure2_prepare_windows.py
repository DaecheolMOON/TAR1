from __future__ import annotations

import argparse
from pathlib import Path

from figure2_common import (
    build_telomere_window_sequence,
    ensure_dir,
    load_manifest,
    sanitize_id_for_filename,
    write_single_fasta,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Figure 2A telomere-proximal window FASTA files")
    parser.add_argument("--manifest", required=True, help="Path to figure2_manifest.tsv")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--outdir", required=True, help="Output FASTA directory")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    manifest = load_manifest(args.manifest)

    for row in manifest.itertuples(index=False):
        seq = build_telomere_window_sequence(
            species=row.species,
            arm=row.arm,
            window_bp=int(row.window_bp),
            data_root=args.data_root,
            orient_telomere_left=True,
        )
        out_name = f"{row.panel}_{row.species}_{sanitize_id_for_filename(row.arm)}_{int(row.window_bp)}bp.fa"
        out_path = Path(outdir) / out_name
        write_single_fasta(
            header=f"{row.arm}|window_bp={int(row.window_bp)}|species={row.species}",
            sequence=seq,
            out_path=out_path,
        )
        print(f"[OK] {row.panel}: wrote {out_path}")


if __name__ == "__main__":
    main()