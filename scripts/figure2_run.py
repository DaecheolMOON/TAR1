from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Figure 2A-E pipeline")
    parser.add_argument("--manifest", required=True, help="Path to figure2_manifest.tsv")
    parser.add_argument("--figure2e-manifest", default="metadata/figure2E_manifest.tsv", help="Path to figure2E manifest")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--base-outdir", required=True, help="Base output directory")
    args = parser.parse_args()

    base_outdir = Path(args.base_outdir)
    fasta_dir = base_outdir / "windows"
    a_dir = base_outdir / "figure2A"
    b_dir = base_outdir / "figure2B"
    c_dir = base_outdir / "figure2C"
    d_dir = base_outdir / "figure2D"
    e_dir = base_outdir / "figure2E"

    py = sys.executable
    script_dir = Path(__file__).resolve().parent
    data_root = args.data_root if args.data_root is not None else "./data"

    run_cmd([
        py,
        str(script_dir / "figure2_prepare_windows.py"),
        "--manifest", args.manifest,
        "--data-root", data_root,
        "--outdir", str(fasta_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure2_scalogram.py"),
        "--manifest", args.manifest,
        "--data-root", data_root,
        "--fasta-dir", str(fasta_dir),
        "--outdir", str(a_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure2_distribution.py"),
        "--data-root", data_root,
        "--outdir", str(b_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure2_density.py"),
        "--summary-csv", str(b_dir / "tar1_counts_summary_all_species.csv"),
        "--outdir", str(c_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure2D_centromere.py"),
        "--data-root", data_root,
        "--outdir", str(d_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure2E_whole_chromosome.py"),
        "--manifest", args.figure2e_manifest,
        "--data-root", data_root,
        "--outdir", str(e_dir),
    ])


if __name__ == "__main__":
    main()