from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all Figure 1 preparation, CWT, and motif steps")
    parser.add_argument("--manifest", required=True, help="Path to figure1_manifest.tsv")
    parser.add_argument("--data-root", default=os.environ.get("TAR1_DATA_ROOT", "./data"), help="Project data root")
    parser.add_argument("--base-outdir", default="results/figure1", help="Base output directory")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    base_out = Path(args.base_outdir)
    fasta_dir = base_out / "fasta"
    cwt_dir = base_out / "cwt"
    motif_dir = base_out / "motif"

    run_cmd([
        sys.executable,
        str(script_dir / "figure1_prepare_loci.py"),
        "--manifest", args.manifest,
        "--data-root", args.data_root,
        "--outdir", str(fasta_dir),
    ])

    run_cmd([
        sys.executable,
        str(script_dir / "figure1_cwt.py"),
        "--manifest", args.manifest,
        "--fasta-dir", str(fasta_dir),
        "--outdir", str(cwt_dir),
    ])

    run_cmd([
        sys.executable,
        str(script_dir / "figure1_motif.py"),
        "--manifest", args.manifest,
        "--fasta-dir", str(fasta_dir),
        "--outdir", str(motif_dir),
    ])

    print("[DONE] Figure 1 outputs generated under", base_out)


if __name__ == "__main__":
    main()
