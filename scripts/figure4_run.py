from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Figure 4 pipeline")
    parser.add_argument("--data-root", default="./data")
    parser.add_argument("--base-outdir", required=True)
    args = parser.parse_args()

    py = sys.executable
    script_dir = Path(__file__).resolve().parent
    base_outdir = Path(args.base_outdir)

    pos_dir = base_outdir / "figure4_positions"
    cor_dir = base_outdir / "figure4_correlations"

    run_cmd([
        py,
        str(script_dir / "figure4_positions.py"),
        "--data-root", args.data_root,
        "--outdir", str(pos_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure4_correlations.py"),
        "--data-root", args.data_root,
        "--outdir", str(cor_dir),
    ])

    print("[DONE] Figure 4 pipeline completed.")


if __name__ == "__main__":
    main()