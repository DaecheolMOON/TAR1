from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Figure 3A-B pipeline")
    parser.add_argument("--figure3a-manifest", required=True, help="Path to figure3A_manifest.tsv")
    parser.add_argument("--data-root", default=None, help="Project data root")
    parser.add_argument("--base-outdir", required=True, help="Base output directory")
    args = parser.parse_args()

    base_outdir = Path(args.base_outdir)
    a_dir = base_outdir / "figure3A"
    b_dir = base_outdir / "figure3B"

    py = sys.executable
    script_dir = Path(__file__).resolve().parent
    data_root = args.data_root if args.data_root is not None else "./data"

    run_cmd([
        py,
        str(script_dir / "figure3A_atcg_profile.py"),
        "--manifest", args.figure3a_manifest,
        "--data-root", data_root,
        "--outdir", str(a_dir),
    ])

    run_cmd([
        py,
        str(script_dir / "figure3B_cap_association.py"),
        "--data-root", data_root,
        "--outdir", str(b_dir),
    ])


if __name__ == "__main__":
    main()