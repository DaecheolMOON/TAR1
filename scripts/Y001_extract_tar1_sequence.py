#!/usr/bin/env python3
import argparse
import subprocess
import sys
import shutil
import tempfile
import os
import glob


def check_dependencies(tools):
    """Check if required command-line tools are installed and available in PATH."""
    print("--- Checking for required tools... ---")
    all_found = True
    for tool in tools:
        if not shutil.which(tool):
            print(
                f"Error: Required tool '{tool}' is not installed or not in your PATH.",
                file=sys.stderr,
            )
            all_found = False
    if not all_found:
        print("Please install the missing tools and try again.", file=sys.stderr)
        sys.exit(1)
    print("All required tools found.\n")


def process_tbl_to_bed_data(tbl_content):
    """
    Parse nhmmer .tbl output and convert it into BED-like records.

    This function reproduces the logic of:
        grep -v "^#" | awk 'BEGIN{OFS="\t"} { a=$9; b=$10; start=(a<b?a:b)-1; end=(a<b?b:a); print $1, start, end }' \
        | sort -k1,1 -k2,2n
    """
    bed_lines = []
    for line in tbl_content.strip().split("\n"):
        # Skip comment lines (starting with "#")
        if line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 10:
            continue

        target_name = fields[0]
        pos1 = int(fields[8])
        pos2 = int(fields[9])

        # BED format uses 0-based start and 1-based end
        # The start/end values below follow this convention.
        start = min(pos1, pos2) - 1
        end = max(pos1, pos2)

        bed_lines.append((target_name, start, end))

    # Sort by chromosome name, then by start coordinate
    bed_lines.sort(key=lambda x: (x[0], x[1]))

    return [f"{name}\t{start}\t{end}" for name, start, end in bed_lines]


def run_command(command, step_name):
    """Run a subprocess command and handle errors."""
    print(f"--- Running: {step_name} ---")
    print(f"Command: {' '.join(command)}")
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        print(f"--- {step_name} successful. ---\n")
        return result
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred during: {step_name}", file=sys.stderr)
        print(f"Command failed with exit code {e.returncode}", file=sys.stderr)
        print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
        print("\n--- STDOUT ---", file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print("\n--- STDERR ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        # Do not terminate the whole script here; return False so that
        # the caller can decide how to handle failures per file.
        return False
    except FileNotFoundError:
        print(
            f"Error: Command '{command[0]}' not found. Is it installed and in your PATH?",
            file=sys.stderr,
        )
        sys.exit(1)


def run_pipeline(hmm_file, fasta_file, output_fasta_file, cpu, incE):
    """Run the nhmmer -> BED -> bedtools getfasta pipeline for a single FASTA file."""
    print(f"Processing file: {os.path.basename(fasta_file)}")

    # Use a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_tbl_path = os.path.join(temp_dir, "hits.tbl")
        temp_bed_path = os.path.join(temp_dir, "blocks.bed")

        # 1. nhmmer search
        nhmmer_cmd = [
            "nhmmer",
            "--cpu",
            str(cpu),
            "--tblout",
            temp_tbl_path,
            "--incE",
            str(incE),
            hmm_file,
            fasta_file,
        ]
        if not run_command(nhmmer_cmd, "nhmmer search"):
            return False

        # 2. Convert nhmmer .tbl output to BED
        print("--- Running: Process TBL to BED ---")
        try:
            with open(temp_tbl_path, "r") as f:
                tbl_content = f.read()

            bed_data = process_tbl_to_bed_data(tbl_content)

            with open(temp_bed_path, "w") as f:
                f.write("\n".join(bed_data))
            print("--- TBL to BED conversion successful. ---\n")
        except Exception as e:
            print(f"Error processing TBL file: {e}", file=sys.stderr)
            return False

        # 3. bedtools getfasta
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)

        bedtools_cmd = [
            "bedtools",
            "getfasta",
            "-fi",
            fasta_file,
            "-bed",
            temp_bed_path,
            # Note: bedtools getfasta uses the 4th BED column as a name if present.
            # We only have 3 columns here, so we omit the -name option and let
            # bedtools construct headers from coordinates.
            "-fo",
            output_fasta_file,
        ]
        if not run_command(bedtools_cmd, "bedtools getfasta"):
            return False

    print(f"Successfully processed and saved to: {output_fasta_file}\n")
    return True


def main():
    """Main entry point: parse arguments and run the pipeline on all FASTA files."""
    parser = argparse.ArgumentParser(
        description=(
            "Run a nhmmer -> bedtools pipeline on multiple FASTA files in a directory.\n\n"
            "For each input FASTA:\n"
            "  1) search with nhmmer using the given HMM,\n"
            "  2) convert .tbl output to BED,\n"
            "  3) extract matching sequences with bedtools getfasta."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--hmm",
        required=True,
        help="Path to the input profile HMM file (e.g., TAR1.hmm)",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to the input directory containing reference FASTA files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the base directory for the output files.",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=6,
        help="Number of CPU cores to use for nhmmer (default: 6)",
    )
    parser.add_argument(
        "--incE",
        type=float,
        default=1e-5,
        help="Inclusion E-value threshold for nhmmer (default: 1e-5)",
    )
    args = parser.parse_args()

    # 0. Check dependencies
    check_dependencies(["nhmmer", "bedtools"])

    # Check input FASTA directory
    if not os.path.isdir(args.fasta):
        print(
            f"Error: Input FASTA path is not a directory: {args.fasta}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create output base directory
    os.makedirs(args.output, exist_ok=True)

    # Collect FASTA files
    fasta_files = glob.glob(os.path.join(args.fasta, "*.fasta"))
    if not fasta_files:
        print(
            f"No .fasta files found in directory: {args.fasta}",
            file=sys.stderr,
        )
        sys.exit(0)

    print(f"Found {len(fasta_files)} FASTA files to process.\n")

    success_count = 0
    fail_count = 0

    # Run pipeline for each FASTA file
    for fasta_path in fasta_files:
        filename = os.path.basename(fasta_path)

        # Extract sample ID from filename.
        # Example: "A_SRR26842280_B" -> "SRR26842280"
        if filename.startswith("A_") and filename.endswith(
            "_B.fasta"
        ):
            sample_id = (
                filename.replace("A_", "")
                .replace("_B.fasta", "")
            )
        else:
            # Fallback: use filename without extension
            sample_id = os.path.splitext(filename)[0]
            print(
                f"Warning: Filename '{filename}' does not match expected format. "
                f"Using '{sample_id}' as sample ID."
            )

        # Output folder and file path
        sample_output_dir = os.path.join(args.output, sample_id)
        final_output_path = os.path.join(sample_output_dir, "tar1_blocks.fa")

        # Run the full pipeline
        if run_pipeline(args.hmm, fasta_path, final_output_path, args.cpu, args.incE):
            success_count += 1
        else:
            fail_count += 1
            print(f"Pipeline failed for {filename}.\n", file=sys.stderr)

    print("=" * 40)
    print("All tasks completed.")
    print(f"Total files processed: {len(fasta_files)}")
    print(f" Successful: {success_count}")
    print(f" Failed: {fail_count}")
    print("=" * 40)


if __name__ == "__main__":
    main()
