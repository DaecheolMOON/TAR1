#!/usr/bin/env python3
import argparse
from Bio import SeqIO


def extract_chrom_ends(input_fasta, output_fasta, flank_length=30000):
    """
    Extract the first and last `flank_length` bp from each chromosome record
    and write them to a new FASTA file with headers chrXp / chrXq.

    Chromosome label:
        - If the word "chromosome" appears in the record description, the text
          after "chromosome" up to the next comma is used (e.g. "1", "X").
        - Otherwise, the record ID is used.
    """
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            desc = record.description

            # Try to extract chromosome label from the description
            chrom_num = None
            for token in desc.split():
                if token.lower() == "chromosome":
                    chrom_num = desc.split("chromosome")[-1].split(",")[0].strip()
                    break
            if chrom_num is None:
                # Fallback: use record ID
                chrom_num = record.id

            seq = str(record.seq)
            length = len(seq)

            if length < flank_length * 2:
                print(
                    f"Warning: Chromosome {chrom_num} shorter than "
                    f"2x flank length ({length} bp)"
                )
                continue

            # p arm: start
            p_seq = seq[:flank_length]
            # q arm: end
            q_seq = seq[-flank_length:]

            # Write to FASTA
            out_f.write(f">chr{chrom_num}p\n{p_seq}\n")
            out_f.write(f">chr{chrom_num}q\n{q_seq}\n")

            print(f"Processed chr{chrom_num}: total length={length} bp")

    print(f"\nFinished! Output saved to: {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract first and last N bp of each chromosome into a new FASTA file."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input reference genome FASTA file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA file",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=30000,
        help="Flank length in bp (default: 30000)",
    )

    args = parser.parse_args()
    extract_chrom_ends(args.input, args.output, args.flank)
