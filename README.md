# Telomere Sequence Analysis Pipeline

## Introduction

This project provides a comprehensive bioinformatics pipeline for processing telomere sequencing data. The included scripts are designed to analyze TAR1 region periodicity, related epigenetic features, telomere length, periodicity (delta a), and histone modification patterns, based on data from various foundational studies in the field.

## Data and Code Availability

The datasets and code supporting the conclusions of this article are available in the repository at [https://github.com/DaecheolMOON/TAR1](https://github.com/DaecheolMOON/TAR1).

## Citation

If you use this pipeline or its associated code in your research, please cite this GitHub repository as follows:

> Moon, Daecheol. (2025). *Telomere Sequence Analysis Pipeline*. GitHub. Retrieved from [https://github.com/DaecheolMOON/TAR1](https://github.com/DaecheolMOON/TAR1)

## Project Structure

```
TAR1/
├── README.md                # Project documentation
├── environment.yml          # Conda environment configuration
├── scripts/                 # All Python scripts
│   ├── Y001_extract_tar1_sequence.py
│   ├── Y002_450_scale_ver2.py
│   ├── V001_draw_v_mode.py
│   ├── V004_fasta_to_scalogram_ver2.py
│   ├── Q003_bw_fastq_correlation.py
│   └── compute_delta_a.py
├── data/                    # Essential data files (e.g., HMM profiles)
│   └── TAR1.hmm
├── example/                 # Example data for quick testing
│   ├── example_reads.fa
│   └── example_chipseq.bw
└── results/                 # Output directory for analysis results
    ├── tar_fasta_masked/
    ├── PLOT/
    └── summary_tables/
```

## Installation and Setup

#### 1\. Clone the Repository

```bash
git clone https://github.com/DaecheolMOON/TAR1.git
cd TAR1
```

#### 2\. Install External Dependencies

This pipeline requires `HMMER` and `BEDTools`. We recommend installing them via Conda/Bioconda. These tools are essential for sequence alignment and feature extraction as described in their respective publications (Eddy, 2011; Quinlan & Hall, 2010).

```bash
conda install -c bioconda hmmer bedtools
```

#### 3\. Create and Activate Conda Environment

The provided `environment.yml` file contains all necessary Python packages. Use it to create a dedicated Conda environment.

```bash
conda env create -f environment.yml
conda activate telomere_env
```

## Usage

All commands should be executed from the project's root directory (`TAR1/`).

### Step 1: Extract TAR1 Sequences

Extract TAR1 regions from input FASTA files using a Profile HMM.

```bash
python scripts/Y001_extract_tar1_sequence.py \
    --hmm data/TAR1.hmm \
    --fasta /path/to/your/input_fasta_dir/ \
    --output results/tar_fasta_masked
```

### Step 2: Periodicity Analysis and Visualization

#### U-mode and V-mode Analysis

Perform Singular Value Decomposition (SVD) to analyze periodicity patterns (V-modes) and per-sample contributions (U-modes).

```bash
python scripts/V001_draw_v_mode.py \
    --ref /path/to/chm13_tar1.fa \
    --hg002 /path/to/hg002_tar1.fa \
    --kar /path/to/karimian_tar1.fa \
    --output_dir results/PLOT/svd_modes
```

#### Single-Sequence Scalogram Visualization

Generate a scalogram for a single TAR1 sequence using Continuous Wavelet Transform (CWT) to visualize its periodic structure.

```bash
python scripts/V004_fasta_to_scalogram_ver2.py \
    --input_fasta results/tar_fasta_masked/some_sample/tar1_blocks.fa \
    --output_dir results/PLOT/scalograms
```

### Step 3: Detailed Analysis and Result Summarization

#### V-mode 2 Scale Analysis

Identify the two primary peaks in V-mode 2 using an adaptive scale approach.

```bash
python scripts/Y002_450_scale_ver2.py \
    --base_dir results/tar_fasta_masked
```

#### Calculate Delta A

Calculate the distance (delta\_a) between the two peaks identified by the `Y002` script and summarize the results in a CSV file.

```bash
python scripts/compute_delta_a.py results/tar_fasta_masked \
    --output results/summary_tables/delta_a_summary.csv
```

### Step 4: Correlation with Epigenomic Data

Analyze the correlation between the periodic signal of TAR1 sequences and ChIP-seq signals (e.g., H3K9me3) from bigWig files.

```bash
# NOTE: The script Q003_bw_fastq_correlation.py contains hardcoded paths.
# It is recommended to modify the script to accept file paths as command-line arguments.
# python scripts/Q003_bw_fastq_correlation.py
```

## Quick Test with Example Data

Test the main pipeline functionality using the small dataset provided in the `example/` directory.

```bash
# 1. Extract TAR1 sequences from the example FASTA file
python scripts/Y001_extract_tar1_sequence.py \
    --hmm data/TAR1.hmm \
    --fasta example/ \
    --output results/example_tar_masked

# 2. Visualize the extracted example TAR1 sequence
python scripts/V004_fasta_to_scalogram_ver2.py \
    --input_fasta results/example_tar_masked/example_reads/tar1_blocks.fa \
    --output_dir results/PLOT/example_scalogram
```

## References

**Primary Studies**

1.  Karimian, K., Groot, A., Huso, V., Kahidi, R., Tan, K.-T., Sholes, S., Keener, R., McDyer, J. F., Alder, J. K., Li, H., Rechtsteiner, A., & Greider, C. W. (2024). Human telomere length is chromosome end-specific and conserved across individuals. *Science*, *384*(6695), 533–539. DOI: `10.1126/science.ado0431`.
2.  Schmidt, T. T., Tyer, C., Rughani, P., Haggblom, C., Jones, J. R., Dai, X., Frazer, K. A., Gage, F. H., Juul, S., Hickey, S., & Karlseder, J. (2024). High resolution long-read telomere sequencing reveals dynamic mechanisms in aging and cancer. *Nature Communications*, *15*(1). DOI: `10.1038/s41467-024-48917-7`.
3.  Kaufman, D. S., & ATRX-mediated chromatin association of H3K9me3 at subtelomeres is required for maintaining telomere integrity in K562 cells. *Aging*, *8*(3), 481-492. DOI: `10.1080/15592294.2016.1169351`.

**Genomic Assemblies and Methods**
4\.  Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A., ... & Phillippy, A. M. (2022). The complete sequence of a human genome. *Science*, *376*(6588), 44-53. (T2T-CHM13v2.0 Assembly). DOI: `10.1126/science.abj6987`.
5\.  Zimin, A. V., Puiu, D., Luo, M. C., Al-Dhmour, H. M., & Salzberg, S. L. (2017). A new VGP reference genome of the naked mole-rat, *Heterocephalus glaber*. *Scientific data*, *4*(1), 1-5. (Referenced for T2T\_HG002v1.1 Assembly Context). DOI: `10.1038/sdata.2016.25`.
6\.  Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, *26*(6), 841-842. DOI: `10.1093/bioinformatics/btq033`.
7\.  Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS computational biology*, *7*(10), e1002195. DOI: `10.1371/journal.pcbi.1002195`.
8\.  Smit, A. F., Hubley, R., & Green, P. (2013-2015). *RepeatMasker Open-4.0*. [http://www.repeatmasker.org](https://www.google.com/search?q=http://www.repeatmasker.org).
9\.  Hubley, R., Finn, R. D., Clements, J., Eddy, S. R., Jones, T. A., Bao, W., ... & Coggill, P. (2016). The Dfam database of repetitive DNA families. *Nucleic acids research*, *44*(D1), D81-D89. DOI: `10.1093/nar/gkv1272`.
10\. Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). *Biological sequence analysis: Probabilistic models of proteins and nucleic acids*. Cambridge university press. (Theoretical basis for Profile HMMs). DOI: `10.1093/bioinformatics/14.9.755`.