# TAR1 manuscript figure code (Figures 1–4)

This repository contains the figure-generation code used for the TAR1 manuscript.

The current repository supports:

- **Figure 1**: representative TAR1 loci, CWT scalograms, and motif profiles
- **Figure 2**: subtelomeric/centromeric/whole-chromosome localization analyses
- **Figure 3**: AT/CG ratio profiles and chimpanzee cap-association analysis
- **Figure 4**: 29 bp motif-defined positional analyses, delta-a correlation plots, and V-mode profiles

The code is organized figure-by-figure so that each main figure can be regenerated with a small number of commands from a fixed directory structure.

---

## Data availability 
All third-party sequencing datasets analyzed in this study are available from the original sources cited in the manuscript. The pre-calculated TAR1 block package required to reproduce the TAR1 block-based analyses has been deposited at Zenodo and is publicly available at https://doi.org/10.5281/zenodo.18617725.

### Download and extract pre-calculated TAR1 blocks

To reproduce the figure generation, please download the pre-calculated TAR1 block package from the Zenodo link above and extract its contents into the `data/tar1_blocks/` directory.

After extraction, the directory structure must look exactly like this:

```text
data/tar1_blocks/
├── human_CHM13/tar1_blocks.fa
├── chimp/tar1_blocks.fa
├── bonobo/tar1_blocks.fa
├── gorilla/tar1_blocks.fa
├── borang/tar1_blocks.fa
├── sorang/tar1_blocks.fa
├── hg002/tar1_blocks.fa
├── JH*/tar1_blocks.fa
└── SRR*/tar1_blocks.fa
```

## 1. Repository layout

```text
figure1_repo/
├── README.md
├── requirements.txt
├── data/
│   ├── annotation/
│   │   └── chm13v2.0_merged_centromeres.bed
│   ├── csv/
│   │   ├── delta_a_summary_SRR_only.csv
│   │   ├── schmidt_corres.csv
│   │   └── table_s2_donor_telomere_lengths.csv
│   ├── processed/
│   │   └── figure1_loci/
│   ├── reference/
│   │   ├── chm13v2.0.fa
│   │   ├── GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna
│   │   ├── GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna
│   │   ├── GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna
│   │   ├── GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna
│   │   └── GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna
│   ├── repeatmasker_out/
│   │   └── chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out
│   └── tar1_blocks/
│       ├── human_CHM13/tar1_blocks.fa
│       ├── chimp/tar1_blocks.fa
│       ├── bonobo/tar1_blocks.fa
│       ├── gorilla/tar1_blocks.fa
│       ├── borang/tar1_blocks.fa
│       ├── sorang/tar1_blocks.fa
│       ├── hg002/tar1_blocks.fa
│       ├── JH*/tar1_blocks.fa
│       └── SRR*/tar1_blocks.fa
├── metadata/
│   ├── figure1_manifest.tsv
│   ├── figure2_manifest.tsv
│   ├── figure2E_manifest.tsv
│   └── figure3A_manifest.tsv
├── scripts/
│   ├── figure1_common.py
│   ├── figure1_prepare_loci.py
│   ├── figure1_cwt.py
│   ├── figure1_motif.py
│   ├── figure1_run.py
│   ├── figure2_common.py
│   ├── figure2_prepare_windows.py
│   ├── figure2_scalogram.py
│   ├── figure2_distribution.py
│   ├── figure2_density.py
│   ├── figure2D_centromere.py
│   ├── figure2E_whole_chromosome.py
│   ├── figure2_run.py
│   ├── figure3_common.py
│   ├── figure3A_atcg_profile.py
│   ├── figure3B_cap_association.py
│   ├── figure3_run.py
│   ├── figure4_common.py
│   ├── figure4_positions.py
│   ├── figure4_correlations.py
│   ├── figure4_vmode.py
│   └── figure4_run.py
└── results/
    ├── figure1/
    ├── figure2/
    ├── figure3/
    └── figure4/
```

---

## 2. Environment setup

Recommended Python version:
- Python 3.9

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## 3. Required input files

### 3.1 Reference genomes

Place the following files under `data/reference/`:
- `chm13v2.0.fa`
- `GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna`
- `GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna`
- `GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna`
- `GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna`
- `GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna`

### 3.2 RepeatMasker annotation

Place the following file under `data/repeatmasker_out/`:
- `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out`

### 3.3 Centromere annotation

Place the following file under `data/annotation/`:
- `chm13v2.0_merged_centromeres.bed`

### 3.4 TAR1 block FASTA files

Place species/sample-specific TAR1 FASTA files under `data/tar1_blocks/`.

Examples:
- `data/tar1_blocks/human_CHM13/tar1_blocks.fa`
- `data/tar1_blocks/chimp/tar1_blocks.fa`
- `data/tar1_blocks/gorilla/tar1_blocks.fa`
- `data/tar1_blocks/hg002/tar1_blocks.fa`
- `data/tar1_blocks/JH104.F48.NB70/tar1_blocks.fa`
- `data/tar1_blocks/SRR26842280/tar1_blocks.fa`

### 3.5 CSV files used in Figure 4

Place the following files under `data/csv/`:
- `delta_a_summary_SRR_only.csv`
- `schmidt_corres.csv`
- `table_s2_donor_telomere_lengths.csv`

---

## 4. Download large reference files

Large reference FASTA files and the CHM13 RepeatMasker annotation are not tracked in the repository.
Download them directly into the expected folders before running the figure scripts.

### 4.1 Reference genomes

Run the following commands from `data/reference/`:

```bash
cd /Users/daecheolmoon/Downloads/figure1_repo/data/reference

curl -L -O [https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)
curl -L -O [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz)
curl -L -O [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz)
curl -L -O [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/625/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/625/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz)
curl -L -O [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/655/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/655/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz)
curl -L -O [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz)
```

Then unzip them:

```bash
gunzip -f chm13v2.0.fa.gz
gunzip -f GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz
gunzip -f GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz
gunzip -f GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz
gunzip -f GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz
gunzip -f GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz
```

Optional check:

```bash
ls -lh /Users/daecheolmoon/Downloads/figure1_repo/data/reference
```

### 4.2 RepeatMasker output

Run the following command from `data/repeatmasker_out/`:

```bash
cd /Users/daecheolmoon/Downloads/figure1_repo/data/repeatmasker_out

curl -L -O [https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out)
```

Optional check:

```bash
ls -lh /Users/daecheolmoon/Downloads/figure1_repo/data/repeatmasker_out
```

---

## 5. Metadata files

The repository uses fixed manifest tables for representative panels.

- `metadata/figure1_manifest.tsv`: Representative loci for Figure 1.
- `metadata/figure2_manifest.tsv`: Representative subtelomeric windows for Figure 2A.
- `metadata/figure2E_manifest.tsv`: Whole-chromosome examples for Figure 2E.
- `metadata/figure3A_manifest.tsv`: Representative windows for Figure 3A.

These files are used directly by the corresponding run scripts.

---

## 6. Running each figure

All commands below assume that you are in the repository root:

```bash
cd /Users/daecheolmoon/Downloads/figure1_repo
```

---

## 7. Figure 1

### Purpose
Figure 1 generates:
- representative TAR1 FASTA panels
- CWT scalograms
- motif profile CSVs and plots

### Run
```bash
python scripts/figure1_run.py \
  --manifest metadata/figure1_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure1
```

### Main outputs
- `results/figure1/fasta/`
- `results/figure1/cwt/`
- `results/figure1/motif/profile_csv/`
- `results/figure1/motif/profile_plots/`
- `results/figure1/motif/summary/`

### Files to check
- `results/figure1/cwt/figure1_cwt_summary.csv`
- representative FASTA files in `results/figure1/fasta/`
- motif summary in `results/figure1/motif/summary/figure1_motif_profile_summary.csv`

---

## 8. Figure 2

### Purpose
Figure 2 generates:
- subtelomeric CWT window examples
- cumulative TAR1 distribution plots
- normalized density plots
- centromere-focused human analysis
- whole-chromosome CWT examples

### Run
```bash
python scripts/figure2_run.py \
  --manifest metadata/figure2_manifest.tsv \
  --figure2e-manifest metadata/figure2E_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure2
```

### Main outputs
- `results/figure2/windows/`
- `results/figure2/figure2A/`
- `results/figure2/figure2B/`
- `results/figure2/figure2C/`
- `results/figure2/figure2D/`
- `results/figure2/figure2E/`

### Files to check
- `results/figure2/figure2A/figure2A_scalogram_summary.csv`
- `results/figure2/figure2B/tar1_counts_summary_all_species.csv`
- `results/figure2/figure2C/Figure2C_normalized_density_table.csv`
- `results/figure2/figure2D/Figure2D_summary.csv`
- `results/figure2/figure2E/Figure2E_summary.csv`

---

## 9. Figure 3

### Purpose
Figure 3 generates:
- AT/CG ratio profiles for representative windows
- chimpanzee cap-association analysis

### Run
```bash
python scripts/figure3_run.py \
  --figure3a-manifest metadata/figure3A_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure3
```

### Main outputs
- `results/figure3/figure3A/`
- `results/figure3/figure3B/`

### Files to check
- `results/figure3/figure3A/Figure3A_summary.csv`
- `results/figure3/figure3A/events/`
- `results/figure3/figure3A/plots/`
- `results/figure3/figure3B/Figure3B_chimp_cap_association.csv`

---

## 10. Figure 4

### Purpose
Figure 4 generates:
- motif-defined positional analyses
- ratio / age / TMM comparison plots
- telomere-length and delta-a correlation plots
- V-mode scale-loading plots

### Run
```bash
python scripts/figure4_run.py \
  --data-root ./data \
  --base-outdir results/figure4
```

### Main outputs
- `results/figure4/figure4_positions/`
- `results/figure4/figure4_correlations/`

### Files to check

**Position panels**
- `results/figure4/figure4_positions/Figure4_final_position_table.csv`
- `results/figure4/figure4_positions/panelB/Figure4B_cluster0_vs_cluster1.png`
- `results/figure4/figure4_positions/panelC/Figure4C_ratio_boxplot.png`
- `results/figure4/figure4_positions/panelD/Figure4D_TERT_vs_ALT_cluster0.png`
- `results/figure4/figure4_positions/panelE/Figure4E_Cluster_0_Pos_Median.png`
- `results/figure4/figure4_positions/panelF/Figure4F_Cluster_1_Pos_Median.png`
- `results/figure4/figure4_positions/panelK/Figure4K_cluster_density_profiles.png`

**Correlation panels**
- `results/figure4/figure4_correlations/Figure4A_total_and_cancer_correlation.png`
- `results/figure4/figure4_correlations/Figure4A_total_human_tar1_vs_telomere.csv`
- `results/figure4/figure4_correlations/Figure4G_H_correlations.png`
- `results/figure4/figure4_correlations/Figure4GH_summary.csv`

---

## 11. Recommended execution order

If you want to regenerate everything from scratch:

```bash
python scripts/figure1_run.py \
  --manifest metadata/figure1_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure1

python scripts/figure2_run.py \
  --manifest metadata/figure2_manifest.tsv \
  --figure2e-manifest metadata/figure2E_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure2

python scripts/figure3_run.py \
  --figure3a-manifest metadata/figure3A_manifest.tsv \
  --data-root ./data \
  --base-outdir results/figure3

python scripts/figure4_run.py \
  --data-root ./data \
  --base-outdir results/figure4
```

---

## 12. Notes
- The scripts assume the repository root layout shown above.
- The default data root is `./data`.
- Figure panels are driven by fixed manifest files where applicable.
- Large reference FASTA files are not bundled and must be downloaded separately.
- Figure 4 uses the CSV files in `data/csv/` and the pre-existing TAR1 block FASTA structure in `data/tar1_blocks/`.

---

## 13. requirements.txt

A minimal working `requirements.txt` for the current repository is:

```text
numpy
pandas
matplotlib
scipy
biopython
pywavelets
tqdm
```

If your environment already contains these packages, no further installation is required.

---

## 14. Output summary

After successful execution, the main figure results can be checked at:
- `results/figure1/`
- `results/figure2/`
- `results/figure3/`
- `results/figure4/`

For manuscript assembly, the most directly useful files are:
- panel image files (`.png`)
- summary tables (`.csv`)
- representative FASTA files for documented loci/windows