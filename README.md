# TAR1 Subtelomeric Periodicity Analysis Pipeline

This repository contains a modular Python pipeline for:

- extracting telomeric/subtelomeric regions from reference genomes and assemblies,
- detecting TAR1-like blocks using profile HMMs,
- performing k-mer based continuous wavelet transform (CWT) analysis of repeat structure,
- decomposing CWT energy into modes using singular value decomposition (SVD),
- deriving simplified scalar metrics (e.g. delta_a) for downstream statistical analysis.

The code is organized into functional groups:

- **Group 1** – Data preprocessing and TAR1 extraction  
- **Group 2** – CWT analysis and periodicity detection  
- **Group 3** – CWT–SVD mode analysis and delta_a computation  

---

## 1. Installation

### 1.1 Python environment

Create and activate a Python 3 virtual environment (optional but recommended), then install the dependencies:

```bash
pip install -r requirements.txt
````

`requirements.txt`:

```txt
numpy
pandas
PyWavelets
scipy
biopython
matplotlib
```

### 1.2 External tools

Some scripts depend on external command-line tools:

* [`nhmmer`](http://hmmer.org/) (from the HMMER suite)
* [`bedtools`](https://bedtools.readthedocs.io/en/latest/) (for `getfasta`)

Example installation with conda:

```bash
conda install -c bioconda hmmer bedtools
```

Example installation with Homebrew (macOS):

```bash
brew install hmmer
brew install bedtools
```

Make sure these tools are in your `PATH` before running `Y001_extract_tar1_sequence.py`.

---

## 2. Suggested directory layout

The repository can be structured as follows:

```text
TAR1-pipeline/
  README.md
  requirements.txt

  group1_preprocessing/
    extract_chrom_ends.py
    Y001_extract_tar1_sequence.py
    parse_tar1_blocks.py

  group2_cwt/
    cwt_batch_from_fasta.py
    cwt_strong_regions_from_fasta.py
    weak_ranges_from_fasta.py
    aggregate_signal_ranges.py
    V004_fasta_to_scalogram_ver2.py

  group3_modes/
    V001_draw_v_mode.py
    V001_draw_v_mode_primate.py
    compute_delta_a.py

  data/                # (optional) example input files
  results/             # (optional) output directory for figures/CSVs
```

You are free to adjust the directory layout, but the scripts assume relative paths for input/output unless otherwise specified.

---

## 3. Group 1 – Data Preprocessing and TAR1 Extraction

Group 1 prepares the basic sequence data used in downstream CWT and TAR1 analyses.

There are two main branches:

1. **Reference-based chromosome end extraction** → used as input to CWT analysis on subtelomeric regions.
2. **TAR1 HMM-based block extraction and summarization** → used to locate and summarize TAR1 blocks across samples.

### 3.1 Extract chromosome ends from a reference genome

**Script**

* `extract_chrom_ends.py`

**Description**

Extracts the first and last *N* base pairs from each chromosome in a reference FASTA and writes them as separate records:

* `chrXp` – p-arm end
* `chrXq` – q-arm end

These records are used as input for CWT analysis on subtelomeric regions.

**Inputs**

* Reference genome FASTA file (e.g. `chm13v2.0.fa`, `hg38.fa`)

**Outputs**

* FASTA file containing chromosome-arm ends, e.g.:

  * `chm13v2.0_chrom_ends_30000.fa`
  * Headers: `>chr1p`, `>chr1q`, `>chr20q`, …

**Usage**

```bash
python extract_chrom_ends.py \
    --input chm13v2.0.fa \
    --output chm13v2.0_chrom_ends_30000.fa \
    --flank 30000
```

**Dependencies**

* Python 3
* `biopython`

---

### 3.2 Extract TAR1 blocks using a profile HMM

**Script**

* `Y001_extract_tar1_sequence.py`

**Description**

Runs a `nhmmer → BED → bedtools getfasta` pipeline on multiple input FASTA files to extract TAR1-like blocks:

1. Search each FASTA file with `nhmmer` using a TAR1 profile HMM.
2. Parse the `.tbl` output and convert hits into BED intervals.
3. Use `bedtools getfasta` to extract the corresponding sequences.
4. Save TAR1 block FASTA files under per-sample subdirectories.

**Inputs**

* `--hmm`: profile HMM file (e.g. `TAR1.hmm`)
* `--fasta`: directory containing input FASTA files

  * Default filename convention:

    ```text
    A_<SAMPLE_ID>_B.fasta
    ```

    If a file matches this pattern, `<SAMPLE_ID>` is extracted automatically.
    Otherwise, the basename without extension is used as sample ID.
* `--output`: base output directory

**Output structure**

For each input FASTA `<file>.fasta`:

```text
<output>/<sample_id>/tar1_blocks.fa
```

* `<sample_id>` is either `<SAMPLE_ID>` from the filename pattern or the basename without extension.
* `tar1_blocks.fa` contains all TAR1-like blocks detected by `nhmmer`.

**Usage**

```bash
python Y001_extract_tar1_sequence.py \
    --hmm TAR1.hmm \
    --fasta /path/to/fasta_dir \
    --output /path/to/output_dir \
    --cpu 8 \
    --incE 1e-5
```

**Dependencies**

* Python 3 (standard library only)
* External tools: `nhmmer`, `bedtools`

---

### 3.3 Summarize TAR1 block coordinates from FASTA headers

**Script**

* `parse_tar1_blocks.py`

**Description**

Parses TAR1 block FASTA headers of the form:

* `>chr19q:25658-27061`
* or `>::chr19q:25658-27061`

For each chromosome arm:

* merges multiple block ranges into a single `(min_pos, max_pos)` interval,
* computes `span = max_pos - min_pos`,
* and writes a summary table.

**Inputs**

* `--fasta`: TAR1 block FASTA file

  * Headers must match:

    ```text
    >chr<chrom>[p|q]:<start>-<end>
    ```

* `--max-pos` (optional):

  * If given, entries with `start > max_pos` or `end > max_pos` are excluded.

**Outputs**

* CSV summary:

  * Columns: `arm, min_pos, max_pos, span`

**Usage**

```bash
python parse_tar1_blocks.py \
    --fasta tar1_blocks.fa \
    --out tar1_block_ranges.csv

# Optional: restrict to blocks with coordinates <= 100,000
python parse_tar1_blocks.py \
    --fasta tar1_blocks.fa \
    --out tar1_block_ranges_max100k.csv \
    --max-pos 100000
```

**Dependencies**

* Python 3 (standard library only)

---

### 3.4 Data flow in Group 1

1. **Reference-based subtelomeric analysis**

```text
reference.fa
   └─ extract_chrom_ends.py
        └─ reference_chrom_ends.fa (chrXp/chrXq)
             └─ used as input to Group 2 (CWT)
```

2. **TAR1 HMM-based detection and summarization**

```text
TAR1.hmm + sample FASTA directory
   └─ Y001_extract_tar1_sequence.py
        └─ <output>/<sample_id>/tar1_blocks.fa
             └─ (optional header normalization)
                  └─ parse_tar1_blocks.py
                       └─ tar1_block_ranges.csv
```

---

## 4. Group 2 – CWT Analysis and Periodicity Detection

Group 2 performs k-mer based CWT analysis on sequences and detects strong/weak periodic regions.
There are four main usage patterns:

1. Reference-arm CWT analysis and run-level summary
2. Read-level strong-run detection (TSV)
3. Read-level weak-span detection (CSV)
4. Scalogram + SVD mode visualization

### 4.1 Batch CWT analysis on chromosome-arm FASTA

**Script**

* `cwt_batch_from_fasta.py`

**Description**

For each FASTA record (typically `chrXp`, `chrXq` or similar):

1. Build a k-mer repetitiveness signal:
   `signal[i] = frequency of k-mer seq[i:i+k] in the entire sequence`
2. Compute a CWT scalogram (`pywt.cwt`, Morlet wavelet) over a range of scales.
3. Select an “optimal” scale by average peak strength.
4. At the optimal scale, detect strong peaks and weak troughs:

   * filter peaks by prominence and width,
   * keep the top fraction of strongest peaks,
   * group contiguous peak/trough indices into runs.
5. For each run:

   * compute mean AT/CG ratio in a local window,
   * extract a short sample sequence around a representative index,
   * convert k-mer indices to base coordinates.

Saves a scalogram PNG per record and a global CSV summarizing all strong/weak runs.

**Inputs**

* `--fasta`: FASTA file with records (e.g. `chr1p`, `chr1q`, or generic contigs)
* Optional parameters (may vary slightly depending on version):

  * `--k`: k-mer length (default: 20)
  * `--scale_min`, `--scale_max`: CWT scale range (e.g. 200–401)
  * `--peak_keep_frac`: fraction of strongest peaks to retain (e.g. 0.5)
  * `--prominence_frac`: prominence fraction for peak detection (e.g. 0.10)
  * `--width`: minimum peak width (e.g. 5)
  * `--log_each`: if set, write a small log file per record

**Outputs**

Under `--outdir`:

* `<record_id>_scalogram.png` for each FASTA record
* `summary_runs.csv` containing run-level information, including:

  * `record_id`
  * `signal_type` (`"strong"` or `"weak"`)
  * `start_kmer_idx`, `end_kmer_idx`
  * `start_base`, `end_base`
  * `run_kmer_len`
  * `scale_used`
  * `mean_atcg`
  * `sample_seq_start`, `sample_seq_end`, `sample_seq_len`
  * `sample_seq`

**Typical usage (reference arms)**

```bash
# 1) Prepare chromosome-arm FASTA from a reference
python extract_chrom_ends.py \
    --input chm13v2.0.fa \
    --output chm13v2.0_chrom_ends_30000.fa \
    --flank 30000

# 2) Run CWT analysis and strong/weak run detection
python cwt_batch_from_fasta.py \
    --fasta chm13v2.0_chrom_ends_30000.fa \
    --outdir cwt_results_chm13_30kb \
    --k 20 \
    --scale_min 200 \
    --scale_max 401
```

**Dependencies**

* Python 3
* `numpy`, `pandas`, `PyWavelets`, `scipy`, `biopython`, `matplotlib`

---

### 4.2 Strong-run detection per record (TSV summarizer)

**Script**

* `cwt_strong_regions_from_fasta.py`

**Description**

Runs a similar k-mer CWT pipeline as `cwt_batch_from_fasta.py`, but focuses only on **strong (peak) runs** and outputs a compact TSV file:

* detects strong peaks at an optimal CWT scale,
* groups contiguous peaks into runs,
* computes AT/CG ratio and a sample sequence for each run,
* writes one row per strong run.

This is useful for read-level or contig-level analysis where a simple text summary is preferred over images.

**Inputs**

* `--fasta`: FASTA file (any set of reads/contigs/arms)
* `--outdir`: output directory (default: `./TXT`)
* `--outfile`: output TSV filename (default: `strong_regions_results.tsv`)
* Other parameters control `k`, scale range, prominence, etc.

**Outputs**

* `<outdir>/<outfile>` (TSV), with columns:

  * `read_id`
  * `seq_len`
  * `k`
  * `optimal_scale`
  * `run_idx`
  * `kmer_start`, `kmer_end`
  * `rep_kmer_idx`
  * `sample_seq_start`, `sample_seq_end`
  * `mean_AT_over_CG_window`
  * `sample_seq`

**Dependencies**

* Python 3
* `numpy`, `PyWavelets`, `scipy`

---

### 4.3 Weak-signal span detection per record (CSV summarizer)

**Script**

* `weak_ranges_from_fasta.py`

**Description**

Runs k-mer based CWT and focuses on **weak (trough) regions**:

* computes a repetitiveness signal and CWT as above,
* detects troughs (peaks on the negated signal),
* groups contiguous trough indices into runs,
* summarizes all weak runs per record by a single global span:

  * `min_kmer`: minimum k-mer index over all weak runs
  * `max_kmer`: maximum k-mer index over all weak runs

This provides a simple measure of how widely weak regions are spread along each sequence.

**Inputs**

* `--fasta`: input FASTA file (reads, contigs, or reference arms)
* `--outdir`: output directory (default: `./CSV`)
* `--outfile`: output CSV filename (default: `weak_ranges.csv`)

**Outputs**

* `<outdir>/<outfile>` (CSV) with columns:

  * `read_id`
  * `min_kmer`
  * `max_kmer`

Each row summarizes the global weak span of one FASTA record.

**Dependencies**

* Python 3
* `numpy`, `PyWavelets`, `scipy`

---

### 4.4 Aggregate CWT run ranges per arm

**Script**

* `aggregate_signal_ranges.py`

**Description**

Aggregates CWT run-level results (e.g. from `summary_runs.csv`) into **arm-level ranges**. For each record/arm:

1. Optionally filter by `signal_type` (`strong`, `weak`, or `both`).

2. Compute:

   * `min_kmer` over all runs
   * `max_kmer` over all runs
   * `n_segments` = number of runs

3. Optionally adjust orientation for q-arms:

   * for arm names ending with `q`, k-mer indices are mirrored within a fixed window size (`--window-bp`), so that p- and q-arms share a comparable coordinate system.

4. Compute:

   * `span = max - min`
   * `span_per_segment = span / n_segments`

5. Keep only arms with `span_per_segment <= max_span_per_segment`.

**Inputs**

* `--in`: input CSV file (e.g. `summary_runs.csv` from `cwt_batch_from_fasta.py`)
* `--out`: output CSV file
* `--signal`: which signal types to use (`strong`, `weak`, or `both`)
* `--window-bp`: window size (in bp) to use when mirroring q-arms
* `--max-span-per-segment`: threshold for filtering arms
* Column names for record, signal type, start/end indices are automatically detected,
  but can be overridden with `--record-col`, `--signal-col`, `--start-col`, `--end-col`.

**Outputs**

* CSV with columns (at minimum):

  * `arm`
  * `min`
  * `max`
  * `span`
  * `n_segments`

**Typical usage**

```bash
python aggregate_signal_ranges.py \
    --in cwt_results_chm13_30kb/summary_runs.csv \
    --out cwt_results_chm13_30kb/agg_strong_arms.csv \
    --signal strong \
    --window-bp 30000 \
    --max-span-per-segment 2000
```

**Dependencies**

* Python 3
* `pandas`

---

### 4.5 CWT scalogram and SVD mode visualization

**Script**

* `V004_fasta_to_scalogram_ver2.py`

**Description**

For each FASTA record:

1. Build a k-mer repetitiveness signal (k=20).

2. Compute a CWT scalogram over scales 100–400 (Morlet wavelet).

3. Compute the SVD of the absolute CWT matrix:

   * `F = |coeffs|` (shape: `n_scales × n_positions`)
   * `F = U S Vᵀ`

4. Visualize:

   * **Scalogram**: `|coeffs|` vs (position, scale)
   * **SVD basis (scale modes)**: first three columns of `U` vs scale
   * **SVD components (position modes)**: first three rows of `Vᵀ` vs position

This script is focused on visualization (no CSV output), and is useful for qualitative inspection of periodic structures and mode patterns.

**Inputs**

* `--input_fasta`: input FASTA file
* `--output_dir`: directory to save PNG plots

**Outputs**

For each record `id`:

* `<output_dir>/<id>_scalogram.png`
* `<output_dir>/<id>_svd_basis.png`
* `<output_dir>/<id>_svd_components.png`

**Usage**

```bash
python V004_fasta_to_scalogram_ver2.py \
    --input_fasta chm13v2.0_chrom_ends_30000.fa \
    --output_dir svd_plots_chm13_30kb
```

**Dependencies**

* Python 3
* `numpy`, `PyWavelets`, `biopython`, `matplotlib`

---

## 5. Group 3 – CWT–SVD Mode Analysis and Delta-a Computation

Group 3 takes the CWT energy information and performs SVD-based mode analysis and scalar summarization.

### 5.1 CWT energy SVD across human arms and reads

**Script**

* `V001_draw_v_mode.py`

**Description**

Builds CWT energy matrices and applies SVD to compare periodicity modes:

1. **CHM13 reference (by chromosome arm)**
2. **HG002m (by chromosome arm)**
3. **HG002_karimian (by read)**

For each matrix (rows = arms or reads, columns = scales):

* `M = CWT_energy(arm or read)` → shape `(n_arms or n_reads, n_scales)`
* `M = U S Vᵀ` via SVD

  * **U modes**: contribution of each arm/read to each mode
  * **V modes**: scale-dependent periodicity modes

Additionally, computes cosine similarities between V modes of HG002m and HG002_karimian.

**Inputs**

* `--ref`: CHM13 FASTA (chromosome-arm sequences)

  * headers must contain `chr<chrom>[p|q]` (e.g. `chr20q`), which are used as arm IDs.
* `--hg002`: HG002m FASTA (chromosome-arm sequences, same arm naming convention)
* `--kar`: HG002_karimian FASTA (read-level sequences)
* `--output_dir`: directory to store plots

Each FASTA is parsed as:

* `parse_fasta_to_arms` for `--ref`, `--hg002`

  * sequences with length ≥ `MIN_LEN + KMER_K - 1` are retained
  * multiple records for the same arm are concatenated
* `parse_fasta_to_reads` for `--kar`

  * read IDs are derived from headers (optionally using `::` as a delimiter)

**Outputs**

* `<output_dir>/CHM13_U_modes.png`
* `<output_dir>/CHM13_V_modes.png`
* `<output_dir>/HG002m_U_modes.png`
* `<output_dir>/HG002m_V_modes.png`
* `<output_dir>/HG002_karimian_V_modes.png`
* Console output with cosine similarities between V modes of HG002m and Karimian:

  * `Mode1: ...`
  * `Mode2: ...`
  * `Mode3: ...`

**Usage**

```bash
python V001_draw_v_mode.py \
    --ref chm13_subtelomere_arms.fa \
    --hg002 hg002m_subtelomere_arms.fa \
    --kar hg002_karimian_reads.fa \
    --output_dir vmode_plots
```

**Dependencies**

* Python 3
* `numpy`, `PyWavelets`, `scipy`, `matplotlib`

---

### 5.2 Primate TAR1 periodic mode analysis

**Script**

* `V001_draw_v_mode_primate.py`

**Description**

For each primate species (Bornean orangutan, bonobo, chimpanzee), this script:

* parses a TAR1-block FASTA file into individual sequence blocks,
* computes a k-mer based CWT energy vector per block (energy vs scale),
* stacks these vectors into a matrix (blocks × scales),
* performs SVD on this matrix,
* and visualizes:

  * **U modes**: contributions of each TAR1 block to each mode
  * **V modes**: scale-dependent periodicity patterns

**Inputs**

* `--borang`: Bornean orangutan TAR1 FASTA (one or more blocks)
* `--bonobo`: bonobo TAR1 FASTA
* `--chimp`: chimpanzee TAR1 FASTA
* `--output_dir`: directory to store PNG plots

Each FASTA header should contain a unique block identifier; if `'::'` is present, the substring after `'::'` is used as block ID (e.g. `::chr8:13587-15321` → `chr8:13587-15321`).

**Outputs**

For each species `<species>` in `{borang, bonobo, chimp}`:

* `<output_dir>/Primate_<species>_U_modes.png`
* `<output_dir>/Primate_<species>_V_modes.png`

**Usage**

```bash
python V001_draw_v_mode_primate.py \
    --borang orangutan_TAR1.fa \
    --bonobo bonobo_TAR1.fa \
    --chimp chimp_TAR1.fa \
    --output_dir primate_vmode_plots
```

**Dependencies**

* Python 3
* `numpy`, `PyWavelets`, `scipy`, `matplotlib`

---
### 5.3 Read-level Mode 2 peak detection (adaptive scale)

Script  
- `Y002_mode2_scale_scan.py`  *(originally `Y002_450_scale_ver2.py`)*

Description  
- For each sample, this script:
  - takes a TAR1 block FASTA (`tar1_blocks.fa`) containing multiple reads,
  - builds a k-mer-based CWT energy vector per read (energy vs scale),
  - stacks these vectors into a matrix `M` (reads × scales),
  - performs SVD on `M`,
  - extracts the **Mode 2 scale-loading vector**,
  - adaptively extends the scale range until at least two Mode 2 peaks are detected
    (or until a hard ceiling, e.g. 800 bp),
  - saves the Mode 2 peak positions and a scale-loading plot.

The Mode 2 peak positions are later used by `compute_delta_a.py` to compute
the scalar metric `delta_a` per sample.

Inputs  
- `--base_dir`: base directory containing sample subdirectories.
- Options:
  - `--sample-prefix`: optional prefix of sample directory names (default: process all subdirectories).
  - `--fasta-name`: FASTA filename inside each sample directory (default: `tar1_blocks.fa`).
  - `--output-subdir`: subdirectory name where summary and plot will be saved (default: `Y002_430`).

Expected directory structure  
```text
base_dir/
  SAMPLE_001/
    tar1_blocks.fa
  SAMPLE_002/
    tar1_blocks.fa
  ...
```

For each sample directory `<sample>` that matches the prefix filter,
the script creates:

```text
base_dir/<sample>/Y002_430/
  mode2_peak_summary.txt
  mode2_scale_loading.png
```

Outputs

* `mode2_peak_summary.txt`: Mode 2 peak positions in bp:

  * `Mode 2 peaks (bp): p1, p2, ...`
* `mode2_scale_loading.png`: Mode 2 scale-loading plot with peaks marked.

Usage

```bash
python Y002_mode2_scale_scan.py \
    --base_dir /path/to/base_dir \
    --sample-prefix SAMPLE_ \
    --fasta-name tar1_blocks.fa \
    --output-subdir Y002_430
```


---

### 5.4 Delta-a computation from Mode 2 peaks

**Script**

* `compute_delta_a.py`

**Description**

Collects Mode 2 peak positions from per-sample summary files and computes a scalar metric `delta_a` for each sample.

Each sample directory is expected to contain:

* `Y002_430/mode2_peak_summary.txt` with a line of the form:

  * `Mode 2 peaks (bp): p1, p2` (two peaks), or
  * `Mode 2 peaks (bp): p1` (single peak)

`delta_a` is defined as:

* two-peak case:

  * `delta_a = p2 - p1`
* single-peak case:

  * `delta_a = single_peak_max - p1` (default `single_peak_max = 800`)

This script aggregates `delta_a` values from all sample directories and writes them to a single CSV file.

**Inputs**

* positional `base_dir`: base directory containing sample subdirectories
* optional arguments:

  * `--sample-prefix`: prefix of sample directory names (default: `JH`)
  * `--single-peak-max`: reference peak position used in the single-peak case (default: `800` bp)
  * `-o/--output`: output CSV path (default: `<base_dir>/delta_a_summary.csv`)

**Expected directory structure**

```text
base_dir/
  JH001/
    Y002_430/
      mode2_peak_summary.txt
  JH002/
    Y002_430/
      mode2_peak_summary.txt
  ...
```

Only subdirectories whose names start with `--sample-prefix` are processed.

**Outputs**

* CSV with columns:

  * `sample`: sample directory name
  * `delta_a`: computed delta_a value

**Usage**

```bash
python compute_delta_a.py /path/to/base_dir \
    --sample-prefix JH \
    --single-peak-max 800 \
    --output /path/to/base_dir/delta_a_summary.csv
```

**Dependencies**

* Python 3 (standard library only)

---

## 6. Example end-to-end workflows

### 6.1 Reference-based subtelomeric periodicity (CHM13)

1. Extract chromosome-arm ends:

```bash
python extract_chrom_ends.py \
    --input chm13v2.0.fa \
    --output chm13v2.0_chrom_ends_30000.fa \
    --flank 30000
```

2. Run CWT analysis and detect strong/weak runs:

```bash
python cwt_batch_from_fasta.py \
    --fasta chm13v2.0_chrom_ends_30000.fa \
    --outdir cwt_results_chm13_30kb
```

3. Aggregate strong-run ranges per arm:

```bash
python aggregate_signal_ranges.py \
    --in cwt_results_chm13_30kb/summary_runs.csv \
    --out cwt_results_chm13_30kb/agg_strong_arms.csv \
    --signal strong \
    --window-bp 30000 \
    --max-span-per-segment 2000
```

4. Visualize CWT scalograms and SVD modes:

```bash
python V004_fasta_to_scalogram_ver2.py \
    --input_fasta chm13v2.0_chrom_ends_30000.fa \
    --output_dir svd_plots_chm13_30kb
```

---

### 6.2 Human sample mode comparison (CHM13 vs HG002m vs Karimian)

1. Prepare subtelomeric arm FASTA for CHM13 and HG002m (e.g. using `extract_chrom_ends.py`).
2. Prepare read-level FASTA for HG002_karimian.
3. Run V-mode comparison:

```bash
python V001_draw_v_mode.py \
    --ref chm13_subtelomere_arms.fa \
    --hg002 hg002m_subtelomere_arms.fa \
    --kar hg002_karimian_reads.fa \
    --output_dir vmode_plots
```

The script will generate U/V mode plots for CHM13 and HG002m, V modes for Karimian reads, and cosine similarities between HG002m and Karimian V modes.

---

### 6.3 Primate TAR1 mode comparison

1. Prepare TAR1 block FASTA for each primate (e.g. orangutan, bonobo, chimpanzee).
2. Run primate mode analysis:

```bash
python V001_draw_v_mode_primate.py \
    --borang orangutan_TAR1.fa \
    --bonobo bonobo_TAR1.fa \
    --chimp chimp_TAR1.fa \
    --output_dir primate_vmode_plots
```

This will generate U/V mode plots for each species, enabling qualitative comparison of TAR1 periodic structures across primates.

---

### 6.4 delta_a extraction for downstream statistics

Assuming each sample directory contains a TAR1 block FASTA (`tar1_blocks.fa`):

1) Run the read-level Mode 2 peak detection:

```bash
python Y002_mode2_scale_scan.py \
    --base_dir /path/to/base_dir \
    --sample-prefix JH \
    --fasta-name tar1_blocks.fa \
    --output-subdir Y002_430
````

This will create, for each sample directory `JH*`:

```text
/path/to/base_dir/JH*/Y002_430/
  mode2_peak_summary.txt
  mode2_scale_loading.png
```

2. Then compute `delta_a` from the Mode 2 peaks:

```bash
python compute_delta_a.py /path/to/base_dir \
    --sample-prefix JH \
    --single-peak-max 800 \
    --output /path/to/base_dir/delta_a_summary.csv
```

You can then use `delta_a_summary.csv` in downstream statistical analyses
(e.g. correlation with epigenetic marks, doubling time, ALT vs telomerase status, etc.).

````

---

If you need to extend the pipeline (e.g. add more species, change CWT scales, or integrate new summary metrics), the existing scripts are modular and can be adapted by adjusting parameters or adding new CLI options.


