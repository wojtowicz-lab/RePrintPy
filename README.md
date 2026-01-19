# reprint-tool

**Tool to compute DNA Repair FootPrint (RePrint) from COSMIC signatures and visualize the results as publication-ready PDF plots.**

![Tests](https://github.com/marcin119a/RePrint.py/actions/workflows/python-tests.yml/badge.svg)


## Features

- Computes RePrint (DNA repair footprint) from input signature matrices (e.g., COSMIC).
- Generates barplot visualizations for each signature.
- Exports all plots as separate PDF files and optionally PNG files.
- Generates and exports a similarity heatmap (dendrogram + optional heatmap) as PNG.
- Can save the computed RePrint matrix to a CSV/TSV file.
- Command-line interface for easy batch processing.

## Demo

Try the interactive demo on Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1XLEvSDd9vCTwiiHCA0Nk1ASAmoXm5m8A?usp=sharing)

## Installation


Install dependencies:

```bash
pip install numpy pandas  scipy
```

For development and testing:

```bash
pip install pytest pytest-cov
```

## Usage

### 1. Compute and plot RePrints, saving each signature as a separate PDF

```bash
python main.py --input path/to/input_signatures.tsv --output_dir pdfs --prefix reprint_ --sep '\t'
```

- `--input` – path to your input file (TSV/CSV, signatures in columns, mutation types in rows)
- `--output_dir` – directory to save PDF files (default: current directory)
- `--prefix` – prefix for output PDF files (default: `reprint_`)
- `--sep` – column separator in input file (default: tab)
- `--save_reprint` – (optional) path to save the computed RePrint matrix as CSV/TSV

### 2. Also export per-signature PNG files

```bash
python main.py \
  --input path/to/input_signatures.tsv \
  --output_dir pdfs \
  --prefix reprint_ \
  --sep '\t' \
  --export_png
```

This will save both PDFs and PNGs for each signature into `--output_dir`.

### 3. Generate and save similarity heatmap as PNG

You can generate a dendrogram-based similarity view across signatures and save it as PNG:

```bash
python main.py \
  --input path/to/input_signatures.tsv \
  --output_dir pdfs \
  --prefix reprint_ \
  --sep '\t' \
  --analyze_png pdfs/all_signatures.png \
  --analyze_metric rmse \
  --analyze_colorscale Blues \
  --analyze_method complete
```

Options:

- `--analyze_png` – output PNG path for the heatmap/dendrogram figure
- `--analyze_metric` – similarity metric: `rmse` (default) or `cosine`
- `--analyze_colorscale` – Plotly colorscale name (default: `Blues`)
- `--analyze_hide_heatmap` – include only dendrograms (no heatmap)
- `--analyze_method` – clustering linkage method (e.g., `complete`, `average`)

### 4. Save only the computed RePrint matrix (no plots)

```bash
python main.py --input path/to/input_signatures.tsv --save_reprint reprint_matrix.tsv --sep '\t'
```

### 5. Complete workflow: Generate all plots (signatures and RePrints)

The `generate_all_plots.py` script provides a complete workflow that:
- Generates both standard signature plots and RePrint plots
- Saves all data to CSV files
- Creates combined PDF files with all plots
- Generates individual PNG files for RePrints

```bash
python generate_all_plots.py <cosmic_file> [output_base_name]
```

**Example:**
```bash
python generate_all_plots.py tests/data/COSMIC_v3.4_SBS_GRCh37.txt cosmic_v3.4_GRCh37
```

## Input Format

Input should be a tab-separated (or CSV) file with mutation types as rows and signatures as columns, e.g.:

```
Type    Signature_1    Signature_2    ...
A[C>A]A 0.011          0.00068        ...
A[C>A]C 0.0091         0.00061        ...
...
```

## Output

- Individual PDF files for each signature (if using `--output_dir`)
- Or a single PDF file with all plots (if using `--single_pdf`)
- The computed RePrint matrix as a CSV/TSV file (if using `--save_reprint`)

## Example

```bash
python main.py --input tests/data/COSMIC_v2_SBS_GRCh37.txt --prefix reprint_ --sep '\t' --save_reprint output/reprint_matrix.tsv
```

## Project Structure

```
reprint_tool/
    core.py         # Core computation (RePrint)
    plot.py         # Plotting and PDF export (standard bar charts)
    plot_static.py  # RePrint plotting (seamless 3-panel visualizations)
    analyze.py      # Similarity heatmap generation
main.py             # Command-line interface
generate_all_plots.py  # Complete workflow script
tests/              # Unit tests and test data
```

