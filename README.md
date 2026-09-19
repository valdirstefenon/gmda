# GMDA — Genetic and Morphometric Data Analysis Platform

**Version 4.4**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-see%20LICENSE-green.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#installation)

**GMDA** is a comprehensive desktop application for analyzing genetic and morphometric data, designed for researchers in population genetics, evolutionary biology, and related fields. It provides a user-friendly interface for complex statistical analyses including similarity calculations, population structure analysis, and correlation studies.

---

## 📋 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Input File Formats](#-input-file-formats)
  - [Genetic Data](#genetic-data-files)
  - [Morphometric Data](#morphometric-data-files)
- [Genetic Analysis](#-genetic-analysis)
  - [Codominant Markers](#codominant-markers)
  - [Dominant Markers (AFLP)](#dominant-markers-aflp)
  - [Haploid Markers (cpDNA, mtDNA, Y-chromosome)](#haploid-markers-cpdna-mtdna-y-chromosome)
- [Morphometric Analysis](#-morphometric-analysis)
- [Correlation Analysis](#-correlation-analysis)
- [Output Files](#-output-files)
- [Troubleshooting](#-troubleshooting)
- [Best Practices](#-best-practices)
- [Citation](#-citation)
- [Support](#-support)

---

## 🚀 Features

### Genetic Analysis

- **Similarity Calculations** — computes genetic similarity between query samples and reference databases.
- **UPGMA Dendrograms** — hierarchical clustering trees with optional bootstrap support.
- **DAPC Analysis** — Discriminant Analysis of Principal Components for population structure.
- **FST Outlier Detection** — identifies loci under selection using FST distribution.
- **Kinship Analysis** — estimates relatedness between individuals.
- **PCoA** — Principal Coordinate Analysis for visualizing genetic distances.

### Morphometric Analysis

- **Descriptive Statistics** — mean, median, mode with standard deviations.
- **ANOVA** — one-way ANOVA with Tukey HSD post-hoc tests.
- **Effect Sizes** — Cohen's d with 95% confidence intervals.
- **PCA** — Principal Component Analysis with variable loadings.
- **UPGMA Dendrograms** — hierarchical clustering of morphological traits.
- **Bootstrap Support** — optional bootstrap validation for dendrograms.

### Correlation Analysis

- **Mantel Tests** — correlation between genetic and morphometric distance matrices.
- **Dendrogram Comparisons** — Mantel tests between different dendrogram types.
- **PCoA/PCA Comparisons** — correlation between ordination spaces.
- **Statistical Significance** — permutation-based p-values.

---

## 💻 Installation

### Prerequisites

- **Python 3.8 or higher**
- Required Python packages:

```text
numpy>=1.21.0
pandas>=1.3.0
openpyxl>=3.0.9
scipy>=1.7.0
scikit-learn>=0.24.0
matplotlib>=3.4.0
seaborn>=0.11.0
statsmodels>=0.12.0
packaging>=21.0
```

### Install dependencies

```bash
pip install numpy>=1.21.0 pandas>=1.3.0 openpyxl>=3.0.9 scipy>=1.7.0 \
            scikit-learn>=0.24.0 matplotlib>=3.4.0 seaborn>=0.11.0 \
            statsmodels>=0.12.0 packaging>=21.0
```

### Run the application

```bash
python3 GMDA4.4.py
```

---

## ⚡ Quick Start

1. Prepare your input files in **Excel format** (`.xlsx` or `.xls`).
2. Launch GMDA.
3. Load **reference** and **query** datasets (genetic) and/or morphometric data.
4. Set parameters (replicates, bootstrap, population assignment).
5. Run analyses in sequence (see [Best Practices](#-best-practices)).
6. Retrieve results from automatically created output folders.

---

## 📁 Input File Formats

### Genetic Data Files

All genetic files must be in **Excel format** (`.xlsx` or `.xls`).

#### 1. Codominant Markers (microsatellites, SNPs)

| Sample | Pop  | Locus1_A1 | Locus1_A2 | Locus2_A1 | Locus2_A2 | Locus3_A1 | Locus3_A2 |
|--------|------|-----------|-----------|-----------|-----------|-----------|-----------|
| Ind1   | PopA | 150       | 150       | 200       | 204       | 300       | 302       |
| Ind2   | PopA | 150       | 152       | 200       | 200       | 300       | 300       |
| Ind3   | PopB | 148       | 150       | 202       | 204       | 300       | 304       |
| Ind4   | PopB | 150       | 150       | 200       | 206       | 302       | 302       |

**Notes:**

- First column: sample names.
- Second column: population name.
- Each locus requires **two adjacent columns** (one per allele).
- Empty cells indicate missing data.
- Header row starts at row 2 (first row contains column labels).
- Alleles can be numeric or string values.

#### 2. Dominant Markers (AFLP, RAPD, etc.)

| Sample | Pop  | Locus1 | Locus2 | Locus3 |
|--------|------|--------|--------|--------|
| Ind1   | PopA | 1      | 0      | 0      |
| Ind2   | PopA | 0      | 1      | 0      |
| Ind3   | PopB | 1      | 1      | 0      |
| Ind4   | PopB | 0      | 0      | 1      |

**Notes:**

- First column: sample names.
- Second column: population name.
- Each column represents one locus.
- Values: `1` (band present) or `0` (band absent).
- Empty cells allowed for missing data.

#### 3. Haploid Markers (cpDNA, mtDNA, Y-chromosome)

| Sample | Pop  | Locus1 | Locus2 | Locus3 |
|--------|------|--------|--------|--------|
| Ind1   | PopA | 100    | 90     | 200    |
| Ind2   | PopA | 101    | 91     | 202    |
| Ind3   | PopB | 100    | 91     | 198    |
| Ind4   | PopB | 103    | 90     | 198    |

**Notes:**

- First column: sample names.
- Second column: population name.
- Each column represents one locus (single allele per locus).
- Alleles can be numeric or character codes.

### Morphometric Data Files

| Sample  | Trait1 | Trait2 | Trait3 | Trait4 | ... |
|---------|--------|--------|--------|--------|-----|
| Sample1 | 15.2   | 8.4    | 12.1   | 5.3    | ... |
| Sample1 | 15.1   | 8.5    | 12.0   | 5.2    | ... |
| Sample1 | 15.3   | 8.3    | 12.2   | 5.4    | ... |
| Sample2 | 14.8   | 7.9    | 11.5   | 4.8    | ... |
| Sample2 | 14.9   | 8.0    | 11.6   | 4.9    | ... |
| Sample2 | 14.7   | 7.8    | 11.4   | 4.7    | ... |

**Notes:**

- First column **must** be named `Sample`.
- Each sample appears in **consecutive rows** (replicates).
- Specify the number of replicates per sample in the interface.
- Replicates can be averages, but individual measurements are preferred for proper statistics.
- Missing values allowed (will be imputed for dendrogram and PCA).

---

## 🧬 Genetic Analysis

### General Workflow

1. **Load Reference Database** — your reference population data.
2. **Load Query Samples** — unknown samples to analyze.
3. **Select Analysis Options** — choose bootstrap if desired.
4. **Run Analysis** — click the appropriate buttons in sequence.

### Analysis Sequence for Each Marker Type

#### 1. Basic Analysis

**Generate Similarity Matrix & Dendrogram**

- Calculates pairwise similarities and creates a UPGMA tree.
- Shows Top 3 most similar references for each query.
- Generates complete similarity matrix.
- Creates dendrogram with optional bootstrap support.

#### 2. Save Results

**Save Top 3 Results**

- Exports similarity results to a TSV file.
- Results automatically saved in `GMDA_[MarkerType]_Outputs/`.

#### 3. Complete Query Analysis

**Query Data Analysis** — comprehensive analysis including:

- Missing data statistics
- Genetic diversity indices (allele counts, heterozygosity, etc.)
- Similarity matrix for query samples only
- Dendrogram of query samples
- PCoA of query samples
- All results saved automatically

#### 4. Advanced Analyses

**DAPC** — Discriminant Analysis of Principal Components

- Requires population assignment.
- Visualizes population structure.
- Saves coordinates and plot.

**FST Outlier Detection**

- Identifies loci under selection.
- Requires population assignment.
- Provides 95% and 99% confidence thresholds.
- Saves histogram and scatter plots.

**Kinship Analysis**

- Estimates relatedness coefficients.
- Generates kinship matrix heatmap.
- Identifies closely related individuals.

### Genetic Diversity Calculations

| Marker type | Metrics |
|---|---|
| **Codominant** | A (alleles per locus), Ae (effective alleles, 1/Σp²), Ho (observed heterozygosity), He (expected heterozygosity), F (inbreeding coefficient) |
| **Dominant** | A (1 or 2 per locus), Ae (effective alleles), p (dominant allele frequency), q (recessive allele frequency), H (genetic diversity = 2pq) |
| **Haploid** | A (alleles per locus), Ae (effective alleles), h (gene diversity = 1 − Σp²) |

---

## 📊 Morphometric Analysis

### Workflow

1. **Load Data** — select Excel file with morphometric measurements.
2. **Set Replicates** — specify number of replicates per sample.
3. **Run Analysis** — click **Analyze Data**.

### Analysis Outputs

#### Descriptive Statistics

For each morphological trait:

- Mean ± standard deviation
- Median ± median absolute deviation
- Mode (or "multiple" if multimodal)
- Sample size (N)

#### ANOVA and Tukey HSD

- One-way ANOVA comparing all samples.
- Post-hoc Tukey HSD when significant (p < 0.05).
- Letter grouping for mean comparisons.
- Effect sizes (Cohen's d) between all pairs.

#### Multivariate Analyses

- **UPGMA Dendrogram** — hierarchical clustering based on Euclidean distance; optional bootstrap support; complete linkage method.
- **Principal Component Analysis (PCA)** — correlation biplot with variable loadings; explained variance per component; sample scores for visualization.

### Missing Data Handling

- Reports missing data by sample and trait.
- Mean imputation applied for multivariate analyses.
- Missing values excluded from pairwise comparisons in similarity calculations.

---

## 🔗 Correlation Analysis

### Available Mantel Tests

**Dendrogram Comparisons**

- Codominant genetic vs morphometric
- Dominant genetic vs morphometric
- Haploid genetic vs morphometric
- Inter-genetic comparisons (Codominant vs Dominant, etc.)

**Ordination Comparisons**

- Genetic PCoA vs morphometric PCA
- PCoA vs PCoA between genetic marker types

### Test Interpretation

**r (Mantel correlation):** ranges from −1 to 1

| Range | Interpretation |
|---|---|
| > 0.5 | Strong positive correlation |
| 0.3 – 0.5 | Moderate positive correlation |
| 0.1 – 0.3 | Weak positive correlation |
| < 0.1 | No correlation |

**p-value:**

- `p < 0.05` — statistically significant correlation.
- `p < 0.01` — highly significant correlation.

### Prerequisites

- Must run **both** analyses (genetic and morphometric) before correlation.
- Only **common samples** between datasets are used.
- Minimum **3 common samples** required for the test.

---

## 📂 Output Files

All outputs are saved in organized folders based on analysis type.

### Directory Structure

```text
[Input_File_Directory]/
├── GMDA_Genetic_codominant_Outputs/
│   ├── codominant_genetic_similarity_matrix.xlsx
│   ├── codominant_genetic_dendrogram.png
│   ├── codominant_genetic_dendrogram_with_bootstrap.png
│   ├── Top3_codominant.tsv
│   ├── codominant_genetic_indices.tsv
│   ├── codominant_genetic_query_dendrogram.png
│   ├── codominant_genetic_query_PCoA.png
│   ├── codominant_genetic_query_PCoA_coordinates.xlsx
│   ├── codominant_DAPC.png
│   ├── codominant_DAPC_coordinates.xlsx
│   ├── codominant_FST_outliers.png
│   ├── codominant_FST_outliers.xlsx
│   └── codominant_kinship_matrix.png
│
├── GMDA_Genetic_dominant_Outputs/
│   └── [similar structure for dominant markers]
│
├── GMDA_Genetic_haploid_Outputs/
│   └── [similar structure for haploid markers]
│
├── GMDA_Morphometric_Outputs/
│   ├── Morphometric_results.xlsx
│   ├── Morphometric_Dendrogram.png
│   ├── Morphometric_Dendrogram_with_bootstrap.png
│   └── Morphometric_PCA.png
│
└── GMDA_Mantel_correlations_Outputs/
    ├── Mantel_Test_*.pdf
    └── [various Mantel test outputs]
```

### File Types

| Extension | Description |
|---|---|
| `.xlsx` | Excel files with complete results |
| `.png` | Publication-quality figures (300 DPI) |
| `.tsv` | Tab-separated text files for easy import |
| `.pdf` | Mantel test results with scatter plots |

---

## 🔧 Troubleshooting

### Common Issues and Solutions

| # | Issue | Solution |
|---|---|---|
| 1 | **Error reading Excel files** | Ensure files are saved as `.xlsx` or `.xls`. Check that headers are correctly formatted. Verify no special characters in sample names. |
| 2 | **Missing data warnings** | Normal behavior — software handles missing data. Check warnings to ensure missing data is acceptable. Large amounts (>10%) may affect results. |
| 3 | **Dendrogram errors** | Minimum 2 samples required for dendrogram. Minimum 3 samples for PCoA and Mantel tests. Check for completely missing samples. |
| 4 | **ANOVA fails** | Need at least 2 groups with data. Groups with single individuals won't be analyzed. Check replicates per sample. |
| 5 | **Bootstrap errors** | Bootstrap requires sufficient data. May fail with very small datasets. Try disabling bootstrap for small samples. |

---

## ✅ Best Practices

### Data Preparation

- **Clean data** — remove obvious outliers before analysis.
- **Consistent naming** — use identical sample names across datasets.
- **Missing data** — keep missing values as empty cells (not `0` or `-99`).
- **Replicates** — for morphometrics, include all measurements, not just means.

### Analysis Order

1. Run genetic analyses first (codominant, dominant, or haploid).
2. Run morphometric analysis.
3. Run correlation analyses last (requires previous results).

### Interpretation Tips

- Bootstrap support **> 70%** is considered reliable.
- Mantel test **p < 0.05** indicates significant correlation.
- DAPC works best with **2–10 populations**.
- FST outliers above the **99% threshold** are strong evidence of selection.

---

## 📝 Citation

If you use GMDA in your research, please cite:

> STEFENON, V.M.; POLETO, T.; GIRARDELLO, G. M.; THALMAYR, P. **Selection of Pecan Genotypes using GMDA: a Software for Genetic and Morphometric Data Analyses.** *Crop Breeding and Applied Biotechnology*, v. 26, p. e54342625, 2026.

---

## 📧 Support

For questions, bug reports, or feature requests:

- Open an issue on GitHub.
- Contact: **valdir.stefenon@ufsc.br**

---

*GMDA — Genetic and Morphometric Data Analysis Platform*
