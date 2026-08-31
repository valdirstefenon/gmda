# gmda
# GMDA - Genetic and Morphometric Data Analysis Platform

**Version 4.4**

GMDA is a comprehensive desktop application for analyzing genetic and morphometric data, designed for researchers in population genetics, evolutionary biology, and related fields. The software provides a user-friendly interface for performing complex statistical analyses, including similarity calculations, population structure analysis, and correlation studies.

## 📋 Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input File Formats](#input-file-formats)
- [Genetic Analysis](#genetic-analysis)
  - [Codominant Markers](#codominant-markers)
  - [Dominant Markers (AFLP)](#dominant-markers-aflp)
  - [Haploid Markers](#haploid-markers-cpdnay-chromosomemt-dna)
- [Morphometric Analysis](#morphometric-analysis)
- [Correlation Analysis](#correlation-analysis)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## 🚀 Features

### Genetic Analysis
- **Similarity Calculations**: Computes genetic similarity between query samples and reference databases
- **UPGMA Dendrograms**: Generates hierarchical clustering trees with optional bootstrap support
- **DAPC Analysis**: Discriminant Analysis of Principal Components for population structure
- **FST Outlier Detection**: Identifies loci under selection using FST distribution
- **Kinship Analysis**: Estimates relatedness between individuals
- **PCoA**: Principal Coordinate Analysis for visualizing genetic distances

### Morphometric Analysis
- **Descriptive Statistics**: Calculates mean, median, mode with standard deviations
- **ANOVA**: One-way ANOVA with Tukey HSD post-hoc tests
- **Effect Sizes**: Cohen's d with 95% confidence intervals
- **PCA**: Principal Component Analysis with variable loadings
- **UPGMA Dendrograms**: Hierarchical clustering of morphological traits
- **Bootstrap Support**: Optional bootstrap validation for dendrograms

### Correlation Analysis
- **Mantel Tests**: Correlation between genetic and morphometric distance matrices
- **Dendrogram Comparisons**: Mantel tests between different dendrogram types
- **PCoA/PCA Comparisons**: Correlation between ordination spaces
- **Statistical Significance**: Permutation-based p-values

## 💻 Installation

### Prerequisites
- Python 3.8 or higher
- Required Python packages:
numpy>=1.21.0
pandas>=1.3.0
openpyxl>=3.0.9
scipy>=1.7.0
scikit-learn>=0.24.0
matplotlib>=3.4.0
seaborn>=0.11.0
statsmodels>=0.12.0
packaging>=21.0

# Install required packages
pip install numpy>=1.21.0 pandas>=1.3.0 openpyxl>=3.0.9 scipy>=1.7.0 scikit-learn>=0.24.0 matplotlib>=3.4.0 seaborn>=0.11.0 statsmodels>=0.12.0 packaging>=21.0

# Run the application
python3 GMDA4.2.py


📁 Input File Formats
Genetic Data Files
All genetic files should be in Excel format (.xlsx or .xls) with the following structures:

1. Codominant Markers (e.g., microsatellites, SNPs)
Ex.
| Sample    | Pop       | Locus1_A1 | Locus1_A2 | Locus2_A1 | Locus2_A2 | Locus3_A1 | Locus3_A2 |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| Ind1      | PopA      | 150       | 150       | 200       | 204       | 300       | 302       |
| Ind2      | PopA      | 150       | 152       | 200       | 200       | 300       | 300       |
| Ind3      | PopB      | 148       | 150       | 202       | 204       | 300       | 304       |
| Ind4      | PopB      | 150       | 150       | 200       | 206       | 302       | 302       |

Important notes:

First column: Sample names

Second column: Population name

Each locus requires two adjacent columns (one per allele)

Empty cells for missing data (will be handled appropriately)

Header row starts at row 2 (first row contains column labels)

Alleles can be numeric values or strings



3. Dominant Markers (AFLP, RAPD, etc.)
Ex.
| Sample    | Pop       | Locus1  | Locus2  | Locus3  |
|-----------|-----------|---------|-------- |---------|
| Ind1      | PopA      | 1       | 0       | 0       |
| Ind2      | PopA      | 0       | 1       | 0       |
| Ind3      | PopB      | 1       | 1       | 0       |
| Ind4      | PopB      | 0       | 0       | 1       |

Important notes:

First column: Sample names

Second column: Population name

Each column represents one locus

Values: 1 (band present) or 0 (band absent)

Empty cells allowed for missing data

Haploid Markers (cpDNA, mtDNA, Y-chromosome)
Ex.
| Sample    | Pop       | Locus1  | Locus2  | Locus3  |
|-----------|-----------|---------|-------- |---------|
| Ind1      | PopA      | 100     | 90      | 200     |
| Ind2      | PopA      | 101     | 91      | 202     |
| Ind3      | PopB      | 100     | 91      | 198     |
| Ind4      | PopB      | 103     | 90      | 198     |

Important notes:

First column: Sample names

Second column: Population name

Each column represents one locus (single allele per locus)

Alleles can be numeric or character codes

Morphometric Data Files
text
| Sample  | Trait1 | Trait2 | Trait3 | Trait4 | ... |
|---------|--------|--------|--------|--------|-----|
| Sample1 | 15.2   | 8.4    | 12.1   | 5.3    | ... |
| Sample1 | 15.1   | 8.5    | 12.0   | 5.2    | ... |
| Sample1 | 15.3   | 8.3    | 12.2   | 5.4    | ... |
| Sample2 | 14.8   | 7.9    | 11.5   | 4.8    | ... |
| Sample2 | 14.9   | 8.0    | 11.6   | 4.9    | ... |
| Sample2 | 14.7   | 7.8    | 11.4   | 4.7    | ... |

Important notes:

First column MUST be named "Sample"

Each sample appears in consecutive rows (replicates)

Specify the number of replicates per sample in the interface

Replicates can be averages, but individual measurements preferred for proper statistics

Missing values allowed (will be imputed for dendrogram and PCA)

🧬 Genetic Analysis
General Workflow
Load Reference Database: Your reference population data

Load Query Samples: Unknown samples to analyze

Select Analysis Options: Choose bootstrap if desired

Run Analysis: Click appropriate buttons in sequence

Analysis Sequence for Each Marker Type
1. Basic Analysis
Generate Similarity Matrix & Dendrogram: Calculates pairwise similarities and creates UPGMA tree

Shows Top 3 most similar references for each query

Generates complete similarity matrix

Creates dendrogram with optional bootstrap support

2. Save Results
Save Top 3 Results: Exports similarity results to TSV file

Results automatically saved in GMDA_[MarkerType]_Outputs/ folder

3. Complete Query Analysis
Query Data Analysis: Comprehensive analysis including:

Missing data statistics

Genetic diversity indices (allele counts, heterozygosity, etc.)

Similarity matrix for query samples only

Dendrogram of query samples

PCoA of query samples

All results saved automatically

4. Advanced Analyses
DAPC: Discriminant Analysis of Principal Components

Requires population assignment

Visualizes population structure

Saves coordinates and plot

FST Outlier Detection:

Identifies loci under selection

Requires population assignment

Provides 95% and 99% confidence thresholds

Saves histogram and scatter plots

Kinship Analysis:

Estimates relatedness coefficients

Generates kinship matrix heatmap

Identifies closely related individuals

Genetic Diversity Calculations
Codominant Markers
A: Number of alleles per locus

Ae: Effective number of alleles (1/Σp²)

Ho: Observed heterozygosity

He: Expected heterozygosity (gene diversity)

F: Inbreeding coefficient

Dominant Markers
A: Number of alleles (1 or 2 per locus)

Ae: Effective number of alleles

p: Dominant allele frequency

q: Recessive allele frequency

H: Genetic diversity (2pq)

Haploid Markers
A: Number of alleles per locus

Ae: Effective number of alleles

h: Gene diversity (1 - Σp²)

📊 Morphometric Analysis
Workflow
Load Data: Select Excel file with morphometric measurements

Set Replicates: Specify number of replicates per sample

Run Analysis: Click "Analyze Data"

Analysis Outputs
Descriptive Statistics
For each morphological trait:

Mean ± Standard Deviation

Median ± Median Absolute Deviation

Mode (or "multiple" if multimodal)

Sample size (N)

ANOVA and Tukey HSD
One-way ANOVA comparing all samples

Post-hoc Tukey HSD when significant (p < 0.05)

Letter grouping for mean comparisons

Effect sizes (Cohen's d) between all pairs

Multivariate Analyses
UPGMA Dendrogram: Hierarchical clustering based on Euclidean distance

Optional bootstrap support

Complete linkage method

Principal Component Analysis (PCA):

Correlation biplot with variable loadings

Explained variance per component

Sample scores for visualization

Missing Data Handling
Reports missing data by sample and trait

Mean imputation applied for multivariate analyses

Missing values excluded from pairwise comparisons in similarity calculations

🔗 Correlation Analysis
Available Mantel Tests
Dendrogram Comparisons
Codominant genetic vs Morphometric

Dominant genetic vs Morphometric

Haploid genetic vs Morphometric

Inter-genetic comparisons (Codominant vs Dominant, etc.)

Ordination Comparisons
Genetic PCoA vs Morphometric PCA

PCoA vs PCoA between genetic marker types

Test Interpretation
r (Mantel correlation): Range from -1 to 1

0.5: Strong positive correlation

0.3-0.5: Moderate positive correlation

0.1-0.3: Weak positive correlation

<0.1: No correlation

p-value:

p < 0.05: Statistically significant correlation

p < 0.01: Highly significant correlation

Prerequisites
Must run both analyses (genetic and morphometric) before correlation

Only common samples between datasets are used

Minimum 3 common samples required for test

📂 Output Files
All outputs are saved in organized folders based on analysis type:

Directory Structure
text
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
│   ├── [similar structure for dominant markers]
│
├── GMDA_Genetic_haploid_Outputs/
│   ├── [similar structure for haploid markers]
│
├── GMDA_Morphometric_Outputs/
│   ├── Morphometric_results.xlsx
│   ├── Morphometric_Dendrogram.png
│   ├── Morphometric_Dendrogram_with_bootstrap.png
│   ├── Morphometric_PCA.png
│
└── GMDA_Mantel_correlations_Outputs/
    ├── Mantel_Test_*.pdf
    └── [various Mantel test outputs]
File Types
.xlsx: Excel files with complete results

.png: Publication-quality figures (300 DPI)

.tsv: Tab-separated text files for easy import

.pdf: Mantel test results with scatter plots

🔧 Troubleshooting
Common Issues and Solutions
1. Error reading Excel files
Ensure files are saved in .xlsx or .xls format

Check that headers are correctly formatted

Verify no special characters in sample names

2. Missing data warnings
Normal behavior - software handles missing data

Check warnings to ensure missing data is acceptable

Large amounts (>10%) may affect results

3. Dendrogram errors
Minimum 2 samples required for dendrogram

Minimum 3 samples for PCoA and Mantel tests

Check for completely missing samples

4. ANOVA fails
Need at least 2 groups with data

Groups with single individuals won't be analyzed

Check replicates per sample

5. Bootstrap errors
Bootstrap requires sufficient data

May fail with very small datasets

Try disabling bootstrap for small samples

Best Practices
Data Preparation
Clean data: Remove obvious outliers before analysis

Consistent naming: Use identical sample names across datasets

Missing data: Keep missing values as empty cells (not 0 or -99)

Replicates: For morphometrics, include all measurements, not just means

Analysis Order
Run genetic analyses first (codominant, dominant, or haploid)

Run morphometric analysis

Run correlation analyses last (requires previous results)

Interpretation Tips
Bootstrap support >70% considered reliable

Mantel test p < 0.05 indicates significant correlation

DAPC best with 2-10 populations

FST outliers >99% threshold strong evidence of selection

📝 Citation
If you use GMDA in your research, please cite:

STEFENON, V.M.; POLETO, T. ; GIRARDELLO, G. M. ; THALMAYR, P. . Selection of Pecan Genotypes using GMDA: a Software for Genetic and Morphometric Data Analyses. Crop Breeding and Applied Biotechnology, v. 26, p. e54342625, 2026.


📧 Support
For questions, bug reports, or feature requests:

Open an issue on GitHub

Contact: valdir.stefenon@ufsc.br



