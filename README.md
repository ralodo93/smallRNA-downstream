# miRNA Differential Expression Analysis

A modular and reproducible R workflow for exploratory analysis and differential expression of miRNA count data.

The pipeline supports:

* Quality control and exploratory analysis
* Count filtering and TMM normalization
* Principal Component Analysis (PCA)
* Hierarchical clustering
* Expression distribution assessment
* Library size evaluation
* miRNA abundance visualization
* Differential expression analysis using the limma-voom framework
* Multiple experimental designs through configuration files
* Reproducible outputs and intermediate objects

---

# Overview

This repository provides a configurable workflow for the analysis of miRNA count matrices obtained from sequencing or public repositories such as GEO.

The analysis is divided into two independent stages:

1. Exploratory analysis
2. Differential expression analysis

Configuration is handled through a single `config.yml` file, allowing the same codebase to be reused across different datasets and experimental designs.

---

# Repository Structure

```text
.
├── config.yml
│
├── data/
│   └── example/
│       ├── example_counts.csv
│       └── example_metadata.csv
│
├── scripts/
│   ├── 00_download_geo_data.R
│   ├── exploratory_analysis.R
│   └── differential_expression.R
│
├── R/
│   ├── load_raw_data.R
│   ├── dge.R
│   ├── pca.R
│   ├── clustering.R
│   ├── distributions.R
│   ├── library_size.R
│   ├── abundance.R
│   ├── diff_expression.R
│   └── utils.R
│
├── results/
│
└── README.md
```

---

# Input Data

The workflow requires two input files.

## Counts Matrix

A CSV file containing raw miRNA counts.

Example:

| gene_id    | sample1 | sample2 | sample3 |
| ---------- | ------- | ------- | ------- |
| hsa-miR-21 | 125     | 98      | 143     |
| hsa-miR-16 | 455     | 390     | 502     |

Rows correspond to miRNAs and columns correspond to samples.

---

## Metadata

A CSV file describing samples and experimental variables.

Example:

| id   | trimester | fetal_sex |
| ---- | --------- | --------- |
| GSM1 | First     | Female    |
| GSM2 | Third     | Male      |

The `id` column must match the sample names present in the count matrix.

---

# Configuration

All analysis parameters are specified in `config.yml`.

Example:

```yaml
metadata_file: data/example/example_metadata.csv
counts_matrix_file: data/example/example_counts.csv

all_variables:
  - trimester
  - fetal_sex

interest_variable: trimester

interest_pcs:
  - PC1
  - PC2

distances:
  - corr
  - euclidean

clustering_methods:
  - complete
  - average

formulas:
  trimester: "~ 0 + trimester"
  sex_trimester: "~ 0 + trimester + fetal_sex"
```

This design allows exploratory and differential expression analyses to be performed without modifying the source code.

---

# Exploratory Analysis

The exploratory workflow performs the following analyses.

## Filtering

Lowly expressed miRNAs are removed using:

```r
filterByExpr()
```

---

## Normalization

Library size normalization is performed using the TMM method:

```r
calcNormFactors()
```

Log2 CPM values are subsequently computed for downstream analyses.

---

## Library Size Assessment

Raw and filtered library sizes are evaluated for each sample.

Outputs:

```text
results/exploratory_analysis/library_sizes/
```

---

## Expression Distributions

Per-sample expression distributions are visualized using normalized log2 CPM values.

Outputs:

```text
results/exploratory_analysis/distributions/
```

---

## Principal Component Analysis (PCA)

PCA is performed on normalized expression values.

Outputs:

```text
results/exploratory_analysis/pca/
```

---

## Hierarchical Clustering

Hierarchical clustering can be performed using:

* Pearson correlation distance
* Euclidean distance

Supported linkage methods:

* complete
* average

Outputs:

```text
results/exploratory_analysis/clustering/
```

---

## miRNA Abundance

Visualization of the most abundant miRNAs across samples and experimental groups.

Outputs:

```text
results/exploratory_analysis/abundance/
```

---

# Differential Expression Analysis

Differential expression analysis is performed using the limma-voom framework.

Workflow:

```r
voom()
lmFit()
contrasts.fit()
eBayes()
```

The design matrix is automatically generated from formulas defined in `config.yml`.

Multiple models can therefore be tested without modifying the analysis code.

Example:

```yaml
formulas:
  trimester: "~ 0 + trimester"
  sex_trimester: "~ 0 + trimester + fetal_sex"
```

---

# Differential Expression Outputs

For each contrast, the workflow generates:

## Complete Differential Expression Table

```text
DE_*.csv
```

Contains all tested miRNAs.

---

## Significant Differentially Expressed miRNAs

```text
DEG_*.csv
```

Contains significant miRNAs after multiple testing correction.

Significance criteria:

* Adjusted p-value < 0.05
* |logFC| > 1

---

## Volcano Plots

```text
DEG_*_volcano.png
```

Visualization of differential expression results.

---

## Heatmaps

```text
DEG_*_heatmap.png
```

Expression heatmaps of significant miRNAs.

---

# Running the Workflow

## Exploratory Analysis

```bash
Rscript scripts/exploratory_analysis.R
```

---

## Differential Expression Analysis

```bash
Rscript scripts/differential_expression.R
```

---

# GEO Example Dataset

An example GEO dataset can be downloaded using:

```bash
Rscript scripts/00_download_geo_data.R
```

The script generates:

```text
data/example/
├── example_counts.csv
└── example_metadata.csv
```

which can be used directly to reproduce the analyses contained in this repository.

---

# Intermediate Objects

The workflow stores reusable intermediate objects:

```text
results/rds/
├── dge.rds
└── log2_cpm.rds
```

These files allow downstream analyses to be performed without repeating filtering and normalization steps.

---

# Main Dependencies

CRAN:

* tidyverse
* ggplot2
* pheatmap
* ggdendro
* cowplot
* ggpubr
* yaml

Bioconductor:

* edgeR
* limma
* GEOquery

---

# Citation

If you use this repository, please cite:

Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics (2010).

Ritchie ME et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research (2015).

---

# License

This project is distributed under the terms of the MIT License.
