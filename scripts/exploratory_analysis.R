library(tidyverse)
library(edgeR)
library(pheatmap)
library(ggdendro)
library(yaml)

config <- read_yaml("config.yml")

metadata_file      <- config$metadata_file
counts_matrix_file <- config$counts_matrix_file
all_variables      <- config$all_variables
interest_variable  <- config$interest_variable
interest_pcs       <- config$interest_pcs
distances          <- config$distances
clustering_methods <- config$clustering_methods
color_palette      <- config$color_palette

source("R/load_raw_data.R")
source("R/dge.R")
source("R/utils.R")


# 1. Load Data -----------------------------------------------------------------

data <- load_raw_data(metadata_file, counts_matrix_file, all_variables)
metadata <- data$metadata
counts_matrix <- data$counts_matrix

# 2. Prepare DGE Object
dge_obj <- create_dge(counts_matrix, metadata, interest_variable)
dge <- dge_obj$dge
log2_cpm <- dge_obj$log2_cpm

## 2.1 Save intermediate objects

dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)
saveRDS(dge, file = "results/rds/dge.rds")
saveRDS(log2_cpm, file = "results/rds/log2_cpm.rds")


# 5. Exploratory Analysis
dir.create("results/exploratory_analysis/", showWarnings = FALSE)

## 5.1 Library Size QC
dir.create("results/exploratory_analysis/library_sizes", showWarnings = FALSE)
source("R/library_size.R")
library_sizes(counts_matrix, dge, metadata, all_variables, color_palette)


## 5.2 Distribution Plots
dir.create("results/exploratory_analysis/distributions", showWarnings = FALSE)
source("R/distributions.R")
get_distribution(log2_cpm, metadata, all_variables, color_palette)


## 5.3 PCA
dir.create("results/exploratory_analysis/pca", showWarnings = FALSE)
source("R/pca.R")
process_pca(log2_cpm, metadata, all_variables, interest_pcs, color_palette)


## 5.4 Clustering
dir.create("results/exploratory_analysis/clustering", showWarnings = FALSE)
source("R/clustering.R")
process_clustering(log2_cpm, metadata, all_variables, distances, clustering_methods, color_palette)


## 5.5 Abundance
dir.create("results/exploratory_analysis/abundance", showWarnings = FALSE)
source("R/abundance.R")
process_abundance(counts_matrix, log2_cpm, metadata, all_variables, color_palette, n = 20)
