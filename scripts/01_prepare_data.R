library(yaml)
library(tidyverse)
library(edgeR)

config <- read_yaml("config.yml")

metadata <- read_delim(config$metadata_file)
counts <- read_delim(config$counts_matrix_file)

gene_col <- colnames(counts)[1]
counts <- counts |> column_to_rownames(var = gene_col)

metadata <- metadata |> column_to_rownames("id") |>
  mutate(across(everything(), as.factor))

if (!all(colnames(counts) %in% rownames(metadata))) {
  stop("colnames from counts and metadata id columns must be the same")
}

# create dge object
d0 <- DGEList(counts = counts, samples = metadata)

# calculate normalization factors to scale the raw library sizes (TMM)
d0 <- calcNormFactors(d0)

# rowSums(cpms <=1) < 3 , require at least 1 cpm in at least 3 samples to keep.
cutoff <- 3
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]

# Normalization
log2_cpm <- cpm(d, prior.count = 2, log = TRUE)

dir.create("results/rds", showWarnings = FALSE, recursive = TRUE)
save(metadata, counts, d, log2_cpm, file = "results/rds/raw.rdata")
