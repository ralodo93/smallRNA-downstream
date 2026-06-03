# load required libraries
library(GEOquery)
library(dplyr)
library(readr)
library(tibble)
library(rtracklayer)
library(GenomicRanges)

# id of the study
gse_id <- "GSE184860"

# download phenoData
gse_list <- getGEO(gse_id)
geo_object <- gse_list[[1]]
pheno <- pData(geo_object)
raw_files <- getGEOSuppFiles(gse_id)


# read raw files
file_path <- rownames(raw_files)
raw_counts <- read_csv(file_path)
colnames(raw_counts)[1] <- "gene_id"
colnames(raw_counts) <- gsub(".UMIs", "", colnames(raw_counts))

raw_counts <- raw_counts[,c("gene_id", pheno$description.1)]
colnames(raw_counts) <- c("gene_id", pheno$geo_accession)

pheno <- pheno[,c("geo_accession","fetal sex:ch1", "trimester:ch1")]
colnames(pheno) <- c("id", "fetal_sex", "trimester")

raw_counts <- raw_counts[grep("hsa-", raw_counts$gene_id),]

dir.create("data/example", showWarnings = FALSE)
write_csv(raw_counts, file = "data/example/example_counts.csv")
write_csv(pheno, file = "data/example/example_metadata.csv")
