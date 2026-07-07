library(tidyverse)

if (!file.exists("data/example/counts.txt")) {
  download.file(
    "https://raw.githubusercontent.com/sivkri/miRNASeq-DE-Analysis-QC/refs/heads/main/abi_v3.txt",
    "data/example/counts.txt"
  )
  counts <- read.delim("data/example/counts.txt")
  colnames(counts) <- c(
    "gene",
    "abi1_aba_1",
    "abi1_aba_2",
    "abi1_aba_3",
    "abi1_control_1",
    "abi1_control_2",
    "abi1_control_3",
    "wt_aba_1",
    "wt_aba_2",
    "wt_aba_3",
    "wt_control_1",
    "wt_control_2",
    "wt_control_3"
  )
  
  write.table(
    counts,
    file = "data/example/counts.txt",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )
  
  design <- data.frame(id = colnames(counts)[-1]) |>
    mutate(condition = ifelse(grepl("aba", id), "ABA", "Control"),
           genotype = ifelse(grepl("abi", id), "abi1", "wt"))
  
  write.table(
    design, file = "data/example/metadata.txt",
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
  )
  
}
