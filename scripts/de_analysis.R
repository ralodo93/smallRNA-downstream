library(edgeR)
library(limma)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(yaml)
library(aplot)

config <- read_yaml("config.yml")

formulas           <- config$formulas
color_palette      <- config$color_palette


source("R/utils.R")
dir.create("results/differential_expression", showWarnings = FALSE)

dge <- readRDS("results/rds/dge.rds")
design_metadata <- dge$samples

source("R/diff_expression.R")
for (formula_name in names(formulas)){
  formula <- formulas[[formula_name]]
  diff_expression(formula, formula_name, design_metadata, dge, color_palette)
}
