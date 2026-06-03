create_dge <- function(counts, metadata, interest_variable) {
  dge <- DGEList(counts = counts_matrix, samples = metadata)
  
  # Filter out low-expressed genes using edgeR's default filtering criteria
  # Note: metadata[, interest_variable] extracts the grouping factor
  keep <- filterByExpr(dge, group = metadata[, interest_variable])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization
  dge <- calcNormFactors(dge, method = "TMM")
  log2_cpm <- cpm(dge, log = TRUE, prior.count = 2)
  
  # Return the DGEList object and the log2 CPM matrix
  return(list(dge = dge, log2_cpm = log2_cpm))
}