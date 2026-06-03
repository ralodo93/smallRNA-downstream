load_raw_data <- function(metadata_file,
                          counts_matrix_file,
                          all_variables) {
  metadata <- read_csv(metadata_file) %>% as.data.frame()
  
  # check if all variables are present in the metadata
  if (!all(all_variables %in% colnames(metadata))) {
    stop("Not all variables are present in the metadata file.")
  }
  
  counts_matrix <- read_csv(counts_matrix_file) %>% as.data.frame()
  
  # order the columns of the counts matrix according to the order of the samples in the metadata
  counts_matrix <- counts_matrix[, c("gene_id", metadata$id)]
  
  # convert the counts matrix to a numeric matrix
  rownames(counts_matrix) <- counts_matrix$gene_id
  counts_matrix <- as.matrix(counts_matrix[, -1])
  
  # set the row names of the metadata to the sample ids
  rownames(metadata) <- metadata$id
  
  # return the metadata and counts matrix as a list
  return(list(metadata = metadata, counts_matrix = counts_matrix))
}