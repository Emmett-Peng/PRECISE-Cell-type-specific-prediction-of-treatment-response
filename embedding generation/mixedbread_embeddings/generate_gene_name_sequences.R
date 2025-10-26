# Input:  Sparse expression matrix (.rds) where rows = genes, columns = cells
# Output: Text file where each line contains ranked gene names for a cell

library(Matrix)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_path <- ifelse(length(args) >= 1, args[1], "data/gene_expression/sparse_expression_matrix.rds")
output_path <- ifelse(length(args) >= 2, args[2], "data/gene_expression/gene_name_sequences.txt")

# Read expression matrix
expression_matrix <- readRDS(input_path)

gene_names <- rownames(expression_matrix)
cell_names <- colnames(expression_matrix)

# Generate ranked sequences
ranked_sequences <- lapply(seq_len(ncol(expression_matrix)), function(cell_index) {
  cell_data <- expression_matrix[, cell_index]
  nonzero_indices <- which(cell_data > 0)
  if (length(nonzero_indices) > 0) {
    ranked_indices <- nonzero_indices[order(cell_data[nonzero_indices], decreasing = TRUE)]
    paste(gene_names[ranked_indices], collapse = " ")
  } else {
    ""
  }
})

writeLines(ranked_sequences, output_path)

