library(Seurat)
library(Matrix)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)


## ---- Load and Preprocess Gene and Barcode Files ----
genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Remove last 3 metadata rows (pCR, Cell_Type, patient_id)
genes <- genes[1:(nrow(genes) - 3), , drop = FALSE]
print(dim(genes))  # Expected: (25288, 1)

## ---- Load Raw Matrix and Assign Row/Column Names ----
raw_matrix <- readMM("bassez_raw_matrix.mtx")
rownames(raw_matrix) <- genes$V1
colnames(raw_matrix) <- barcodes$V1
print(dim(raw_matrix))


## ---- Batch the Matrix to Create Seurat Objects ----
batch_size <- 30000
num_batches <- ceiling(ncol(raw_matrix) / batch_size)

for (i in 1:num_batches) {
  cat("Processing batch", i, "of", num_batches, "\n")
  start_cell <- ((i - 1) * batch_size) + 1
  end_cell <- min(i * batch_size, ncol(raw_matrix))
  raw_matrix_subset <- raw_matrix[, start_cell:end_cell]
  seurat_batch <- CreateSeuratObject(counts = raw_matrix_subset)
  saveRDS(seurat_batch, file = paste0("07-29-30000_batch_", i, ".rds"))
}

# Merge all batches
batch_files <- list.files(pattern = "07-29-30000_batch_.*.rds")
batch_list <- lapply(batch_files, readRDS)
seurat_obj <- merge(batch_list[[1]], y = batch_list[-1], add.cell.ids = 1:length(batch_list))
saveRDS(seurat_obj, file = "07-29-30000_merged.rds")


## ---- Normalize, Scale, and Dimensional Reduction ----
seurat_obj <- readRDS("07-29-30000_merged.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)
ElbowPlot(seurat_obj)


## ---- Clustering and Visualization ----
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.125)

cluster_assignments <- data.frame(Cell = colnames(seurat_obj),
                                  Cluster = Idents(seurat_obj))
write.csv(cluster_assignments, "07-29-12PC-30000_14clusters.csv", row.names = FALSE)

# Run and visualize UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:12)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, raster = FALSE) +
  ggtitle("UMAP of Clustering Results (Resolution = 0.125)") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))


## ---- Identify Marker Genes ----
seurat_obj <- JoinLayers(seurat_obj, layers = c("counts", "data"))
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
write.csv(markers, "07-29-12PC-30000_14clusters_markers.csv", row.names = FALSE)


## ---- Annotate Clusters with SingleR ----
# Build a binary presence matrix of genes across clusters
pseudo_bulk <- all_markers %>%
  select(cluster, gene) %>%
  distinct() %>%
  mutate(pseudo_expr = 1) %>%
  pivot_wider(
    names_from = cluster,
    values_from = pseudo_expr,
    values_fill = list(pseudo_expr = 0)
  )

# Convert to numeric matrix with genes as rownames
pseudo_bulk_mat <- as.data.frame(pseudo_bulk)
rownames(pseudo_bulk_mat) <- pseudo_bulk_mat$gene
pseudo_bulk_mat$gene <- NULL
pseudo_bulk_mat[] <- lapply(pseudo_bulk_mat, as.numeric)
pseudo_bulk_mat <- as.matrix(pseudo_bulk_mat)

# Load reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR annotation
singleR_result <- SingleR(test = pseudo_bulk_mat, ref = ref, labels = ref$label.main)
table(singleR_result$labels)

# Save cluster-to-cell type predictions
cluster_annotations <- data.frame(
  cluster = colnames(pseudo_bulk_mat),
  predicted_cell_type = singleR_result$labels
)
write.csv(cluster_annotations, "07-29-12PC-14clusters_pred_celltypes.csv", row.names = FALSE)


## ---- Map Predicted Cell Types to Seurat Object ----
predicted <- read.csv("07-29-12PC-14clusters_pred_celltypes.csv")
predicted$cluster <- as.character(predicted$cluster)

# Restore cluster identities if overwritten
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data) ||
    length(unique(seurat_obj$seurat_clusters)) == 1) {
  cluster_info <- read.csv("07-29-12PC-30000_14clusters.csv")
  rownames(cluster_info) <- cluster_info$Cell
  seurat_obj$seurat_clusters <- cluster_info[colnames(seurat_obj), "Cluster"]
}

seurat_obj$seurat_clusters <- as.character(seurat_obj$seurat_clusters)
matched_indices <- match(seurat_obj$seurat_clusters, predicted$cluster)
seurat_obj$predicted_cell_type <- predicted$predicted_cell_type[matched_indices]
table(seurat_obj$predicted_cell_type, useNA = "ifany")


## ---- UMAP Visualization by Predicted Cell Types ----
seurat_obj$predicted_cell_type <- dplyr::recode(seurat_obj$predicted_cell_type,
                                                "B_cell" = "B cells",
                                                "Embryonic_stem_cells" = "Embryonic stem cells",
                                                "Endothelial_cells" = "Endothelial cells",
                                                "Epithelial_cells" = "Epithelial cells",
                                                "Erythroblast" = "Erythroblasts",
                                                "Fibroblasts" = "Fibroblasts",
                                                "Macrophage" = "Macrophages",
                                                "Monocyte" = "Monocytes",
                                                "NK_cell" = "NK cells",
                                                "T_cells" = "T cells",
                                                "NA" = "Unassigned",
                                                "iPS_cells" = "iPS cells")

custom_colors <- c(
  "T cells"              = "#018bb2",  
  "NK cells"             = "#806e49",  
  "Epithelial cells"     = "#cc77ab", 
  "Fibroblasts"          = "#936bd6",  
  "Keratinocytes"        = "#9dbfd0",  
  "Macrophages"          = "#5cbd4d",  
  "Monocytes"            = "#265734",  
  "Embryonic stem cells" = "#ee9a34", 
  "iPS cells"            = "#ea7f5a",  
  "Endothelial cells"    = "#f2cf63", 
  "B cells"              = "#6c0297",
  "Neurons"              = "#591834",
  "Erythroblasts"        = "#e63973",  
  "Osteoblasts"          = "#1f4e8c",  
  "Hepatocytes"          = "#3eb489",  
  "Astrocytes"           = "#8e8e8e"    
)

DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "predicted_cell_type",
  cols = custom_colors,
  raster = FALSE
) + ggtitle("UMAP of Predicted Cell Types from 14 Clusters") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )


## ---- Merge Marker Genes with Predicted Cell Types ----
all_markers$cluster <- as.character(all_markers$cluster)
predicted$cluster <- as.character(predicted$cluster)
markers_with_celltypes <- left_join(all_markers, predicted, by = "cluster")
write.csv(markers_with_celltypes, "07-29-12PC-14clusters_markers_with_celltypes.csv", row.names = FALSE)


## ---- Filter and Rank Marker Genes ----
markers_cta <- read.csv("07-29-12PC-14clusters_markers_with_celltypes.csv")

# Filter marker genes based on log2FC, adjusted p-value, and expression threshold
filtered_markers <- markers_cta %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05, pct.1 >= 0.25)

# Identify cell types with no genes passing the filter
missing_cell_types <- setdiff(
  unique(markers_cta$predicted_cell_type),
  unique(filtered_markers$predicted_cell_type)
)

# Select top 100 genes for each missing cell type based on ranking criteria
top_genes <- markers_cta %>%
  filter(predicted_cell_type %in% missing_cell_types) %>%
  arrange(predicted_cell_type, p_val_adj, desc(avg_log2FC), desc(pct.1)) %>%
  group_by(predicted_cell_type) %>%
  slice_head(n = 100) %>%
  ungroup()

# Combine filtered and top-ranked genes
final_markers <- bind_rows(filtered_markers, top_genes) %>%
  select(-X)
write.csv(final_markers, "bassez_raw_07-29-12PC-14clusters_markers_with_CTA_filtered_ranked.csv", row.names = FALSE)
