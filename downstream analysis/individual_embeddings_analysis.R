library(arrow)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(tidyr)
library(readr)


## ---- Load and Filter Embeddings (Mixedbread) ----
# Read embedding matrix and remove excluded patients
embeddings <- read_feather("bassez_mxbai_emb_768.feather")
removed_idx <- scan("removed_patient_indices.txt", quiet = TRUE)
embeddings <- embeddings[-removed_idx, ]


## ---- Load Alternative Embeddings (scGPT) ----
# Uncomment if using scGPT embeddings instead of Mixedbread
# embeddings <- read_feather("bassez_scgpt_emb_768.feather")


## ---- Initialize Seurat Object from Embeddings ----
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Convert to matrix and align with barcodes
embeddings_mat <- as.matrix(embeddings)
rownames(embeddings_mat) <- barcodes$V1
embeddings_sparse <- as(embeddings_mat, "dgCMatrix")

# Create Seurat object with embeddings as counts (transpose to match format)
seurat_obj <- CreateSeuratObject(counts = t(embeddings_sparse), project = "BreastCancer", assay = "RNA")
rownames(embeddings_matrix) <- colnames(seurat_obj)
colnames(embeddings_matrix) <- 1:ncol(embeddings_matrix)
seurat_obj[["embeddings"]] <- CreateDimReducObject(embeddings = embeddings_matrix, key = "embed_", assay = "RNA")


## ---- Perform Clustering using Louvain Algorithm ----
options(future.globals.maxSize = 10 * 1024^3)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "embeddings", dims = 1:768)

# Try a range of resolutions to determine the optimal resolution
resolutions <- seq(0.005, 0.1, by = 0.005)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  print(paste("Resolution:", res, "Number of clusters:", length(unique(seurat_obj$seurat_clusters))))
}

# Perform clustering at the optimal resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.230)


## ---- Run UMAP and Visualize Clusters----
seurat_obj <- RunUMAP(seurat_obj, reduction = "embeddings", dims = 1:768)

DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  raster = FALSE
) + ggtitle("UMAP of Clustering Results (Resolution = 0.230)") + 
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size=12) 
  )

table(seurat_obj$seurat_clusters)


## ---- Load Raw Matrix and Assign Clusters ----
# Load count matrix and corresponding gene/barcode info
raw_matrix <- readMM("cell_gene_pcr_matrix_extended.mtx")
genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_matrix) <- genes$V1
colnames(raw_matrix) <- barcodes$V1


# Cluster assignment metadata
cluster_assignments <- data.frame(Cell = colnames(seurat_obj), Cluster = Idents(seurat_obj))
#write.csv(cluster_assignments, "bassez_scgpt_06-27_9clusters.csv", row.names = FALSE)
write.csv(cluster_assignments, "bassez_mxbai_06-28_14clusters.csv", row.names = FALSE)

# Create Seurat object from raw counts and add cluster info
seurat_raw <- CreateSeuratObject(counts = raw_matrix)
rownames(cluster_assignments) <- cluster_assignments$Cell
seurat_raw$seurat_clusters <- cluster_assignments[colnames(seurat_raw), "Cluster"]


## ---- Normalize Data and Identify Marker Genes ----
seurat_raw <- NormalizeData(seurat_raw)
seurat_raw <- FindVariableFeatures(seurat_raw)

markers <- FindAllMarkers(seurat_raw, only.pos = TRUE, group.by = "seurat_clusters")
#write.csv(markers, "bassez_scgpt_06-27_9clusters_all_markers.csv", row.names = FALSE)
write.csv(markers, "bassez_mxbai_06-28_14clusters_all_markers.csv", row.names = FALSE)


## ---- Create Pseudo-Bulk Matrix and Annotate with SingleR ----
# Load previously saved marker genes
#all_markers <- read.csv("bassez_scgpt_06-27_9clusters_all_markers.csv")
all_markers <- read.csv("bassez_mxbai_06-28_14clusters_all_markers.csv")

# Load reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Create pseudo-bulk matrix
pseudo_bulk <- all_markers %>%
  select(cluster, gene) %>%
  distinct() %>%
  mutate(pseudo_expr = 1) %>%
  pivot_wider(names_from = cluster, values_from = pseudo_expr, values_fill = list(pseudo_expr = 0))

pseudo_bulk_mat <- as.matrix(as.data.frame(pseudo_bulk[-1]))
rownames(pseudo_bulk_mat) <- pseudo_bulk$gene

# Annotate clusters with SingleR
singleR_result <- SingleR(test = pseudo_bulk_mat, ref = ref, labels = ref$label.main)
table(singleR_result$labels)
cluster_annotations <- data.frame(cluster = colnames(pseudo_bulk_mat), predicted_cell_type = singleR_result$labels)
#write.csv(cluster_annotations, "bassez_scgpt_06-27_9clusters_all_markers_pred_celltypes.csv", row.names = FALSE)
write.csv(cluster_annotations, "bassez_mxbai_06-28_14clusters_all_markers_pred_celltypes.csv", row.names = FALSE)


## ---- Summarize Predicted Cell Types ----
#predicted <- read.csv("bassez_scgpt_06-27_9clusters_all_markers_pred_celltypes.csv")
predicted <- read.csv("bassez_mxbai_06-28_14clusters_all_markers_pred_celltypes.csv")
cluster_sizes <- table(seurat_raw$seurat_clusters)
predicted$cluster <- as.character(predicted$cluster)
predicted$cluster_size <- cluster_sizes[predicted$cluster]

predicted_summary <- predicted %>%
  group_by(predicted_cell_type) %>%
  summarise(predicted_cells = sum(cluster_size)) %>%
  arrange(desc(predicted_cells))
print(predicted_summary)



## ---- Add Annotations and Visualize by Cell Type ----
seurat_raw$predicted_cell_type <- predicted$predicted_cell_type[match(seurat_raw$seurat_clusters, predicted$cluster)]
seurat_raw[["embeddings"]] <- seurat_obj[["embeddings"]]
seurat_raw <- RunUMAP(seurat_raw, reduction = "embeddings", dims = 1:768)


## ---- UMAP Visualization by Predicted Cell Types ----
seurat_raw$predicted_cell_type <- dplyr::recode(seurat_raw$predicted_cell_type,
                                                "B_cell" = "B cells",
                                                "Embryonic_stem_cells" = "Embryonic stem cells",
                                                "Endothelial_cells" = "Endothelial cells",
                                                "Epithelial_cells" = "Epithelial cells",
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
  seurat_raw,
  reduction = "umap",
  group.by = "predicted_cell_type",
  #label = TRUE,
  cols = custom_colors,
  raster = FALSE
) + ggtitle("UMAP of Predicted Cell Types from 9 Clusters") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))


## ---- Cluster Frequency Table ----
table(seurat_raw$predicted_cell_type, useNA = "ifany")
