This repository contains all code, scripts, and results for the research project titled **PRECISE: Cell-Type Inference from Single-Cell Embeddings Enhances Prediction of Breast Cancer Immunotherapy Response.**  

---

## Repository Structure


### `embedding generation/`  
Contains scripts used to generate embeddings.

- **`scgpt_embeddings.ipynb`** - Colab Python notebook used to generate scGPT embeddings and save them for downstream processing in R. Embedding concatenation of scGPT and Mixedbread to generate scMIXE was performed in this notebook.
- **`mixedbread_embeddings.ipynb`** - 

---


### `downstream analysis/`  
Contains scripts for clustering, marker gene identification, and cell-type annotation, followed by shared downstream steps for marker filtering and gene set enrichment analysis. All embedding types (raw data, individual embeddings, and concatenated scMIXE embeddings) undergo the same filtering and enrichment pipeline.

- **`raw_data_analysis.R`** – Using raw data ([Bassez et al. data](https://lambrechtslab.sites.vib.be/en/single-cell))
- **`individual_embedding_analysis.R`** – Using scGPT or Mixedbread embeddings
- **`scMIXE_analysis.R`** – Using scMIXE (concatenated) embeddings
- **`marker_gene_filtering.R`** - Filters marker genes by log2FC, adjusted p-value, and expression threshold. Ensures each predicted cell type is represented by selecting top-ranked genes if none pass the filter.
- **`gene_set_enrichment_analysis.R`** - Performs GO enrichment analysis using Enrichr on filtered marker genes from each embedding strategy.

---

### `treatment outcome prediction/`
Contains scripts used to 



## Foundation Models Used
- [scGPT](https://github.com/bowang-lab/scGPT)
- [Mixedbread](https://www.mixedbread.com/docs/inference/embedding)

---

## Environment  
Code was developed and run in:
- R (version ≥ 4.2.0)
- Python 3.10+
- Google Colab Pro (L4 GPU for embedding extraction)

--- 

## Acknowledgements
Thank you to Dr. Pingzhao Hu for his guidance and support throughout this project. 


