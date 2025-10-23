This repository contains all code, scripts, and results for the research project titled **PRECISE: Cell-Type Inference from Single-Cell Embeddings Enhances Prediction of Breast Cancer Immunotherapy Response.**  

---

## Repository Structure


### `embedding generation/`  
Contains all scripts used to generate embeddings

-
- 

---


### `downstream analysis/`  
Contains all scripts used to perform clustering, identify markers, and annotate cell types. 

- **`raw_data_analysis.R`** – Using raw data ([Bassez et al. data](https://lambrechtslab.sites.vib.be/en/single-cell))
- **`individual_embedding_analysis.R`** – Using scGPT or Mixedbread embeddings
- **`scMIXE_analysis.R`** – Using scMIXE (concatenated) embeddings

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


