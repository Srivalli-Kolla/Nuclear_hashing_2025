# Nuclear_hashing_2025 🧬🔬

**Single-nucleus hashing experiment for separating nuclei using sample-specific hashtags**

## 📌 Overview

This repository contains the analysis pipeline and results for the **nuclear hashing experiment**. The experiment explores the feasibility of demultiplexing pooled single-nucleus data based on sample-specific hashtags. Two experimental groups were used:

- **MCMV**: Male mice infected with murine cytomegalovirus
- **Non-infectious**: Female control mice

The nuclei were tagged with distinct hashtags and pooled together for sequencing. The central goal was to verify whether nuclei can be accurately separated by their hashtags and further analyzed for biological insights.

---

## 🎯 Project Aim

- **Identify the hashtag associated with each nucleus**
- **Filter out hashtag doublets** and retain high-confidence singlets
- **Perform quality control** and biological doublet detection
- **Annotate cell types** and define **refined clusters** for downstream analysis

---

## 🧬 Bioinformatics Workflow

1. **Hashtag demultiplexing**  
   Performed using [`HTODemux`] to:
   - Assign hashtags to nuclei
   - Filter out hashtag doublets
   - Retain only hashtag singlets

2. **Basic QC & Doublet Detection**  
   - Quality control filtering (e.g., low counts, high mitochondrial content)
   - Detection of potential **biological doublets** using [`Scrublet`]

3. **Annotation and Clustering**
   - Dimensionality reduction (PCA, UMAP)
   - Clustering and **cell type annotation**
   - Refined sub-cluster identification

---