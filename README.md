# soil-microbiome-milkweed
Soil Microbiome Comparison Across Milkweed Ranges

This repository contains code, analysis, and results from a comparative soil microbiome study of **Asclepias syriaca** and **Asclepias speciosa**, focused on their respective native ranges in South Dakota.

## 🌱 Project Overview

The aim of this study is to compare the **bulk soil microbiome communities** associated with A. syriaca and A. speciosa across two distinct geographic locations:
- **Outdoor Campus – Sioux Falls, SD** (A. syriaca home range)
- **Custer State Park – Custer, SD** (A. speciosa home range)

We used **16S rRNA (V4)** and **ITS1 region** amplicon sequencing to investigate whether each milkweed species may have a **soil microbiome advantage** within its home range.

Although originally intended to target the rhizosphere, the available samples consist of **bulk soil surrounding the root zones**.

## 🚧 Project Status: In Progress

This analysis is currently ongoing. Documentation, results, and scripts are being actively developed and refined.

---

## 🧪 Methods

- **Sample Type:** Bulk soil
- **Sequencing Platform:** Illumina NextSeq 2000
- **Amplicons:** 
  - 16S rRNA (V4) for bacteria and archaea
  - ITS1 for fungi
- **Processing Tool:** QIIME2
  - Denoising with DADA2
  - Taxonomic classification with pre-trained classifiers
  - Rarefaction, alpha and beta diversity metrics
  - Ordination using Principal Coordinates Analysis (PCoA)
  - Differential abundance and indicator species analysis

---

# 📁 Repository Structure

```
soil-microbiome-milkweed/
├── data/           # Metadata and mapping files
├── scripts/        # QIIME2 processing and R scripts
├── results/        # Plots, tables, and key outputs
├── docs/           # Project overview and supplemental info
└── README.md
```
---

## 📊 Analyses Performed

- Quality filtering and ASV generation (DADA2)
- Taxonomy assignment (16S and ITS)
- Alpha diversity (Shannon, Observed ASVs, Faith's PD)
- Beta diversity (PCoA with Bray-Curtis and UNIFRAC distances)
- PERMANOVA tests
- Indicator species analysis
- Differential abundance testing (DESeq2, ANCOM)

---

## 📎 Example Outputs

Plots and tables showcasing community composition, diversity differences, and indicator taxa will be included in the `results/` folder.

---

## 🤝 Acknowledgments

This work was supported by the **South Dakota Biomedical Research Infrastructure Network (SD-BRIN)**.

Special thanks to:
- **Dr. Steven Matzner**
- Student researchers: **Jordan Hastad, Jack Erickson, and Pedro Da Rocha Borin**

---

## 📜 License

This project is released for educational and collaborative use. Please cite appropriately.

---

Created by **Kylene Guse**, 2025

