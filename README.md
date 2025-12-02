# AM-scRNA-Seq Analysis

This repository contains reproducible R scripts for analyzing single-cell RNA-seq
data from airway macrophages (AM) and related immune subsets. The codebase is 
organized to generate publication-ready figures and to perform CellChat and 
polarization analyses on curated Seurat objects.

All scripts are **safe for public release**. No patient identifiers, sample IDs, 
or sensitive metadata are included.

---

##  Repository Structure

```

AM-scRNA-Seq/
├── R/
│   ├── figure1\_alldata\_plot.R
│   ├── cellchat\_analysis.R
│   ├── myeloid\_macrophage\_analysis.R
│
├── data/        # Input data (.qs, .csv, .xlsx); not tracked in Git
├── out/         # Output figures (auto-generated)
├── .gitignore
└── README.md

````

---

##  Requirements

The following R packages are needed:

- Seurat  
- dplyr  
- ggplot2  
- patchwork  
- ggsci  
- tidyr  
- openxlsx  
- CellChat  
- Startrac  
- ggbreak (optional)  
- qs  

Install missing packages using:

```r
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "ggsci",
  "tidyr", "openxlsx", "qs"
))

# From GitHub:
# remotes::install_github("sqjin/CellChat")
# remotes::install_github("wu-yc/Startrac")
````

---

##  Running the Scripts

Place all required input files (Seurat objects, CellChat objects, gene lists,
DEG tables, etc.) inside the `data/` directory.

Then run any script using:

```bash
Rscript R/figure1_alldata_plot.R
Rscript R/cellchat_analysis.R
Rscript R/myeloid_macrophage_analysis.R
```

All output figures will be saved automatically to `out/`.

---

##  Script Overview

### **1. `figure1_alldata_plot.R`**

Generates:

* UMAPs for the full dataset

  * split by selected metadata columns
* Marker DotPlot for selected gene sets

This script corresponds to the "All Data" overview figure.

---

### **2. `cellchat_analysis.R`**

Performs multiple CellChat analyses:

* Pathway heatmap for one condition
* RankNet comparison between two CellChat objects
* Circle plot comparison for a chosen pathway
* Infoflow scatter plot comparing outgoing communication strength

All cell identities are anonymized automatically.

---

### **3. `myeloid_macrophage_analysis.R`**

Contains:

* Myeloid subset UMAPs (full and split)
* Macrophage M1/M2 polarization (module scores)
* ROE heatmap using Startrac
* DEG barplot (labels anonymized when enabled)

Designed for downstream functional inspection of myeloid populations.

---

##  Figures

Example outputs include:

* UMAP distributions
* Marker DotPlots
* Communication networks
* Heatmaps (pathway and ROE)
* DEG barplots

All figures are saved to the `out/` folder.

---

##  Data Privacy

This repository:

* Does **not** contain raw sequencing data
* Does **not** expose sensitive metadata
* Contains only anonymized cell-type identifiers
* Uses generic filenames (e.g., `seurat_object.qs`, `cellchat_A.qs`)

Safe for public GitHub use.

