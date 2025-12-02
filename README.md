# AM-scRNA-Seq

## Plotting scripts (Seurat + CellChat)

This repository contains R scripts to generate figures from precomputed **Seurat** and **CellChat** objects.
Input data files are expected to be prepared locally and are **not included** in this repository.

## Repository structure

```text
.
├── R/               # scripts
├── data/            # local inputs (NOT tracked by git)
└── out/             # local outputs (NOT tracked by git)
````

## Requirements

* R (>= 4.2 recommended)

### R packages used

Depending on which scripts you run:

* Core: `Seurat`, `qs`, `dplyr`, `ggplot2`, `patchwork`, `ggsci`
* Tables: `openxlsx`, `tidyr`
* CellChat: `CellChat`
* Optional:

  * `SCP` (only needed if you want `CellStatPlot`)
  * `Startrac` (only needed for the Roe/o-e heatmap)
  * `ggbreak` (only needed for DEG bar plot with y-axis break)

Example install (CRAN packages):

```r
install.packages(c("qs","dplyr","ggplot2","patchwork","ggsci","openxlsx","tidyr","ggbreak"))
```

Notes:

* Some packages (e.g. Seurat / CellChat / Startrac) may require additional installation steps depending on your system.

## Local inputs (place under `data/`)

These scripts expect local files under `data/` with neutral filenames. Example set:

* `data/seurat_object.qs` (Seurat object for global plots)
* `data/myeloid_object.qs` (Seurat object for myeloid-only plots)
* `data/cellchat_A.qs` and `data/cellchat_B.qs` (two CellChat objects for comparison)
* `data/cellchat_summary.xlsx` (summary table for a simple pathway heatmap)
* `data/deg_summary.xlsx` (DEG counts table for bar plot)
* `data/seurat_subset.qs` (subset Seurat object for module score + Roe plot)
* `data/gene_sets.csv` (optional/private gene sets file; recommended to keep untracked)

### `gene_sets.csv` format (example)

```csv
set,gene
A,GENE1
A,GENE2
B,GENE3
B,GENE4
```

## Outputs

All scripts write figures into `out/` (PDF). Example outputs:

* `out/figure1_umap_split_1.pdf`
* `out/figure1_umap_split_2.pdf`
* `out/figure1_dotplot_markers.pdf`
* `out/figure1_myeloid_umap.pdf`
* `out/cellchat_heatmap_single_group.pdf`
* `out/cellchat_ranknet_comparison.pdf`
* `out/circle_network_comparison.pdf`
* `out/infoflow_scatter_top_pathways.pdf`
* `out/roe_heatmap_ggplot.pdf`
* `out/deg_barplot.pdf`

Exact filenames depend on which scripts you run.

## Scripts and how to run

### 1) Seurat plots (UMAP / DotPlot / optional cell stats)

Script: `R/figure1_plots.R`

Config fields:

* `REDUCTION`: reduction name to plot
* `LABEL_COL`: metadata column used as labels
* `SPLIT_COL1`, `SPLIT_COL2`: metadata columns used for splitting panels
* (optional) `GROUP_BY` if using `SCP::CellStatPlot`

Run:

```r
source("R/figure1_plots.R")
```

### 2) CellChat heatmap from a summary table (single group)

Typically included in `R/figure1_plots.R` or separated as its own script.

Config fields:

* `COND_COL`: condition/group column name
* `COND_VAL`: a single group value to plot
* `CELL_COL`: cell-type column name
* `META_N`: number of metadata columns before pathway columns

Run:

```r
source("R/figure1_plots.R")
```

### 3) CellChat RankNet comparison

Typically included in `R/figure1_plots.R` or separated as `R/cellchat_ranknet.R`.

Run:

```r
source("R/figure1_plots.R")
```

### 4) CellChat circle plot comparison

Script: `R/cellchat_circle_compare.R`

Config fields:

* `PATHWAY`: pathway name to visualize
* `SOURCE_IDX_A`: source node index set for object A
* `SOURCE_IDX_B`: source node index set for object B

Run:

```r
source("R/cellchat_circle_compare.R")
```

### 5) CellChat infoflow scatter (top pathways by strength)

Script: `R/cellchat_infoflow_scatter.R`

Config fields:

* `SRC`: source node index
* `TGT_1`, `TGT_2`: target node indices for two panels
* `TOP_N`: number of top pathways to select

Run:

```r
source("R/cellchat_infoflow_scatter.R")
```

## Privacy / Reproducibility notes

* This repository does not include raw data or large binary objects.
* All local inputs go under `data/` and are ignored by git.
* All outputs go under `out/` and are ignored by git.
* Scripts avoid embedding private paths and use configurable metadata column names and index-based selection where possible.



就非常稳了。
```
