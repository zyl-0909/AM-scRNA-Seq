# ============================================================
#   figure1_alldata_plot.R
#
#   Description:
#     UMAP visualizations and marker DotPlot for the full
#     single-cell dataset ("All-data").
#
#   Inputs:
#     - data/seurat_object.qs
#
#   Outputs:
#     - umap_split_group.pdf
#     - umap_split_sample.pdf
#     - dotplot_markers.pdf
#
#   User-configurable parameters:
#     REDUCTION   : reduction name (e.g., "umap")
#     LABEL_COL   : cell annotation column
#     SPLIT_COL1  : first facet column
#     SPLIT_COL2  : second facet column
#     MARKERS     : marker gene list
#
#   Notes:
#     - No original sample IDs or patient labels included.
#     - Safe for public GitHub release.
#     - Suggested to load: theme_publication()
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggsci)
  library(qs)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj <- qs::qread(file.path(DATA_DIR, "seurat_object.qs"))

REDUCTION  <- "umap"
LABEL_COL  <- "celltype"
SPLIT_COL1 <- "group"
SPLIT_COL2 <- "sample_type"

MARKERS <- c(
  "ACTA2","TAGLN","TPM2","PMEL","TYRP1","MLANA",
  "LUM","COL1A2","COL1A1","COL3A1","KRT5","KRT14","KRT1",
  "LYZ","HLA-DRA","HLA-DPB1","CD14",
  "PECAM1","VWF","CLDN5","PTPRC","CD3D","CD8A","CCL21",
  "TPSB2","CPA3","S100B","MPZ","SOX10"
)

p1 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL, split.by = SPLIT_COL1)
p2 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL, split.by = SPLIT_COL2)

ggsave(file.path(OUT_DIR, "umap_group.pdf"),  p1, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "umap_sample.pdf"), p2, width = 10, height = 5)

Idents(obj) <- LABEL_COL
p_dot <- DotPlot(obj, features = rev(MARKERS))

ggsave(file.path(OUT_DIR, "dotplot_markers.pdf"), p_dot, width = 8, height = 6)

# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
