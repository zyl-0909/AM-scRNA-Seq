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
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
  library(qs)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------
obj <- qs::qread(file.path(DATA_DIR, "seurat_object.qs"), nthreads = 4)

# ------------------------------------------------------------
# User-configurable columns
# ------------------------------------------------------------
REDUCTION  <- "umap"          # any available reduction in the object
LABEL_COL  <- "celltype"      # cell annotation column
SPLIT_COL1 <- "group"         # first split column
SPLIT_COL2 <- "sample_type"   # second split column

# Marker genes (safe)
MARKERS <- c(
  "ACTA2","TAGLN","TPM2",
  "PMEL","TYRP1","MLANA",
  "LUM","COL1A2","COL1A1","COL3A1",
  "KRT5","KRT14","KRT1",
  "LYZ","HLA-DRA","HLA-DPB1","CD14",
  "PECAM1","VWF","CLDN5",
  "PTPRC","CD3D","CD8A","CCL21",
  "TPSB2","CPA3",
  "S100B","MPZ","SOX10"
)

# ------------------------------------------------------------
# Small checks
# ------------------------------------------------------------
stopifnot(REDUCTION %in% Reductions(obj))
stopifnot(all(c(LABEL_COL, SPLIT_COL1, SPLIT_COL2) %in% colnames(obj@meta.data)))

# ------------------------------------------------------------
# Color palette
# ------------------------------------------------------------
make_cols <- function(x) {
  lv <- sort(unique(x))
  cols <- colorRampPalette(pal_npg("nrc")(10))(length(lv))
  setNames(cols, lv)
}
label_cols <- make_cols(obj@meta.data[[LABEL_COL]])

# ------------------------------------------------------------
# 1) UMAP: split-by group
# ------------------------------------------------------------
p_umap_1 <- DimPlot(
  obj, reduction = REDUCTION,
  group.by = LABEL_COL, cols = label_cols,
  label = FALSE, raster = FALSE, split.by = SPLIT_COL1
)

ggsave(file.path(OUT_DIR, "figure1_umap_split_group.pdf"),
       p_umap_1, width = 10, height = 5)

# ------------------------------------------------------------
# 2) UMAP: split-by sample_type
# ------------------------------------------------------------
p_umap_2 <- DimPlot(
  obj, reduction = REDUCTION,
  group.by = LABEL_COL, cols = label_cols,
  label = FALSE, raster = FALSE, split.by = SPLIT_COL2
)

ggsave(file.path(OUT_DIR, "figure1_umap_split_sample.pdf"),
       p_umap_2, width = 10, height = 5)

# ------------------------------------------------------------
# 3) DotPlot for markers 
# ------------------------------------------------------------
Idents(obj) <- LABEL_COL

p_dot <- DotPlot(obj, features = rev(MARKERS)) +
  scale_color_gradient(low = "white", high = "#C93430") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  labs(color = "Avg Expr", size = "Pct")

ggsave(file.path(OUT_DIR, "figure1_dotplot_markers.pdf"),
       p_dot, width = 8, height = 6)

# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
