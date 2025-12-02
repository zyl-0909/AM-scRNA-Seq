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

obj <- qs::qread(file.path(DATA_DIR, "seurat_object.qs"), nthreads = 4)

# ==== user config ====
REDUCTION <- "umap"
LABEL_COL <- "celltype"
SPLIT_COL1 <- "group"
SPLIT_COL2 <- "sample_type"

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

stopifnot(REDUCTION %in% Reductions(obj))
stopifnot(all(c(LABEL_COL, SPLIT_COL1, SPLIT_COL2) %in% colnames(obj@meta.data)))

# colors
make_cols <- function(x) {
  lv <- sort(unique(x))
  cols <- colorRampPalette(pal_npg("nrc")(10))(length(lv))
  setNames(cols, lv)
}
label_cols <- make_cols(obj@meta.data[[LABEL_COL]])

# ==== UMAP ====
p_umap_1 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL,
                    cols = label_cols, split.by = SPLIT_COL1)
p_umap_2 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL,
                    cols = label_cols, split.by = SPLIT_COL2)

ggsave(file.path(OUT_DIR, "umap_split_group.pdf"),  p_umap_1, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "umap_split_sample_type.pdf"), p_umap_2, width = 10, height = 5)

# ==== DotPlot ====
Idents(obj) <- LABEL_COL
p_dot <- DotPlot(obj, features = rev(MARKERS)) +
  scale_color_gradient(low = "white", high = "#C93430") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "dotplot_markers.pdf"),
       p_dot, width = 8, height = 6)
     
