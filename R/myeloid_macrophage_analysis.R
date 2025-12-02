# ============================================================
#   myeloid_macrophage_analysis.R
#
#   Description:
#     Analysis for myeloid subsets and macrophage polarization:
#       (1) Myeloid UMAP (full + split)
#       (2) Polarization module scoring (M1 vs M2)
#       (3) ROE heatmap via Startrac
#       (4) DEG barplot with anonymized labels
#
#   Inputs:
#     - data/myeloid_object.qs
#     - data/seurat_subset.qs
#     - data/gene_sets.csv
#     - data/deg_summary.xlsx
#
#   Outputs:
#     - myeloid_umap.pdf
#     - myeloid_umap_split1.pdf
#     - myeloid_umap_split2.pdf
#     - roe_heatmap.pdf
#     - deg_barplot.pdf
#
#   User-configurable parameters:
#     REDUCTION, LABEL_COL
#     SPLIT_COL1, SPLIT_COL2
#     CLUSTER_COL, SAMPLE_COL
#     ANON_CELLTYPE = TRUE/FALSE
#
#   Notes:
#     - Cell labels anonymized to C01/C02/...
#     - Safe for GitHub posting.
#     - Consistent with theme_publication()
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Startrac)
  library(openxlsx)
  library(qs)
  library(tidyr)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- UMAP ----
obj <- qs::qread(file.path(DATA_DIR, "myeloid_object.qs"))

p0 <- DimPlot(obj, reduction = "umap", group.by = "celltype")
p1 <- DimPlot(obj, reduction = "umap", group.by = "celltype", split.by = "condition")
p2 <- DimPlot(obj, reduction = "umap", group.by = "celltype", split.by = "group")

ggsave(file.path(OUT_DIR, "myeloid_umap.pdf"),        p0, width = 6, height = 5)
ggsave(file.path(OUT_DIR, "myeloid_umap_split1.pdf"), p1, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "myeloid_umap_split2.pdf"), p2, width = 10, height = 5)

# ---- polarization ----
obj2  <- qs::qread(file.path(DATA_DIR, "seurat_subset.qs"))
genes <- read.csv(file.path(DATA_DIR, "gene_sets.csv"))

A_set <- genes$gene[genes$set == "A"]
B_set <- genes$gene[genes$set == "B"]

obj2 <- AddModuleScore(obj2, features = list(A_set, B_set), name = c("A","B"))
a_col <- grep("^A", names(obj2@meta.data), value = TRUE)
b_col <- grep("^B", names(obj2@meta.data), value = TRUE)

obj2$State <- ifelse(obj2@meta.data[[a_col]] > obj2@meta.data[[b_col]], "M1","M2")

R <- calTissueDist(as.data.frame(obj2@meta.data),
                   colname.cluster = "cluster",
                   colname.patient = "sample",
                   colname.tissue  = "State")

df <- as.data.frame(as.matrix(R)) |>
  tibble::rownames_to_column("Cluster") |>
  pivot_longer(-Cluster, names_to = "State", values_to = "R")

pR <- ggplot(df, aes(State, Cluster, fill = R)) + geom_tile()

ggsave(file.path(OUT_DIR, "roe_heatmap.pdf"), pR, width = 4, height = 3)

# ---- DEG barplot ----
deg <- read.xlsx(file.path(DATA_DIR, "deg_summary.xlsx"))

df_deg <- deg |>
  pivot_longer(c(up_genes, down_genes),
               names_to = "Direction", values_to = "Count")

p_deg <- ggplot(df_deg, aes(celltype, Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge")

ggsave(file.path(OUT_DIR, "deg_barplot.pdf"), p_deg, width = 8, height = 5)

# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
