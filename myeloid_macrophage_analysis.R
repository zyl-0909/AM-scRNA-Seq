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

# Myeloid UMAP + macrophage polarization + ROE + DEG barplot

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(Startrac)
  library(tidyr)
  library(openxlsx)
  library(ggbreak)
  library(qs)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1) Myeloid UMAP
obj <- qs::qread(file.path(DATA_DIR, "myeloid_object.qs"))

REDUCTION  <- "umap"
LABEL_COL  <- "celltype"
SPLIT_COL1 <- "condition"
SPLIT_COL2 <- "group"

make_cols <- function(x) {
  lv <- sort(unique(x))
  cols <- colorRampPalette(pal_npg("nrc")(10))(length(lv))
  setNames(cols, lv)
}
label_cols <- make_cols(obj@meta.data[[LABEL_COL]])

p0 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL, cols = label_cols)
p1 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL, cols = label_cols, split.by = SPLIT_COL1)
p2 <- DimPlot(obj, reduction = REDUCTION, group.by = LABEL_COL, cols = label_cols, split.by = SPLIT_COL2)

ggsave(file.path(OUT_DIR, "myeloid_umap.pdf"),        p0, width = 6,  height = 5)
ggsave(file.path(OUT_DIR, "myeloid_umap_split1.pdf"), p1, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "myeloid_umap_split2.pdf"), p2, width = 10, height = 5)

# 2) Macrophage polarization + ROE
obj2  <- qs::qread(file.path(DATA_DIR, "seurat_subset.qs"))
genes <- read.csv(file.path(DATA_DIR, "gene_sets.csv"))

setA <- genes$gene[genes$set == "A"]
setB <- genes$gene[genes$set == "B"]

obj2 <- AddModuleScore(obj2, features = list(setA, setB), name = c("A", "B"))
a_col <- tail(grep("^A", colnames(obj2@meta.data), value = TRUE), 1)
b_col <- tail(grep("^B", colnames(obj2@meta.data), value = TRUE), 1)

obj2$State <- ifelse(obj2@meta.data[[a_col]] > obj2@meta.data[[b_col]], "M1", "M2")

CLUSTER_COL <- "cluster"
SAMPLE_COL  <- "sample"

Roe <- calTissueDist(
  as.data.frame(obj2@meta.data),
  byPatient       = FALSE,
  colname.cluster = CLUSTER_COL,
  colname.patient = SAMPLE_COL,
  colname.tissue  = "State",
  method          = "chisq"
)

df <- as.data.frame(as.matrix(Roe)) |>
  tibble::rownames_to_column("Cluster") |>
  pivot_longer(-Cluster, names_to = "State", values_to = "R")

clv <- sort(unique(df$Cluster))
mp  <- setNames(sprintf("C%03d", seq_along(clv)), clv)
df$Cluster <- unname(mp[df$Cluster])

df$Symbol <- dplyr::case_when(
  df$R > 1.5 ~ "+++",
  df$R > 1   ~ "++",
  df$R >= .5 ~ "+",
  df$R > 0   ~ "+/-",
  TRUE       ~ "-"
)

pR <- ggplot(df, aes(State, Cluster, fill = R)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Symbol), size = 3) +
  scale_fill_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027", midpoint = 1) +
  theme_classic()

ggsave(file.path(OUT_DIR, "roe_heatmap.pdf"), pR, width = 3, height = 2.2)

# 3) DEG barplot
deg <- read.xlsx(file.path(DATA_DIR, "deg_summary.xlsx"))

CELL_COL <- "celltype"
UP_COL   <- "up_genes"
DOWN_COL <- "down_genes"
ANON_CELLTYPE <- TRUE

deg_long <- deg %>%
  transmute(
    Cell = .data[[CELL_COL]],
    up   = .data[[UP_COL]],
    down = .data[[DOWN_COL]]
  ) %>%
  pivot_longer(c(up, down), names_to = "Direction", values_to = "Count")

if (ANON_CELLTYPE) {
  lv <- sort(unique(deg_long$Cell))
  mp <- setNames(sprintf("C%03d", seq_along(lv)), lv)
  deg_long$Cell <- unname(mp[deg_long$Cell])
}

order_df <- deg_long %>%
  group_by(Cell) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  arrange(desc(total))

deg_long$Cell <- factor(deg_long$Cell, levels = order_df$Cell)

p_deg <- ggplot(deg_long, aes(Cell, Count, fill = Direction)) +
  geom_col(position = "dodge", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, "deg_barplot.pdf"), p_deg, width = 8, height = 5)
