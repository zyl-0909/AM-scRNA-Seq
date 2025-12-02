# ============================================================
#   cellchat_analysis.R
#
#   Description:
#     CellChat analyses including:
#       (1) Pathway heatmap for a selected condition
#       (2) RankNet comparison between two CellChat objects
#       (3) Circle network visualization for a selected pathway
#       (4) Infoflow scatter comparison
#
#   Inputs:
#     - data/cellchat_summary.xlsx
#     - data/cellchat_A.qs
#     - data/cellchat_B.qs
#
#   Outputs:
#     - cellchat_heatmap.pdf
#     - cellchat_ranknet.pdf
#     - cellchat_circle_compare.pdf
#     - cellchat_scatter.pdf
#
#   User-configurable parameters:
#     COND_COL, COND_VAL  : filtering for heatmap
#     CELL_COL, META_N    : cell-type / metadata setup
#     PATHWAY             : pathway for circle plot
#     SOURCE_IDX_A / B    : source node indices
#     SRC, TGT1, TGT2     : source & target sets for scatter
#
#   Notes:
#     - All cell labels are anonymized.
#     - Designed for GitHub publication.
#     - Recommended theme: theme_publication()
#
# ============================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(qs)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. CellChat Heatmap (Single Condition)
# ============================================================
df <- read.xlsx(file.path(DATA_DIR, "cellchat_summary.xlsx"), sheet = "Sheet1")

# ---- user configuration ----
COND_COL <- "condition"
COND_VAL <- "A"
CELL_COL <- "celltype"
META_N   <- 5

# ---- sanity checks ----
stopifnot(all(c(COND_COL, CELL_COL) %in% colnames(df)))
stopifnot(ncol(df) > META_N)

sub <- df %>% filter(.data[[COND_COL]] == COND_VAL)
stopifnot(nrow(sub) > 0)

pathway_cols <- colnames(sub)[(META_N + 1):ncol(sub)]

long_df <- sub %>%
  select(all_of(CELL_COL), all_of(pathway_cols)) %>%
  group_by(.data[[CELL_COL]]) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-all_of(CELL_COL), names_to = "pathway", values_to = "strength")

p_heat <- ggplot(long_df, aes(x = .data[[CELL_COL]], y = pathway, fill = strength)) +
  geom_tile() +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CellChat Pathway Heatmap", x = NULL, y = NULL, fill = "Strength")

ggsave(file.path(OUT_DIR, "cellchat_heatmap.pdf"), p_heat, width = 10, height = 8)


# ============================================================
# 2. RankNet Comparison
# ============================================================
obj_A <- qs::qread(file.path(DATA_DIR, "cellchat_A.qs"))
obj_B <- qs::qread(file.path(DATA_DIR, "cellchat_B.qs"))

merged <- mergeCellChat(list(A = obj_A, B = obj_B), add.names = TRUE)

p_rank1 <- rankNet(merged, mode = "comparison", stacked = TRUE,  do.stat = TRUE)
p_rank2 <- rankNet(merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)

ggsave(file.path(OUT_DIR, "cellchat_ranknet.pdf"),
       p_rank1 + p_rank2, width = 12, height = 5)


# ============================================================
# 3. Circle Plot Comparison
# ============================================================
PATHWAY <- "SPP1"
SOURCE_IDX_A <- 1
SOURCE_IDX_B <- c(1, 2)

nodes_A <- levels(obj_A@idents)
nodes_B <- levels(obj_B@idents)
nodes <- sort(intersect(nodes_A, nodes_B))

# sanity checks
stopifnot(length(nodes) > 1)
stopifnot(PATHWAY %in% dimnames(obj_A@netP$prob)[[3]])
stopifnot(PATHWAY %in% dimnames(obj_B@netP$prob)[[3]])
stopifnot(max(SOURCE_IDX_A, SOURCE_IDX_B) <= length(nodes))

anon_map <- setNames(sprintf("C%03d", seq_along(nodes)), nodes)
anon_nodes <- unname(anon_map[nodes])

matA <- obj_A@netP$prob[nodes, nodes, PATHWAY, drop = FALSE]
matB <- obj_B@netP$prob[nodes, nodes, PATHWAY, drop = FALSE]

rownames(matA) <- colnames(matA) <- anon_nodes
rownames(matB) <- colnames(matB) <- anon_nodes

col_use <- setNames(hue_pal()(length(anon_nodes)), anon_nodes)

p_circle_A <- netVisual_circle(
  matA,
  sources.use = anon_nodes[SOURCE_IDX_A],
  targets.use = anon_nodes,
  color.use   = col_use,
  title.name  = "Condition A"
)

p_circle_B <- netVisual_circle(
  matB,
  sources.use = anon_nodes[SOURCE_IDX_B],
  targets.use = anon_nodes,
  color.use   = col_use,
  title.name  = "Condition B"
)

pdf(file.path(OUT_DIR, "cellchat_circle_compare.pdf"), width = 6, height = 6)
p_circle_A
p_circle_B
dev.off()


# ============================================================
# 4. Infoflow Scatter Plots
# ============================================================
SRC   <- 1
TGT1  <- c(2, 3)
TGT2  <- c(4)
TOP_N <- 20

infoflow_scatter <- function(src_idx, tgt_idx, tag) {
  # shared nodes and pathways
  nodes  <- sort(intersect(levels(obj_A@idents), levels(obj_B@idents)))
  paths  <- intersect(dimnames(obj_A@netP$prob)[[3]], dimnames(obj_B@netP$prob)[[3]])

  stopifnot(src_idx <= length(nodes))
  stopifnot(all(tgt_idx <= length(nodes)))

  src_name <- nodes[src_idx]
  tgt_name <- nodes[tgt_idx]

  score_fun <- function(obj) {
    sapply(paths, function(pw) {
      m <- obj@netP$prob[, , pw]
      sum(m[src_name, tgt_name, drop = FALSE], na.rm = TRUE)
    })
  }

  sA <- score_fun(obj_A)
  sB <- score_fun(obj_B)

  # rank by total strength in A + B
  total <- sA + sB
  top   <- names(sort(total, decreasing = TRUE))[seq_len(min(TOP_N, length(total)))]

  df <- data.frame(
    Pathway = top,
    A = sA[top],
    B = sB[top]
  )

  ggplot(df, aes(x = A, y = B)) +
    geom_point(size = 3) +
    geom_text(aes(label = Pathway), vjust = -0.5, size = 3) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("Infoflow comparison", tag),
      x = "Condition A",
      y = "Condition B"
    )
}

p_info1 <- infoflow_scatter(SRC, TGT1, "1")
p_info2 <- infoflow_scatter(SRC, TGT2, "2")

ggsave(
  file.path(OUT_DIR, "cellchat_scatter.pdf"),
  p_info1 + p_info2,
  width = 10, height = 5
)

# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
