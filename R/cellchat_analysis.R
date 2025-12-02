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
  library(qs)
})

DATA_DIR <- "data"
OUT_DIR  <- "out"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- heatmap ----
df <- read.xlsx(file.path(DATA_DIR, "cellchat_summary.xlsx"), sheet = "Sheet1")
COND_COL <- "condition"
COND_VAL <- "A"
CELL_COL <- "celltype"
META_N   <- 5

sub <- df %>% filter(.data[[COND_COL]] == COND_VAL)
paths <- colnames(sub)[(META_N + 1):ncol(sub)]

long_df <- sub %>%
  select(all_of(CELL_COL), all_of(paths)) %>%
  group_by(.data[[CELL_COL]]) %>%
  summarise(across(.fns = sum), .groups = "drop") %>%
  pivot_longer(-all_of(CELL_COL))

p_heat <- ggplot(long_df, aes(.data[[CELL_COL]], name, fill = value)) +
  geom_tile()

ggsave(file.path(OUT_DIR, "cellchat_heatmap.pdf"), p_heat, width = 10, height = 8)

# ---- rankNet ----
A <- qs::qread(file.path(DATA_DIR, "cellchat_A.qs"))
B <- qs::qread(file.path(DATA_DIR, "cellchat_B.qs"))

merged <- mergeCellChat(list(A = A, B = B), add.names = TRUE)

p_rank1 <- rankNet(merged, mode = "comparison", stacked = TRUE)
p_rank2 <- rankNet(merged, mode = "comparison", stacked = FALSE)

ggsave(file.path(OUT_DIR, "cellchat_ranknet.pdf"), p_rank1 + p_rank2, width = 12, height = 5)

# ---- circle plot ----
PATHWAY <- "SPP1"
nodes <- sort(intersect(levels(A@idents), levels(B@idents)))
matA <- A@netP$prob[nodes, nodes, PATHWAY]
matB <- B@netP$prob[nodes, nodes, PATHWAY]

pdf(file.path(OUT_DIR, "cellchat_circle.pdf"), width = 6, height = 6)
netVisual_circle(matA)
netVisual_circle(matB)
dev.off()

# ---- scatter ----
SRC <- 1
TGT <- 2

score_fun <- function(obj) {
  paths <- dimnames(obj@netP$prob)[[3]]
  sapply(paths, function(pw) obj@netP$prob[SRC, TGT, pw])
}

sA <- score_fun(A)
sB <- score_fun(B)

df_scatter <- data.frame(Pathway = names(sA), A = sA, B = sB)

p_scatter <- ggplot(df_scatter, aes(A, B)) + geom_point()

ggsave(file.path(OUT_DIR, "cellchat_scatter.pdf"), p_scatter, width = 6, height = 5)

# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
