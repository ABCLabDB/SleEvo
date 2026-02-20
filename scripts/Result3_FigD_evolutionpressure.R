############################################################
# Result3_Figure3D_NREM_heatmap.R
#
# Selection signature heatmap for
# NREM ratio–associated circadian genes.
#
# Input:
#   data/Heatmap/NREM_ratio.tsv
#   data/Result3/NREM_ratio_Anova_Result.tsv
#
# Output:
#   figures/Result3/Figure3D/NREM_ratio_Heatmap.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

## =========================================================
## 1. Libraries
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(Cairo)
  library(gridExtra)
})

## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

HEATMAP_FILE <- file.path(PROJECT_DIR,
                          "data", "Heatmap",
                          "NREM_ratio.tsv")

PERM_FILE <- file.path(PROJECT_DIR,
                       "data", "Result3",
                       "NREM_ratio_Anova_Result.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures", "Result3",
                     "Figure3D")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Load data
## =========================================================

heat_df <- fread(HEATMAP_FILE) |> as.data.frame()
perm_df <- fread(PERM_FILE)    |> as.data.frame()

sig_genes <- perm_df |>
  subset(adj.P < 0.05)

if (nrow(sig_genes) == 0) {
  stop("No significant genes (adj.P < 0.05)")
}

## =========================================================
## 4. Filter for significant genes
## =========================================================

heat_df <- heat_df[heat_df$Gene %in% sig_genes$Gene, ]

if (nrow(heat_df) == 0) {
  stop("No overlapping genes between heatmap and perm_test.")
}

heat_df <- heat_df[order(heat_df$Cluster, decreasing = TRUE), ]
heat_df <- heat_df[order(heat_df$Gene), ]

rownames(heat_df) <- heat_df$Gene

## =========================================================
## 5. Matrix construction
## =========================================================

mat <- heat_df[, c("Ps", "Pi", "TajimasD", "dN_dS")]
colnames(mat) <- c("Ps", "π", "Tajima's D", "dN/dS")

## Optional: top 30 genes only
mat <- head(mat, 30)

## =========================================================
## 6. Heatmaps
## =========================================================

heatmap_Tajima <- pheatmap(
  t(mat[, "Tajima's D", drop = FALSE]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#053061", "white", "#67001f"))(100),
  border_color = "gray80",
  main = "Tajima's D",
  silent = TRUE
)

heatmap_dNdS <- pheatmap(
  t(mat[, "dN/dS", drop = FALSE]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
  border_color = "gray80",
  main = "dN/dS",
  silent = TRUE
)

heatmap_pi <- pheatmap(
  t(mat[, "π", drop = FALSE]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#E67E22", "#B03A2E"))(100),
  border_color = "gray80",
  main = "Nucleotide Diversity (π)",
  silent = TRUE
)

## =========================================================
## 7. Save figure
## =========================================================

CairoPDF(file.path(OUT_DIR,
                   "Figure3D_NREM_ratio_Heatmap.pdf"),
         width = 24, height = 8)

grid.arrange(
  heatmap_Tajima[[4]],
  heatmap_dNdS[[4]],
  heatmap_pi[[4]],
  nrow = 3
)

dev.off()

message("Figure3D heatmap saved successfully.")
