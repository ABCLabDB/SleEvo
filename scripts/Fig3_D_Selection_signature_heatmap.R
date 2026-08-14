############################################################
# Fig3_D_Selection_signature_heatmap.R
#
# Figure 3D. Selection signature scatter plot
#   - x: dN/dS, y: Tajima's D
#   - Four quadrants labeled as in the manuscript figure:
#       Neutral drift / Balancing selection /
#       Purifying selection / Positive selection
#
# Input:
#   - data/Result2/Total_sleep_time_Anova_Result.tsv
#   - data/Heatmap/Total_sleep_time.tsv
#
# Output:
#   - Result/Fig3/Fig3_D_Selection_signature.jpg
#   - Result/Fig3/Fig3_D_Selection_signature.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


## =========================================================
## 1. Libraries
## =========================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

ANOVA_FILE <- file.path(PROJECT_DIR, "data", "Result2",
                        "Total_sleep_time_Anova_Result.tsv")
HEATMAP_FILE <- file.path(PROJECT_DIR, "data", "Heatmap",
                          "Total_sleep_time.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Load data
## =========================================================

sleep_candidate <- fread(ANOVA_FILE) |> as.data.frame()
sleep_candidate <- sleep_candidate[sleep_candidate$adj.P < 0.05, ]
sleep_gene <- sleep_candidate$Gene

evolution_score_DF <- fread(HEATMAP_FILE) |> as.data.frame()
evolution_score_DF <- evolution_score_DF[
  evolution_score_DF$Gene %in% sleep_gene, ]
rownames(evolution_score_DF) <- NULL

evolution_score_DF1 <- evolution_score_DF %>%
  mutate(
    dN_dS_Group = ifelse(dN_dS > 1, "dN/dS > 1", "dN/dS <= 1"),
    TajimasD_Group = ifelse(TajimasD > 0, "Tajima's D > 0", "Tajima's D <= 0")
  )

highlight_genes <- c("ARNTL2", "CPT1A", "ATF4", "PPARA", "CRY1")


## =========================================================
## 4. Plot (panel D style with quadrant labels)
## =========================================================

p1 <- ggplot(evolution_score_DF1, aes(x = dN_dS, y = TajimasD)) +

  ## quadrant backgrounds
  annotate("rect", xmin = -Inf, xmax = 1,   ymin = 0,    ymax = Inf,
           fill = "#ab5852", alpha = 0.1) +   # Neutral drift
  annotate("rect", xmin = 1,   xmax = Inf, ymin = 0,    ymax = Inf,
           fill = "#7593af", alpha = 0.1) +   # Balancing selection
  annotate("rect", xmin = -Inf, xmax = 1,   ymin = -Inf, ymax = 0,
           fill = "#d69e49", alpha = 0.1) +   # Purifying selection
  annotate("rect", xmin = 1,   xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#eadaa0", alpha = 0.1) +   # Positive selection

  ## quadrant labels (Fig3 D style)
  annotate("text", x = -Inf, y = Inf,
           label = "Neutral drift",
           hjust = -0.05, vjust = 1.5,
           color = "#8a3f3a", size = 4.2, fontface = "bold") +
  annotate("text", x = Inf, y = Inf,
           label = "Balancing selection",
           hjust = 1.05, vjust = 1.5,
           color = "#4f6b86", size = 4.2, fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf,
           label = "Purifying selection",
           hjust = -0.05, vjust = -0.8,
           color = "#a87728", size = 4.2, fontface = "bold") +
  annotate("text", x = Inf, y = -Inf,
           label = "Positive selection",
           hjust = 1.05, vjust = -0.8,
           color = "#b39d4a", size = 4.2, fontface = "bold") +

  ## non-highlighted genes (semi-transparent)
  geom_point(
    data = evolution_score_DF1[!evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(color = interaction(dN_dS_Group, TajimasD_Group)),
    size = 6.5,
    alpha = 0.5
  ) +

  ## highlighted genes (fully opaque)
  geom_point(
    data = evolution_score_DF1[evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(color = interaction(dN_dS_Group, TajimasD_Group)),
    size = 6.5,
    alpha = 1
  ) +

  geom_text_repel(
    data = evolution_score_DF1[!evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(label = Gene),
    size = 4,
    max.overlaps = 2,
    color = "black"
  ) +

  geom_text_repel(
    data = evolution_score_DF1[evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(label = Gene),
    size = 4.5,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf
  ) +

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60",
             linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60",
             linewidth = 0.6) +

  scale_color_manual(
    values = c(
      "dN/dS > 1.Tajima's D > 0"  = "#7593af",
      "dN/dS > 1.Tajima's D <= 0" = "#eadaa0",
      "dN/dS <= 1.Tajima's D > 0" = "#ab5852",
      "dN/dS <= 1.Tajima's D <= 0"= "#d69e49"
    ),
    name = "Type"
  ) +
  labs(
    title = "",
    x = "dN/dS",
    y = "Tajima's D"
  ) +
  theme_minimal() +
  theme(
    title = element_text(size = 16),
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )


## =========================================================
## 5. Save
## =========================================================

ggsave(file.path(OUT_DIR, "Fig3_D_Selection_signature.jpg"),
       plot = p1, width = 10, height = 6, units = "in", dpi = 300)
ggsave(file.path(OUT_DIR, "Fig3_D_Selection_signature.pdf"),
       plot = p1, width = 10, height = 6)

message("Saved: Result/Fig3/Fig3_D_Selection_signature.jpg|.pdf")
