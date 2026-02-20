############################################################
# Result3_Figure3C_phylogenetic_signal.R
#
# Visualization of phylogenetic signal (Blomberg's K and
# Moran's I) for NREM ratio–associated circadian genes.
#
# Input:
#   data/Result3/NREM_ratio_phylogeneticsignal.tsv
#   data/Result3/NREM_ratio_Anova_Result.tsv
#
# Output:
#   figures/Result3/Figure3C/Phylogenetic_Signal_barplot.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

## =========================================================
## 1. Libraries
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

SIGNAL_FILE <- file.path(PROJECT_DIR,
                         "data", "Result3",
                         "NREM_ratio_phylogeneticsignal.tsv")

PERM_FILE <- file.path(PROJECT_DIR,
                       "data", "Result3",
                       "NREM_ratio_Anova_Result.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures", "Result3",
                     "Figure3C")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Load data
## =========================================================

DF <- fread(SIGNAL_FILE) |> as.data.frame()
perm_test <- fread(PERM_FILE) |> as.data.frame()

sig_genes <- perm_test %>%
  filter(adj.P < 0.05)

## Keep only significant genes from enrichment test
DF <- DF %>%
  filter(Gene %in% sig_genes$Gene)

## =========================================================
## 4. Filter by phylogenetic significance
## =========================================================

DF_filtered <- DF %>%
  filter(P.K < 0.05, P.I < 0.05) %>%   # both significant
  arrange(`Moran’s I`) %>%             # 낮은 I 값 순
  slice_head(n = 30)

if (nrow(DF_filtered) == 0) {
  stop("No genes satisfy P.K < 0.05 and P.I < 0.05.")
}

## Rename for plotting
colnames(DF_filtered)[colnames(DF_filtered) == "K"] <- "Blomberg’s K"
colnames(DF_filtered)[colnames(DF_filtered) == "I"] <- "Moran’s I"

## Preserve ordering
DF_filtered$Gene <- factor(
  DF_filtered$Gene,
  levels = DF_filtered$Gene
)

## =========================================================
## 5. Long format
## =========================================================

df_long <- DF_filtered %>%
  select(Gene, `Blomberg’s K`, `Moran’s I`) %>%
  pivot_longer(
    cols = c(`Blomberg’s K`, `Moran’s I`),
    names_to = "Signal",
    values_to = "Value"
  )

## =========================================================
## 6. Plot
## =========================================================

Phylogenetic_Signal_plot <- ggplot(
  df_long,
  aes(x = Gene, y = Value, color = Signal)
) +
  geom_line(aes(group = Signal),
            linewidth = 1.2,
            alpha = 0.4,
            color = "grey60") +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Blomberg’s K" = "#7593af",
    "Moran’s I"    = "#730220"
  )) +
  theme_bw(base_size = 12) +
  labs(
    x = "",
    y = "Phylogenetic signal score",
    color = "Signal"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

## =========================================================
## 7. Save
## =========================================================

ggsave(file.path(OUT_DIR,
                 "Figure3C_Phylogenetic_Signal.pdf"),
       Phylogenetic_Signal_plot,
       width = 6, height = 3.5)

ggsave(file.path(OUT_DIR,
                 "Figure3C_Phylogenetic_Signal.jpg"),
       Phylogenetic_Signal_plot,
       width = 6, height = 3.5,
       dpi = 300)

message("Figure3C saved successfully.")
