############################################################
# Fig3_C_Phylogenetic_signal.R
#
# Figure 3C. Phylogenetic signal across genes
#   - Blomberg’s K and Moran’s I (line + point)
#   - Genes: alphabetical order, first 20 (Fig3 C style)
#   Extracted / restyled from Result2_visualization.R (Figure 2E)
#
# Input:
#   - data/Result2/Totalsleeptime_Phylogenetic_signal.tsv
#
# Output:
#   - Result/Fig3/Fig3_C_Phylogenetic_signal.jpg
#   - Result/Fig3/Fig3_C_Phylogenetic_signal.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


## =========================================================
## 1. Libraries
## =========================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

SIGNAL_FILE <- file.path(
  PROJECT_DIR, "data", "Result2",
  "Totalsleeptime_Phylogenetic_signal.tsv"
)

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Load & prepare signal data
## =========================================================

signal_df <- fread(SIGNAL_FILE) |> as.data.frame()
signal_df[is.na(signal_df)] <- 0
colnames(signal_df)[c(2, 6)] <- c("Blomberg’s K", "Moran’s I")

## Fig3 C: alphabetical sort, take 20 genes (no P-value filter)
signal_top20 <- signal_df %>%
  arrange(Gene) %>%
  slice_head(n = 20)

signal_long <- signal_top20 |>
  dplyr::select(Gene, `Blomberg’s K`, `Moran’s I`) |>
  tidyr::pivot_longer(
    cols = c(`Blomberg’s K`, `Moran’s I`),
    names_to = "Signal",
    values_to = "Value"
  )

signal_long$Gene <- factor(signal_long$Gene, levels = unique(signal_top20$Gene))
signal_long$Signal <- factor(
  signal_long$Signal,
  levels = c("Blomberg’s K", "Moran’s I")
)


## =========================================================
## 4. Plot (Fig3 C style)
## =========================================================

p_signal <- ggplot(signal_long,
                   aes(x = Gene, y = Value, color = Signal, group = Signal)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 3.2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Blomberg’s K" = "#7593af",
      "Moran’s I"    = "#d69e49"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(y = "Phylogenetic signal score", x = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    plot.margin = margin(6, 10, 6, 6)
  )


## =========================================================
## 5. Save
## =========================================================

ggsave(file.path(OUT_DIR, "Fig3_C_Phylogenetic_signal.jpg"),
       p_signal, width = 7, height = 3.5, units = "in", dpi = 300)
ggsave(file.path(OUT_DIR, "Fig3_C_Phylogenetic_signal.pdf"),
       p_signal, width = 7, height = 3.5)

message("Genes (n=", nrow(signal_top20), "): ",
        paste(signal_top20$Gene, collapse = ", "))
message("Saved: Result/Fig3/Fig3_C_Phylogenetic_signal.jpg|.pdf")
