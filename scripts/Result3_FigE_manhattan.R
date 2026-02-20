############################################################
# Result3_Figure3E_NREM_manhattan.R
#
# Manhattan plot for NREM ratio–associated SNPs
#
# Input:
#   data/Result3/NREM_ratio_SNP.tsv
#
# Output:
#   figures/Result3/Figure3E/Manhattan_plot.pdf
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
  library(ggplot2)
  library(ggrepel)
})

## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

SNP_FILE <- file.path(PROJECT_DIR,
                      "data", "Result3",
                      "NREM_ratio_SNP.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures", "Result3",
                     "Figure3E")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Load data
## =========================================================

TS_site1 <- fread(SNP_FILE) |> as.data.table()

if (!all(c("CHR", "BP", "P") %in% colnames(TS_site1))) {
  stop("Input file must contain CHR, BP, and P columns.")
}

TS_site1[, CHR := as.numeric(CHR)]

## =========================================================
## 4. Cumulative genomic position
## =========================================================

don <- TS_site1 %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  left_join(TS_site1, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(
    BPcum = BP + tot
  )

## X-axis chromosome centers
axisdf <- don %>%
  group_by(CHR) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2,
            .groups = "drop")

## =========================================================
## 5. Chromosome colors (auto alternating)
## =========================================================

n_chr <- length(unique(don$CHR))

chr_colors <- rep(c("#1e3d58", "#43b0f1"),
                  length.out = n_chr)

names(chr_colors) <- sort(unique(don$CHR))

## =========================================================
## 6. Manhattan plot
## =========================================================

P_threshold <- 0.05

P1 <- ggplot(don, aes(x = BPcum, y = -log10(P))) +
  geom_point(aes(color = factor(CHR)),
             alpha = 0.8,
             size = 1.2) +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(
    label = axisdf$CHR,
    breaks = axisdf$center
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(-log10(don$P), na.rm = TRUE) + 1)
  ) +
  geom_hline(yintercept = -log10(P_threshold),
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  geom_text_repel(
    data = subset(don,
                  P < P_threshold &
                  Gene %in% c("ARNTL2", "PPARA")),
    aes(label = SNP),
    size = 3,
    max.overlaps = 10
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  )

## =========================================================
## 7. Save
## =========================================================

ggsave(file.path(OUT_DIR,
                 "Figure3E_NREM_ratio_Manhattan.pdf"),
       P1,
       width = 14, height = 4)

ggsave(file.path(OUT_DIR,
                 "Figure3E_NREM_ratio_Manhattan.jpg"),
       P1,
       width = 14, height = 4,
       dpi = 300)

message("Manhattan plot saved successfully.")
