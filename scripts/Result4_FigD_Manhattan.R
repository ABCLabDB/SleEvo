############################################################
# Result4_Fig4D_Manhattan.R
#
# Manhattan plot for sleep timing association
#
# Input:
#   data/Result4/Sleeptiming_Manhattan_Dataset.tsv
#
# Output:
#   figures/Result4/Figure4D_Manhattan.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

PROJECT_DIR <- getwd()

## =========================================================
## 1. Paths
## =========================================================

DATA_FILE <- file.path(PROJECT_DIR,
                       "data","Result4",
                       "Sleeptiming_Manhattan_Dataset.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures","Result4")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. Load data
## =========================================================

TS <- fread(DATA_FILE) |> as.data.frame()

TS <- TS %>%
  mutate(
    CHR = as.numeric(CHR),
    FDR = P,
    Label = ifelse(-log10(P) > 3, SNP, NA_character_)
  )

## =========================================================
## 3. Cumulative genomic position
## =========================================================

chr_info <- TS %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(offset = cumsum(chr_len) - chr_len)

TS_cum <- TS %>%
  left_join(chr_info, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + offset)

axisdf <- TS_cum %>%
  group_by(CHR) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2)

## =========================================================
## 4. Chromosome colors
## =========================================================

base_colors <- c("#1e3d58","#43b0f1")
chr_levels <- sort(unique(TS_cum$CHR))

chr_colors <- rep(base_colors,
                  length.out = length(chr_levels))

names(chr_colors) <- chr_levels

## =========================================================
## 5. Plot
## =========================================================

highlight_genes <- c("CAVIN3","ATF5","PER1","CRY2")

p <- ggplot(TS_cum,
            aes(x = BPcum,
                y = -log10(P))) +

  geom_point(aes(color = factor(CHR)),
             alpha = 0.8,
             size = 1.5) +

  scale_color_manual(values = chr_colors) +

  scale_x_continuous(
    breaks = axisdf$center,
    labels = axisdf$CHR
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +

  geom_text_repel(
    data = subset(TS_cum,
                  -log10(P) > -log10(0.05) &
                  Gene %in% highlight_genes),
    aes(label = SNP),
    size = 3.5,
    max.overlaps = 15
  ) +

  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_text(face = "bold")
  ) +

  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  )

## =========================================================
## 6. Save
## =========================================================

ggsave(file.path(OUT_DIR,
                 "Figure4D_Manhattan.pdf"),
       p,
       width = 9,
       height = 4)

ggsave(file.path(OUT_DIR,
                 "Figure4D_Manhattan.jpg"),
       p,
       width = 9,
       height = 4,
       dpi = 300)

message("Figure4D Manhattan plot saved.")
