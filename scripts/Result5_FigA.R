############################################################
# Result5_FigA.R
#
# Sleep frequency–associated genes
# Cochran test lollipop plot
#
# Input:
#   data/Result5/Sleep_frequency_Cochran.tsv
#
# Output:
#   figures/Result5/Figure5A_Lollipop.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------

cochran_df <- fread(
  "data/Result5/Sleep_frequency_Cochran.tsv"
) %>% as.data.frame()

# ----------------------------------------------------------
# 2. Data processing
# ----------------------------------------------------------

cochran_df$log_Pvalue <- -log10(cochran_df$P_Value)

threshold <- -log10(0.05)

cochran_df <- cochran_df %>%
  mutate(
    Significant =
      ifelse(log_Pvalue > threshold,
             "Sleep frequency related",
             "Not significant")
  ) %>%
  arrange(desc(log_Pvalue)) %>%
  slice(1:32)

# ----------------------------------------------------------
# 3. Plot
# ----------------------------------------------------------

p1 <- ggplot(
  cochran_df,
  aes(x = reorder(Gene, -log_Pvalue),
      y = log_Pvalue,
      color = Significant)
) +
  geom_segment(
    aes(xend = reorder(Gene, -log_Pvalue),
        y = 0,
        yend = log_Pvalue),
    linewidth = 5
  ) +
  geom_point(size = 5) +
  geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    color = "red",
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Sleep frequency related" = "#ffa600",
      "Not significant" = "grey85"
    )
  ) +
  labs(
    x = NULL,
    y = expression(-log[10](P)),
    color = "Significance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x =
      element_text(angle = 90,
                   hjust = 1,
                   vjust = 1,
                   size = 12),
    axis.line =
      element_line(linewidth = 0.8)
  )

# ----------------------------------------------------------
# 4. Save
# ----------------------------------------------------------

OUT_DIR <- "figures/Result5/"
dir.create(OUT_DIR,
           recursive = TRUE,
           showWarnings = FALSE)

ggsave(
  file.path(OUT_DIR,
            "Figure5A_Lollipop.pdf"),
  p1,
  width = 8,
  height = 5
)

ggsave(
  file.path(OUT_DIR,
            "Figure5A_Lollipop.jpg"),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)
