############################################################
# Result6_visualization_DEG_boxplot.R
# Boxplot visualization for DEG genes (Result6 Figure B)
# Linear model: Expression ~ Group + Species:Sequencer
# Author: Sleep Evolution Project
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(broom)
library(stringr)

############################################################
# Base paths (GitHub structure)
############################################################

BASE_DATA  <- "data/Result6"
BASE_PLOT  <- "data/Result6/Figure6B_particle"

dir.create(BASE_PLOT, showWarnings = FALSE, recursive = TRUE)

############################################################
# General plotting function
############################################################

plot_deg_boxplot <- function(meta_file,
                             count_file,
                             group_column,
                             genes,
                             coef_name,
                             colors,
                             suffix) {

  # Load metadata
  meta <- fread(file.path(BASE_DATA, meta_file)) |> as.data.frame()
  rownames(meta) <- meta$Run

  # Load expression matrix
  count <- fread(file.path(BASE_DATA, count_file)) |> as.data.frame()
  rownames(count) <- count[,1]
  count <- count[,-1]

  # Loop through selected genes
  for (gene in genes) {

    if (!gene %in% rownames(count)) next

    Expression <- as.numeric(count[gene, ])

    df <- data.frame(
      Sample     = meta$Run,
      Group      = meta[[group_column]],
      Sequencer  = meta$Instrument,
      Species    = meta$Species,
      Expression = Expression
    )

    # Remove NA / zero-only cases
    df <- df[!is.na(df$Group), ]
    df <- df[!is.na(df$Expression), ]

    if (all(df$Expression == 0)) next

    # Linear model
    fit <- lm(Expression ~ Group + Species:Sequencer, data = df)
    lm_summary <- summary(fit)

    if (!coef_name %in% rownames(lm_summary$coefficients)) next

    coef_info <- lm_summary$coefficients[coef_name, ]
    p_value   <- coef_info["Pr(>|t|)"]
    logFC     <- coef_info["Estimate"]

    message(paste("Gene:", gene,
                  "| logFC:", round(logFC, 3),
                  "| p:", signif(p_value, 3)))

    # Plot
    p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(width = 0.6, size = 1, color = "black", outlier.shape = NA) +
      geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 3) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(y = bquote("Expression of " ~ italic(.(gene)))) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )

    ggsave(
      filename = file.path(BASE_PLOT,
                           paste0(gene, "_", suffix, ".pdf")),
      plot = p,
      width = 3.5,
      height = 3
    )
  }
}

############################################################
# 1. Total Sleep Time
############################################################

plot_deg_boxplot(
  meta_file  = "Total_sleep_time_Metadata.tsv",
  count_file = "Total_sleep_time_Count.tsv",
  group_column = "Type",
  genes = c("ATF4", "ADCY1", "PER1", "MYBBP1A"),
  coef_name = "Groupshort sleep",
  colors = c("Short Sleep" = "#939ca3",
             "Long Sleep"  = "#677bab"),
  suffix = "tst"
)

############################################################
# 2. NREM Ratio
############################################################

plot_deg_boxplot(
  meta_file  = "NREM_ratio_Metadata.tsv",
  count_file = "NREM_ratio_Count.tsv",
  group_column = "Type",
  genes = c("CPT1A", "GHRL", "NPAS2", "ADCY1",
            "PER1", "MYBBP1A", "NR1H3"),
  coef_name = "GroupLow NREM ratio",
  colors = c("Low NREM ratio"  = "#d69e49",
             "High NREM ratio" = "#838469"),
  suffix = "nrem"
)

############################################################
# 3. Sleep Timing
############################################################

plot_deg_boxplot(
  meta_file  = "Sleep_timing_Metadata.tsv",
  count_file = "Sleep_timing_Count.tsv",
  group_column = "Sleep_Timing",
  genes = c("MAPK9", "PER1", "ADCY1", "MYBBP1A", "NR1H3"),
  coef_name = "GroupSleep at night",
  colors = c("Sleep at daytime" = "#c5e8e1",
             "Sleep at night"   = "#3eb59d"),
  suffix = "timing"
)

############################################################
# 4. Sleep Frequency
############################################################

plot_deg_boxplot(
  meta_file  = "Sleep_frequency_Metadata.tsv",
  count_file = "Sleep_frequency_Count.tsv",
  group_column = "Sleep_Times",
  genes = c("PER1", "ADCY1"),
  coef_name = "GroupOnce",
  colors = c("Once"          = "#b8cdab",
             "More than two" = "#476066"),
  suffix = "frequency"
)
