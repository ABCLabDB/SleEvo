############################################################
# Fig5_H_DEG_boxplot.R
#
# Figure 5H. DEG boxplots (expression ~ Group + Species:Sequencer)
#   Phenotypes: total sleep time, NREM ratio, sleep timing,
#   sleep frequency.
#
# Input (data/Result6/):
#   - Total_sleep_time_{Metadata,Count}.tsv
#   - NREM_ratio_{Metadata,Count}.tsv
#   - Sleep_timing_{Metadata,Count}.tsv
#   - Sleep_frequency_{Metadata,Count}.tsv
#
# Output:
#   - Result/Fig5/Fig5_H_DEG_boxplot/*.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()
BASE_DATA   <- file.path(PROJECT_DIR, "data", "Result6")
OUT_DIR     <- file.path(PROJECT_DIR, "Result", "Fig5",
                         "Fig5_H_DEG_boxplot")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


## =========================================================
## 2. Plot helper
## =========================================================

plot_deg_boxplot <- function(meta_file,
                             count_file,
                             group_column,
                             genes,
                             coef_name,
                             colors,
                             keep_groups = NULL,
                             suffix) {

  meta <- fread(file.path(BASE_DATA, meta_file)) |> as.data.frame()
  if (!"Run" %in% names(meta)) {
    stop("Metadata missing Run column: ", meta_file)
  }
  rownames(meta) <- meta$Run

  count <- fread(file.path(BASE_DATA, count_file)) |> as.data.frame()
  rownames(count) <- count[, 1]
  count <- count[, -1, drop = FALSE]

  for (gene in genes) {
    if (!gene %in% rownames(count)) {
      message("Skip (gene missing in count): ", gene, " [", suffix, "]")
      next
    }

    Expression <- as.numeric(count[gene, meta$Run])

    df <- data.frame(
      Sample     = meta$Run,
      Group      = meta[[group_column]],
      Sequencer  = meta$Instrument,
      Species    = meta$Species,
      Expression = Expression,
      stringsAsFactors = FALSE
    )

    df <- df[!is.na(df$Group) & df$Group != "" & !is.na(df$Expression), ]
    if (!is.null(keep_groups)) {
      df <- df[df$Group %in% keep_groups, ]
    }
    if (nrow(df) < 4 || all(df$Expression == 0)) next

    ## factor levels = colour order; first level is LM reference
    df$Group <- factor(df$Group, levels = names(colors))

    fit <- tryCatch(
      lm(Expression ~ Group + Species:Sequencer, data = df),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    lm_summary <- summary(fit)
    if (!coef_name %in% rownames(lm_summary$coefficients)) {
      message("Skip (coef missing): ", gene, " | wanted ", coef_name,
              " | available: ",
              paste(rownames(lm_summary$coefficients), collapse = ", "))
      next
    }

    coef_info <- lm_summary$coefficients[coef_name, ]
    p_value   <- coef_info["Pr(>|t|)"]
    logFC     <- coef_info["Estimate"]

    message("Gene: ", gene,
            " | logFC: ", round(logFC, 3),
            " | p: ", signif(p_value, 3),
            " | ", suffix)

    p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(width = 0.6, linewidth = 1, color = "black",
                   outlier.shape = NA) +
      geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 3) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(
        y = bquote("Expression of " ~ italic(.(gene))),
        title = paste0(gene, "  (", suffix, ")"),
        subtitle = paste0("logFC = ", round(logFC, 3),
                          "   p = ", signif(p_value, 3))
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.border = element_rect(color = "black", fill = NA,
                                    linewidth = 1)
      )

    stem <- paste0("Fig5_H_", gene, "_", suffix)
    ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")),
           p, width = 3.5, height = 3.4, bg = "white")
    ggsave(file.path(OUT_DIR, paste0(stem, ".jpg")),
           p, width = 3.5, height = 3.4, dpi = 300, bg = "white")
  }
}


## =========================================================
## 3. Total sleep time
## =========================================================

plot_deg_boxplot(
  meta_file    = "Total_sleep_time_Metadata.tsv",
  count_file   = "Total_sleep_time_Count.tsv",
  group_column = "Type",
  genes        = c("ATF4", "ADCY1", "PER1", "MYBBP1A"),
  coef_name    = "Groupshort sleep",
  keep_groups  = c("long sleep", "short sleep"),
  ## first level = reference for lm
  colors = c("long sleep"  = "#677bab",
             "short sleep" = "#939ca3"),
  suffix = "tst"
)


## =========================================================
## 4. NREM ratio
## =========================================================

plot_deg_boxplot(
  meta_file    = "NREM_ratio_Metadata.tsv",
  count_file   = "NREM_ratio_Count.tsv",
  group_column = "Type",
  genes = c("CPT1A", "GHRL", "NPAS2", "ADCY1",
            "PER1", "MYBBP1A", "NR1H3"),
  coef_name   = "Grouplow ratio",
  keep_groups = c("high ratio", "low ratio"),
  colors = c("high ratio" = "#838469",
             "low ratio"  = "#d69e49"),
  suffix = "nrem"
)


## =========================================================
## 5. Sleep timing
## =========================================================

plot_deg_boxplot(
  meta_file    = "Sleep_timing_Metadata.tsv",
  count_file   = "Sleep_timing_Count.tsv",
  group_column = "Sleep_Timing",
  genes = c("MAPK9", "PER1", "ADCY1", "MYBBP1A", "NR1H3"),
  coef_name   = "GroupSleep at night",
  keep_groups = c("Sleep at daytime", "Sleep at night"),
  colors = c("Sleep at daytime" = "#c5e8e1",
             "Sleep at night"   = "#3eb59d"),
  suffix = "timing"
)


## =========================================================
## 6. Sleep frequency
## =========================================================

plot_deg_boxplot(
  meta_file    = "Sleep_frequency_Metadata.tsv",
  count_file   = "Sleep_frequency_Count.tsv",
  group_column = "Sleep_Times",
  genes = c("PER1", "ADCY1"),
  coef_name   = "GroupOnce",
  keep_groups = c("More than two", "Once"),
  colors = c("More than two" = "#476066",
             "Once"          = "#b8cdab"),
  suffix = "frequency"
)

message("Saved DEG boxplots under: Result/Fig5/Fig5_H_DEG_boxplot/")
