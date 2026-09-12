############################################################
# Fig4_A_Sleeptiming_enrichment.R
#
# Figure 4A. Fisher enrichment lollipop plots
#   (top)    Sleep timing related
#   (bottom) Sleep frequency related
#
# Input:
#   - data/Result4/Sleeptiming_Fisher_Result.tsv
#   - data/Result5/Sleep_frequency_Fisher.tsv
#
# Output:
#   - Result/Fig4/Fig4_A_Sleeptiming_enrichment.jpg|.pdf
#   - Result/Fig4/Fig4_A_Sleepfrequency_enrichment.jpg|.pdf
#   - Result/Fig4/Fig4_A_enrichment_combined.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

TIMING_FILE <- file.path(
  PROJECT_DIR, "data", "Result4",
  "Sleeptiming_Fisher_Result.tsv"
)
FREQ_FILE <- file.path(
  PROJECT_DIR, "data", "Result5",
  "Sleep_frequency_Fisher.tsv"
)

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

threshold <- -log10(0.05)


## =========================================================
## 2. Shared lollipop helper
## =========================================================

make_enrichment_lollipop <- function(df, sig_label, sig_color,
                                    n_sig, n_total, top_n = 50) {
  plot_df <- df |>
    arrange(desc(logP)) |>
    slice_head(n = top_n)

  plot_df$Significant <- ifelse(
    plot_df$logP > threshold, sig_label, "No association"
  )
  plot_df$Significant <- factor(
    plot_df$Significant,
    levels = c(sig_label, "No association")
  )

  fill_cols <- c(
    setNames(sig_color, sig_label),
    "No association" = "#d0d4d8"
  )

  ggplot(plot_df, aes(x = reorder(Gene, -logP), y = logP,
                      color = Significant)) +
    geom_segment(aes(x = reorder(Gene, -logP),
                     xend = reorder(Gene, -logP),
                     y = 0, yend = logP),
                 linewidth = 6) +
    geom_point(size = 6) +
    geom_hline(yintercept = threshold, linetype = "dashed",
               colour = "#c0392b", linewidth = 0.7) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0(n_sig, "/", n_total, " circadian genes"),
             hjust = 1.05, vjust = 1.6,
             size = 3.8, colour = "grey25") +
    scale_color_manual(values = fill_cols, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = expression(-log[10](italic(P)))) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 11),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = 10, face = "italic",
                                 colour = "grey20"),
      axis.text.y = element_text(size = 11, colour = "grey25"),
      axis.title.y = element_text(size = 12, colour = "grey15"),
      axis.line = element_line(colour = "grey30", linewidth = 0.5),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.4),
      plot.margin = margin(8, 12, 6, 8)
    )
}


## =========================================================
## 3. Sleep timing (Fig4A top)
## =========================================================

timing <- fread(TIMING_FILE) |> as.data.frame()
timing_df <- timing |>
  transmute(Gene = Gene, logP = -log10(P_Value))

n_timing_total <- nrow(timing_df)
n_timing_sig   <- sum(timing_df$logP > threshold, na.rm = TRUE)

p_timing <- make_enrichment_lollipop(
  df = timing_df,
  sig_label = "Sleep timing related",
  sig_color = "#008585",
  n_sig = n_timing_sig,
  n_total = n_timing_total,
  top_n = 50
)


## =========================================================
## 4. Sleep frequency (Fig4A bottom; from Result5 Fig5A)
## =========================================================

freq <- fread(FREQ_FILE) |> as.data.frame()
freq_df <- freq |>
  transmute(Gene = Gene, logP = -log10(P_Value))

n_freq_total <- nrow(freq_df)
n_freq_sig   <- sum(freq_df$logP > threshold, na.rm = TRUE)

p_freq <- make_enrichment_lollipop(
  df = freq_df,
  sig_label = "Sleep frequency related",
  sig_color = "#ffa600",
  n_sig = n_freq_sig,
  n_total = n_freq_total,
  top_n = 32
)


## =========================================================
## 5. Save individual + combined panels
## =========================================================

ggsave(file.path(OUT_DIR, "Fig4_A_Sleeptiming_enrichment.jpg"),
       p_timing, width = 12, height = 3.8, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Fig4_A_Sleeptiming_enrichment.pdf"),
       p_timing, width = 12, height = 3.8, bg = "white")

ggsave(file.path(OUT_DIR, "Fig4_A_Sleepfrequency_enrichment.jpg"),
       p_freq, width = 12, height = 3.8, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Fig4_A_Sleepfrequency_enrichment.pdf"),
       p_freq, width = 12, height = 3.8, bg = "white")

p_combined <- p_timing / p_freq +
  plot_layout(heights = c(1, 1))

ggsave(file.path(OUT_DIR, "Fig4_A_enrichment_combined.jpg"),
       p_combined, width = 12, height = 7.2, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Fig4_A_enrichment_combined.pdf"),
       p_combined, width = 12, height = 7.2, bg = "white")

message("Sleep timing: ", n_timing_sig, " / ", n_timing_total)
message("Sleep frequency: ", n_freq_sig, " / ", n_freq_total)
message("Saved: Result/Fig4/Fig4_A_Sleeptiming_enrichment.jpg|.pdf")
message("Saved: Result/Fig4/Fig4_A_Sleepfrequency_enrichment.jpg|.pdf")
message("Saved: Result/Fig4/Fig4_A_enrichment_combined.jpg|.pdf")
