############################################################
# Fig5_D_LOEUF.R
#
# Figure 5D. Gene-level intolerance (LOEUF)
#   1) Permutation test: target circadian genes vs random
#      background genes (one-sided: smaller mean LOEUF)
#   2) Panel D-style horizontal bar (same layout as Fig5_C):
#      LOEUF < 0.6 (green) vs rest (grey)
#
# Input:
#   - data/Cross_validation/TableS22.txt
#   - data/Cross_validation/TableS24.txt
#   - data/Cross_validation/gnomAD_LOEUF_result_all.tsv
#
# Output:
#   - Result/Fig5/Fig5_D_LOEUF_permutation.jpg|.pdf
#   - Result/Fig5/Fig5_D_LOEUF.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)

set.seed(1234)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

TABLE_S22 <- file.path(PROJECT_DIR, "data", "Cross_validation",
                       "TableS22.txt")
TABLE_S24 <- file.path(PROJECT_DIR, "data", "Cross_validation",
                       "TableS24.txt")
GNOMAD_FILE <- file.path(PROJECT_DIR, "data", "Cross_validation",
                         "gnomAD_LOEUF_result_all.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Load data
## =========================================================

supple_22    <- fread(TABLE_S22) |> as.data.frame()
supple_24    <- fread(TABLE_S24) |> as.data.frame()
all_genes_df <- fread(GNOMAD_FILE) |> as.data.frame()

## target (sleep / circadian) gene set
useGene <- unique(supple_22$Gene)
n_circadian <- length(useGene)


## =========================================================
## 3. Background pool = all genes EXCEPT the target set
## =========================================================

## support both formats:
##   - unique gene table (gene, LOEUF)
##   - duplicated-row table (take every 2nd / canonical row)
if ("gene" %in% names(all_genes_df)) {
  gene_col <- "gene"
} else if ("Gene" %in% names(all_genes_df)) {
  gene_col <- "Gene"
} else {
  stop("gnomAD LOEUF table must contain a gene column")
}

if (!"LOEUF" %in% names(all_genes_df)) {
  stop("gnomAD LOEUF table must contain a LOEUF column")
}

## if odd/even duplicate pattern is present, keep even rows
if (nrow(all_genes_df) > 2 * length(unique(all_genes_df[[gene_col]]))) {
  all_genes_df <- all_genes_df[seq(2, nrow(all_genes_df), by = 2), ]
}

all_genes_df$LOEUF <- as.numeric(all_genes_df$LOEUF)

idx <- match(useGene, all_genes_df[[gene_col]])
idx <- idx[!is.na(idx)]
all_genes_df_LOEUF <- all_genes_df[-idx, , drop = FALSE]


## =========================================================
## 4. Target gene LOEUF (TableS24)
## =========================================================

gene_sym_col <- if ("Gene symbol" %in% names(supple_24)) {
  "Gene symbol"
} else if ("Gene_symbol" %in% names(supple_24)) {
  "Gene_symbol"
} else {
  stop("TableS24 must contain Gene symbol")
}

supple_LOEUF <- supple_24[match(useGene, supple_24[[gene_sym_col]]), ]
id_col <- if ("Gene ID" %in% names(supple_LOEUF)) "Gene ID" else "Gene_ID"
supple_LOEUF <- supple_LOEUF[!is.na(supple_LOEUF[[id_col]]), ]
rownames(supple_LOEUF) <- NULL
supple_LOEUF$LOEUF <- as.numeric(supple_LOEUF$LOEUF)


## =========================================================
## 5. Permutation test
##    H1: target genes have a SMALLER mean LOEUF (more constrained)
## =========================================================

obs_mean <- mean(supple_LOEUF$LOEUF, na.rm = TRUE)
n_sample <- sum(!is.na(supple_LOEUF$LOEUF))

bg <- all_genes_df_LOEUF$LOEUF
bg <- bg[!is.na(bg)]

n_perm <- 10000
perm_means <- replicate(
  n_perm,
  mean(sample(bg, size = n_sample, replace = FALSE))
)

## one-sided p: how often a random set is as / more constrained
p_perm <- (sum(perm_means <= obs_mean) + 1) / (n_perm + 1)

message(sprintf("Observed mean LOEUF (target, n=%d): %.4f", n_sample, obs_mean))
message(sprintf("Null mean LOEUF (random sampling):  %.4f", mean(perm_means)))
message(sprintf("Permutation P (target < random):    %.5f", p_perm))


## =========================================================
## 6. Null distribution plot
## =========================================================

perm_df <- data.frame(perm_means = perm_means)

p_perm_plot <- ggplot(perm_df, aes(x = perm_means)) +
  geom_histogram(bins = 50, fill = "grey80", colour = "white") +
  geom_vline(xintercept = obs_mean, colour = "#2CA25F",
             linewidth = 1, linetype = "dashed") +
  annotate("text", x = obs_mean, y = Inf, hjust = -0.05, vjust = 1.5,
           colour = "#2CA25F", fontface = "bold",
           label = paste0("Observed = ", round(obs_mean, 3),
                          "\nP = ", signif(p_perm, 3))) +
  labs(x = paste0("Mean LOEUF of ", n_sample, " random genes (null)"),
       y = "Count",
       title = "LOEUF permutation test (10,000 permutations)") +
  theme_classic(base_size = 13)

ggsave(file.path(OUT_DIR, "LOEUF_permutation.jpg"),
       p_perm_plot, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "LOEUF_permutation.pdf"),
       p_perm_plot, width = 6, height = 4, bg = "white")


## =========================================================
## 7. Panel D bar: LOEUF < 0.6 among circadian genes
##    (same geom_rect layout as Fig5_C_CADD.R)
## =========================================================

n_hi <- sum(supple_LOEUF$LOEUF < 0.6, na.rm = TRUE)
total_n <- n_circadian   # 127 circadian genes (Fig5 D)

## p label: match panel style (e.g. 0.0001)
p_lab <- if (p_perm < 0.0001) {
  format(p_perm, scientific = TRUE, digits = 2)
} else {
  format(round(p_perm, 4), scientific = FALSE)
}

col_hi <- "#2CA25F"   # green (LOEUF < 0.6)
col_lo <- "#D9D9D9"   # rest (grey)

bar_y <- 1.1
bar_h <- 0.8

rect_df <- data.frame(
  xmin = c(0, n_hi),
  xmax = c(n_hi, total_n),
  ymin = bar_y,
  ymax = bar_y + bar_h,
  grp  = c("hi", "lo")
)

p_bar <- ggplot() +
  ## two-segment bar
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = grp),
            colour = NA) +
  ## label inside green segment (two lines, Panel D style)
  annotate("text", x = n_hi / 2, y = bar_y + bar_h / 2,
           label = paste0(n_hi, " circadian genes\n(LOEUF < 0.6)"),
           colour = "white", fontface = "bold", size = 4.2,
           lineheight = 0.95) +
  ## top bracket + total count
  annotate("segment", x = 0, xend = total_n, y = 2.25, yend = 2.25,
           colour = "black", linewidth = 0.5) +
  annotate("text", x = total_n / 2, y = 2.42,
           label = paste0(total_n, " circadian genes"),
           fontface = "bold", size = 4.8) +
  ## "Highly constrained" callout (down, then arrow right)
  annotate("segment", x = n_hi * 0.15, xend = n_hi * 0.15,
           y = bar_y, yend = 0.55, colour = "black", linewidth = 0.5) +
  annotate("segment", x = n_hi * 0.15, xend = n_hi * 0.6,
           y = 0.55, yend = 0.55, colour = "black", linewidth = 0.5,
           arrow = arrow(length = unit(0.10, "inches"), type = "closed")) +
  annotate("text", x = n_hi * 0.65, y = 0.55, hjust = 0,
           label = paste0("Highly constrained (P = ", p_lab, ")"),
           fontface = "bold", size = 4.2) +
  scale_fill_manual(values = c(hi = col_hi, lo = col_lo), guide = "none") +
  coord_cartesian(xlim = c(0, total_n), ylim = c(0, 3), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(6, 30, 6, 6))

ggsave(file.path(OUT_DIR, "Fig5_D_LOEUF.pdf"),
       plot = p_bar, width = 6.5, height = 2.4)
ggsave(file.path(OUT_DIR, "Fig5_D_LOEUF.jpg"),
       plot = p_bar, width = 6.5, height = 2.4, dpi = 300, bg = "white")

message("Circadian genes: ", total_n,
        " | LOEUF < 0.6: ", n_hi,
        " | P = ", signif(p_perm, 4))
message("Saved: Result/Fig5/Fig5_D_LOEUF.jpg|.pdf")
message("Saved: Result/Fig5/Fig5_D_LOEUF_permutation.jpg|.pdf")
