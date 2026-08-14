############################################################
# Fig3_E_NREMratio.R
#
# Figure 3E. NREM-ratio phylogeny (ARNTL2)
#   Same design as Fig3_AB_Sleepduration_tree.R:
#   rectangular cladogram + dotted leaders + right-side bars.
#   Colors: High / Others / Low NREM ratio.
#
# Input:
#   - data/Fasta/ARNTL2_muscle.fasta
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result3/NREM_ratio_Anova_Result.tsv
#
# Output:
#   - Result/Fig3/Fig3_E_ARNTL2.jpg
#   - Result/Fig3/Fig3_E_ARNTL2.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


## ---- one-time install (Bioconductor) --------------------------------------
# install.packages("BiocManager")
# BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))
# install.packages(c("ggnewscale", "ggplot2", "dplyr", "data.table", "ape"))


## =========================================================
## 1. Libraries
## =========================================================

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggplot2)
library(dplyr)
library(data.table)
library(ape)
library(Biostrings)


## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
ANOVA_FILE <- file.path(PROJECT_DIR, "data", "Result3",
                        "NREM_ratio_Anova_Result.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. NREM-ratio colors
## =========================================================

fill_colors <- c(
  "High Ratio" = "#2c6e49",
  "Others"     = "grey",
  "Low Ratio"  = "#f5d7b0"
)


## =========================================================
## 4. Load metadata and define NREM categories
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

meta <- meta |>
  filter(!is.na(Percentage_of_NREM_time_per_day)) |>
  arrange(Percentage_of_NREM_time_per_day)

prop_rate <- 0.33

bottom_species <- meta$Species_name_ensembl[
  1:ceiling(prop_rate * nrow(meta))
]
top_species <- meta$Species_name_ensembl[
  (nrow(meta) - ceiling(prop_rate * nrow(meta)) + 1):nrow(meta)
]

meta$Type <- "Others"
meta$Type[meta$Species_name_ensembl %in% bottom_species] <- "Low Ratio"
meta$Type[meta$Species_name_ensembl %in% top_species]    <- "High Ratio"


## =========================================================
## 5. Load ANOVA / gene list
## =========================================================

anova_Result <- fread(ANOVA_FILE) |> as.data.frame()


## ---------------------------------------------------------------------------
##  Design helper (same layout as Fig3_AB)
## ---------------------------------------------------------------------------
plot_nrem_tree <- function(tree, use_meta, gene_name, cladogram = TRUE) {

  dat <- use_meta %>%
    transmute(label = Species_name_ensembl,
              Type  = factor(Type, levels = names(fill_colors)),
              Value = Percentage_of_NREM_time_per_day)

  p <- ggtree(tree,
              layout        = "rectangular",
              branch.length = if (cladogram) "none" else "branch.length",
              ladderize     = TRUE,
              linewidth     = 0.6,
              colour        = "#3A4A63") %<+% dat

  ## ---- knobs (same as Fig3_AB) ------------------------------------------
  offset        <- 1.00
  pwidth        <- 0.35
  label_size    <- 4.3
  label_colour  <- "#1f2937"
  char_w_factor <- 0.022

  xmax      <- max(p$data$x, na.rm = TRUE)
  bar_start <- xmax * (1 + offset)
  tips_df   <- subset(p$data, isTip)

  tips_df$lab_end <- pmin(
    xmax * 1.02 + nchar(as.character(tips_df$label)) * xmax * char_w_factor
    + xmax * 0.01,
    bar_start - xmax * 0.01)

  p +
    geom_tiplab(colour = label_colour, size = label_size, offset = xmax * 0.02) +
    geom_segment(data = tips_df,
                 aes(x = lab_end, xend = bar_start, y = y, yend = y),
                 linetype = "dotted", linewidth = 0.4, colour = "#c3c7cd",
                 inherit.aes = FALSE) +

    geom_fruit(geom        = geom_col,
               mapping     = aes(x = Value, fill = Type),
               orientation = "y",
               width       = 0.7,
               offset      = offset,
               pwidth      = pwidth,
               axis.params = list(axis       = "x",
                                  text.size  = 3.0,
                                  title      = "NREM (% per day)",
                                  title.size = 3.4,
                                  nbreak     = 6,
                                  vline      = FALSE)) +
    scale_fill_manual(values = fill_colors, name = "NREM ratio", drop = FALSE) +
    ggtree::hexpand(0.03) +

    theme(plot.title      = element_text(size = 15, face = "bold",
                                         colour = "#2c6e49"),
          legend.position = c(0.96, 0.86),
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 11),
          plot.margin     = margin(6, 50, 6, 6)) +
    ggtitle(gene_name)
}


## ===========================================================================
##  Loop — ARNTL2 only
## ===========================================================================

target_genes <- c("ARNTL2")

for (i in seq_along(anova_Result$Gene)) {
  gene <- anova_Result$Gene[i]

  if (!is.null(target_genes) && !(gene %in% target_genes)) next

  fasta_path <- file.path(DATA_FASTA, paste0(gene, "_muscle.fasta"))
  if (!file.exists(fasta_path)) {
    message("Skip (FASTA missing): ", gene)
    next
  }

  cds_muscle <- readDNAStringSet(fasta_path)
  names(cds_muscle) <- sapply(strsplit(names(cds_muscle), ":"), `[`, 2)

  keep <- intersect(names(cds_muscle), meta$Species_symbol_name_ensembl)
  cds_muscle <- cds_muscle[names(cds_muscle) %in% keep]
  names(cds_muscle) <- meta$Species_name_ensembl[
    match(names(cds_muscle), meta$Species_symbol_name_ensembl)]

  use_meta <- meta[match(names(cds_muscle), meta$Species_name_ensembl), ]
  rownames(use_meta) <- NULL

  dna  <- as.DNAbin(cds_muscle)
  dm   <- dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE)
  tree <- njs(dm)
  tree$edge.length[tree$edge.length < 0] <- 0

  p <- plot_nrem_tree(tree, use_meta, gene, cladogram = TRUE)

  ggsave(filename = file.path(OUT_DIR, paste0("Fig3_E_", gene, ".jpg")),
         plot = p, width = 11, height = 12, dpi = 300, bg = "white",
         limitsize = FALSE)
  ggsave(filename = file.path(OUT_DIR, paste0("Fig3_E_", gene, ".pdf")),
         plot = p, width = 11, height = 12, limitsize = FALSE)

  message("Saved: Fig3_E_", gene)
}
