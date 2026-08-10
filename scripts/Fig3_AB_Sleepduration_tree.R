############################################################
# Fig3_AB_Sleepduration_tree.R
#
# Sleep-phylogeny figure in the "tree of life" reference design.
# Rectangular cladogram + right-side dotted leaders + sleep-time bars.
# Uses ggtree + ggtreeExtra.
#
# Input:
#   - data/Fasta/<Gene>_muscle.fasta
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result2/Total_sleep_time_Anova_Result.tsv
#
# Output:
#   - Result/Fig3/Fig3_A_<Gene>.jpg
#   - Result/Fig3/Fig3_A_<Gene>.pdf
#
# Required working directory: repository root (SleEvo / Sleep_Evolution)
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
## 2. Project paths (relative to repository root)
## =========================================================

PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
ANOVA_FILE <- file.path(PROJECT_DIR, "data", "Result2",
                        "Total_sleep_time_Anova_Result.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Sleep-type colors (also drives the legend)
## =========================================================

fill_colors <- c(
  "Long Sleep"  = "#1864ab",
  "Others"      = "#B8BCC2",
  "Short Sleep" = "#c92a2a"
)


## =========================================================
## 4. Load metadata and define sleep categories
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

meta <- meta |>
  filter(!is.na(Total_sleep_time_per_day)) |>
  arrange(Total_sleep_time_per_day)

prop_rate <- 0.33

bottom_species <- meta$Species_name_ensembl[
  1:ceiling(prop_rate * nrow(meta))
]
top_species <- meta$Species_name_ensembl[
  (nrow(meta) - ceiling(prop_rate * nrow(meta)) + 1):nrow(meta)
]

meta$Type <- "Others"
meta$Type[meta$Species_name_ensembl %in% bottom_species] <- "Short Sleep"
meta$Type[meta$Species_name_ensembl %in% top_species]    <- "Long Sleep"


## =========================================================
## 5. Load ANOVA results
## =========================================================

anova_Result <- fread(ANOVA_FILE) |> as.data.frame()


## ---------------------------------------------------------------------------
##  Design helper
## ---------------------------------------------------------------------------
plot_sleep_tree <- function(tree, use_meta, gene_name, cladogram = TRUE) {

  dat <- use_meta %>%
    transmute(label = Species_name_ensembl,
              Type  = factor(Type, levels = names(fill_colors)),
              Sleep = Total_sleep_time_per_day)

  p <- ggtree(tree,
              layout        = "rectangular",
              branch.length = if (cladogram) "none" else "branch.length",
              ladderize     = TRUE,
              linewidth     = 0.6,
              colour        = "#3A4A63") %<+% dat

  ## ---- knobs -------------------------------------------------------------
  ##  offset        : width of the lane before the bars.
  ##                  BIGGER  -> tree takes a smaller share (branches look shorter)
  ##  pwidth        : bar length
  ##  label_size    : species-name font size
  ##  label_colour  : single colour for ALL species names
  ##  char_w_factor : data-units per character, used to start the dotted line
  ##                  just to the RIGHT of each name. If the dots start too
  ##                  early (overlap text) raise it; if too late lower it.
  offset        <- 1.00
  pwidth        <- 0.35
  label_size    <- 4.3
  label_colour  <- "#1f2937"
  char_w_factor <- 0.022

  xmax      <- max(p$data$x, na.rm = TRUE)     # tree width in data units
  bar_start <- xmax * (1 + offset)             # where the bars begin
  tips_df   <- subset(p$data, isTip)

  ## right edge of each label (approx), then start the dots a hair after it
  tips_df$lab_end <- pmin(
    xmax * 1.02 + nchar(as.character(tips_df$label)) * xmax * char_w_factor
    + xmax * 0.01,
    bar_start - xmax * 0.01)

  p +
    ## species names — all one colour
    geom_tiplab(colour = label_colour, size = label_size, offset = xmax * 0.02) +
    ## dotted leader: starts at the RIGHT of each name and runs to the bar
    geom_segment(data = tips_df,
                 aes(x = lab_end, xend = bar_start, y = y, yend = y),
                 linetype = "dotted", linewidth = 0.4, colour = "#c3c7cd",
                 inherit.aes = FALSE) +

    ## ---- sleep-time bars (colored by type; carries the legend) ----
    geom_fruit(geom        = geom_col,
               mapping     = aes(x = Sleep, fill = Type),
               orientation = "y",
               width       = 0.7,
               offset      = offset,
               pwidth      = pwidth,
               axis.params = list(axis       = "x",
                                  text.size  = 3.0,
                                  title      = "Total sleep (h/day)",
                                  title.size = 3.4,
                                  nbreak     = 6,
                                  vline      = FALSE)) +
    scale_fill_manual(values = fill_colors, name = "Sleep type", drop = FALSE) +
    ggtree::hexpand(0.03) +

    theme(plot.title      = element_text(size = 15, face = "bold",
                                         colour = "#194a7a"),
          legend.position = c(0.96, 0.86),
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 11),
          plot.margin     = margin(6, 50, 6, 6)) +
    ggtitle(gene_name)
}


## ===========================================================================
##  Loop
## ===========================================================================

## Only draw these genes; set to NULL to plot all.
target_genes <- c("ADRB1", "ATF4")

for (i in seq_along(anova_Result$Gene)) {
  gene <- anova_Result$Gene[i]

  if (!is.null(target_genes) && !(gene %in% target_genes)) next
  if (anova_Result$Cluster[i] == 2) next

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

  p <- plot_sleep_tree(tree, use_meta, gene, cladogram = TRUE)

  ggsave(filename = file.path(OUT_DIR, paste0("Fig3_A_", gene, ".jpg")),
         plot = p, width = 11, height = 12, dpi = 300, bg = "white",
         limitsize = FALSE)
  ggsave(filename = file.path(OUT_DIR, paste0("Fig3_A_", gene, ".pdf")),
         plot = p, width = 11, height = 12, limitsize = FALSE)

  message("Saved: Fig3_A_", gene)
}

## ---------------------------------------------------------------------------
##  Notes
##  - Tree still too wide?  raise `offset` (1.2, 1.4, ...) — it directly
##    shrinks the tree's share of the width.
##  - Dots not starting exactly at the name's right edge? tweak `char_w_factor`
##    (bigger = dots start further right).
##  - Legend keys are colored squares from the bar fill (Long / Others / Short).
## ---------------------------------------------------------------------------
