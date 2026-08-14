############################################################
# Fig4_B_PER1_sleeptiming_tree.R
#
# Figure 4B. Sleep-timing phylogeny (PER1)
#   Rectangular cladogram + dotted leaders + categorical
#   heatmap column (geom_rect; cell width/gap controlled).
#
# Input:
#   - data/Fasta/PER1_muscle.fasta
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result4/Sleeptiming_Cochran_Result.tsv
#
# Output:
#   - Result/Fig4/Fig4_B_PER1_sleeptiming_tree.jpg
#   - Result/Fig4/Fig4_B_PER1_sleeptiming_tree.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


## ---- one-time install (Bioconductor) --------------------------------------
# BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))
# install.packages(c("ggplot2", "dplyr", "data.table", "ape"))


library(data.table)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(ape)
library(Biostrings)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
COCHRAN_FILE <- file.path(PROJECT_DIR, "data", "Result4",
                          "Sleeptiming_Cochran_Result.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Colors
## =========================================================

fill_colors <- c(
  "Sleep at night"   = "#2c6e49",
  "Sleep at anytime" = "grey85",
  "Sleep at daytime" = "#d68c45"
)
branch_colour <- "#595959"
na_colour     <- "#9aa0a6"


## =========================================================
## 3. Metadata (Sleep Timing)
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()
meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

## match original script naming
meta$Sleep_Timing <- meta$Sleep_timing_per_day
meta$Sleep_Timing[meta$Species_name_ensembl == "Human"] <- "Sleep at night"
meta <- meta[!is.na(meta$Sleep_Timing), ]
rownames(meta) <- NULL

meta$Type <- "Sleep at anytime"
meta$Type[meta$Sleep_Timing == "Sleep at night"]   <- "Sleep at night"
meta$Type[meta$Sleep_Timing == "Sleep at daytime"] <- "Sleep at daytime"


## =========================================================
## 4. Gene list (Cochran / knee)
## =========================================================

knee <- fread(COCHRAN_FILE) |> as.data.frame()
knee$gene.cluster.idx <- as.numeric(
  sapply(strsplit(as.character(knee$gene.cluster.idx), "_"), `[`, 2)
)
rownames(knee) <- NULL


## =========================================================
## 5. Design helper
## =========================================================

plot_timing_tree <- function(tree, use_meta, gene_name,
                             cladogram = TRUE, color_branches = FALSE) {

  dat <- use_meta %>%
    transmute(label = Species_name_ensembl,
              Type  = factor(Type, levels = names(fill_colors)))

  if (color_branches) {
    p <- ggtree(tree, layout = "rectangular",
                branch.length = if (cladogram) "none" else "branch.length",
                ladderize = TRUE, linewidth = 0.9) %<+% dat
    ntip  <- length(tree$tip.label); total <- ntip + tree$Nnode
    types <- names(fill_colors)
    tip_type <- as.character(dat$Type)[match(tree$tip.label, dat$label)]
    counts <- matrix(0, total, length(types), dimnames = list(NULL, types))
    for (j in seq_len(ntip)) if (!is.na(tip_type[j])) counts[j, tip_type[j]] <- 1
    po <- reorder(tree, "postorder")$edge
    for (r in seq_len(nrow(po))) {
      counts[po[r, 1], ] <- counts[po[r, 1], ] + counts[po[r, 2], ]
    }
    maj <- apply(counts, 1, function(z) {
      if (all(z == 0)) NA_character_ else names(z)[which.max(z)]
    })
    p$data$branch_type <- factor(maj[p$data$node], levels = types)
    p <- p + aes(colour = branch_type)
    branch_scale <- scale_colour_manual(
      values = fill_colors, na.value = na_colour, guide = "none"
    )
  } else {
    p <- ggtree(tree, layout = "rectangular",
                branch.length = if (cladogram) "none" else "branch.length",
                ladderize = TRUE, linewidth = 0.6,
                colour = branch_colour) %<+% dat
    branch_scale <- NULL
  }

  ## ---- knobs -------------------------------------------------------------
  offset        <- 0.85
  label_size    <- 4.3
  label_colour  <- "#1f2937"
  char_w_factor <- 0.022
  cell_w_frac   <- 0.06
  cell_h        <- 0.8

  xmax      <- max(p$data$x, na.rm = TRUE)
  bar_start <- xmax * (1 + offset)
  tile_w    <- xmax * cell_w_frac
  tips_df   <- subset(p$data, isTip)
  tips_df$lab_end <- pmin(
    xmax * 1.02 + nchar(as.character(tips_df$label)) * xmax * char_w_factor
    + xmax * 0.01,
    bar_start - xmax * 0.01)

  p +
    branch_scale +
    geom_tiplab(colour = label_colour, size = label_size,
                offset = xmax * 0.02) +
    geom_segment(data = tips_df,
                 aes(x = lab_end, xend = bar_start, y = y, yend = y),
                 linetype = "dotted", linewidth = 0.4, colour = "#c3c7cd",
                 inherit.aes = FALSE) +

    ## categorical HEATMAP (manual rect; cell_h controls the GAP)
    geom_rect(data = tips_df,
              aes(xmin = bar_start, xmax = bar_start + tile_w,
                  ymin = y - cell_h / 2, ymax = y + cell_h / 2,
                  fill = Type),
              inherit.aes = FALSE, colour = NA) +
    scale_fill_manual(values = fill_colors, name = "Sleep timing",
                      drop = FALSE) +
    annotate("text", x = bar_start + tile_w / 2, y = min(tips_df$y) - 1.4,
             label = "Sleep timing", size = 3.2, colour = "grey30") +
    ggtree::hexpand(0.05) +

    theme(plot.title      = element_text(size = 15, face = "bold",
                                         colour = "#2c6e49"),
          legend.position = "right",
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 11),
          plot.margin     = margin(6, 10, 6, 6)) +
    ggtitle(gene_name)
}


## =========================================================
## 6. Loop — PER1 only
## =========================================================

target_genes <- c("PER1")

for (i in seq_along(knee$gene.idx)) {
  gene <- knee$gene.idx[i]

  if (!is.null(target_genes) && !(gene %in% target_genes)) next
  if (knee$gene.cluster.idx[i] == 2) next

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

  p <- plot_timing_tree(tree, use_meta, gene,
                        cladogram = TRUE, color_branches = FALSE)

  ggsave(file.path(OUT_DIR, paste0("Fig4_B_", gene, "_sleeptiming_tree.jpg")),
         p, width = 11, height = 12, dpi = 300, bg = "white",
         limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0("Fig4_B_", gene, "_sleeptiming_tree.pdf")),
         p, width = 11, height = 12, limitsize = FALSE)

  message("Saved: Fig4_B_", gene, "_sleeptiming_tree.jpg|.pdf")
}

## ---------------------------------------------------------------------------
##  Notes
##  - cell_h  : heatmap cell height. <1 leaves a GAP between cells.
##  - cell_w_frac : heatmap cell width.
##  - color_branches = TRUE -> branches colored by clade-majority timing.
## ---------------------------------------------------------------------------
