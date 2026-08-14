############################################################
# Fig4_F_Sleepfrequency_tree.R
#
# Figure 4F. Sleep-frequency (sleep COUNT) phylogeny
#   Rectangular cladogram + dotted leaders + categorical
#   heatmap column (geom_rect; cell width/gap controlled).
#   Genes: ATF5, NFIL3 (side-by-side combined panel).
#
# Input:
#   - data/Fasta/{ATF5,NFIL3}_muscle.fasta
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result5/Sleep_frequency_Cochran.tsv
#
# Output:
#   - Result/Fig4/Fig4_F_ATF5.jpg|.pdf
#   - Result/Fig4/Fig4_F_NFIL3.jpg|.pdf
#   - Result/Fig4/Fig4_F_ATF5_NFIL3_combined.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


## ---- one-time install (Bioconductor) --------------------------------------
# BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))
# install.packages(c("ggplot2", "dplyr", "data.table", "ape", "patchwork"))


library(data.table)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(ape)
library(Biostrings)
library(patchwork)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
COCHRAN_FILE <- file.path(PROJECT_DIR, "data", "Result5",
                          "Sleep_frequency_Cochran.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Colors (Number of sleep times per day)
## =========================================================

fill_colors <- c(
  "More than twice" = "#476066",
  "Once"            = "#b8cdab"
)
branch_colour <- "#595959"
na_colour     <- "#9aa0a6"


## =========================================================
## 3. Metadata (sleep COUNT)
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()
meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

meta$Number_of_sleep_times_per_day[
  meta$Number_of_sleep_times_per_day == "once"] <- "Once"
meta <- meta[!is.na(meta$Number_of_sleep_times_per_day), ]
rownames(meta) <- NULL

meta$Type <- "More than twice"
meta$Type[meta$Number_of_sleep_times_per_day == "Once"] <- "Once"


## =========================================================
## 4. Gene list (Cochran / knee)
## =========================================================

knee <- fread(COCHRAN_FILE) |> as.data.frame()
## align with Fig4_B / original naming
knee$gene.idx <- knee$Gene
knee$gene.cluster.idx <- as.numeric(knee$Cluster)
rownames(knee) <- NULL


## =========================================================
## 5. Design helper
## =========================================================

plot_count_tree <- function(tree, use_meta, gene_name,
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
    scale_fill_manual(values = fill_colors, name = "Number of sleep",
                      drop = FALSE) +
    annotate("text", x = bar_start + tile_w / 2, y = min(tips_df$y) - 1.4,
             label = "Number of sleep", size = 3.2, colour = "grey30") +
    ggtree::hexpand(0.05) +

    theme(plot.title      = element_text(size = 15, face = "bold",
                                         colour = "#476066"),
          legend.position = "right",
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 11),
          plot.margin     = margin(6, 10, 6, 6)) +
    ggtitle(gene_name)
}


## =========================================================
## 6. Build trees for ATF5 / NFIL3
## =========================================================

target_genes <- c("ATF5", "NFIL3")
plot_list <- list()

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

  p <- plot_count_tree(tree, use_meta, gene,
                       cladogram = TRUE, color_branches = FALSE)

  ggsave(file.path(OUT_DIR, paste0("Fig4_F_", gene, ".jpg")),
         p, width = 11, height = 12, dpi = 300, bg = "white",
         limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0("Fig4_F_", gene, ".pdf")),
         p, width = 11, height = 12, limitsize = FALSE)

  plot_list[[gene]] <- p
  message("Saved: Fig4_F_", gene, ".jpg|.pdf")
}


## =========================================================
## 7. Side-by-side combined panel
## =========================================================

missing <- setdiff(target_genes, names(plot_list))
if (length(missing) > 0) {
  stop("Missing plots for combined panel: ", paste(missing, collapse = ", "))
}

p_combined <- (
  (plot_list[["ATF5"]] + theme(legend.position = "none")) |
    plot_list[["NFIL3"]]
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(file.path(OUT_DIR, "Fig4_F_ATF5_NFIL3_combined.jpg"),
       p_combined, width = 18, height = 12, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(OUT_DIR, "Fig4_F_ATF5_NFIL3_combined.pdf"),
       p_combined, width = 18, height = 12, limitsize = FALSE)

message("Saved: Fig4_F_ATF5_NFIL3_combined.jpg|.pdf")

## ---------------------------------------------------------------------------
##  Notes
##  - Categories: "Once" (#b8cdab) vs "More than twice" (#476066).
##  - cell_h : heatmap cell height. <1 leaves a gap between cells.
##  - cell_w_frac : heatmap cell width.  offset : distance from the names.
##  - color_branches = TRUE -> branches colored by clade-majority category.
## ---------------------------------------------------------------------------
