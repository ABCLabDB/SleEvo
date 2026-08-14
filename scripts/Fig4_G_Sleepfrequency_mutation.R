############################################################
# Fig4_G_Sleepfrequency_mutation.R
#
# Figure 4G. Sleep-frequency SNP alignment panels
#   Nucleotide window (±4) with trait dots (Once / More than twice).
#   Targets only:
#     - ATF5   Chr19:49932628  (CDS_49932628)
#     - CARTPT Chr5:71719366   (CDS_71719366)
#
# Input:
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result1/circadian_Gene_list.tsv
#   - data/Result5/SNP.tsv
#   - data/Circadian_gene_Nucleotide_Matrix/{ATF5,CARTPT}.tsv
#
# Output:
#   - Result/Fig4/Fig4_G_ATF5_Chr19_49932628.jpg|.pdf
#   - Result/Fig4/Fig4_G_CARTPT_Chr5_71719366.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

META_FILE <- file.path(PROJECT_DIR, "data", "Result1",
                       "species_sleep_metadata.txt")
CHR_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                       "circadian_Gene_list.tsv")
SNP_FILE  <- file.path(PROJECT_DIR, "data", "Result5", "SNP.tsv")
SEQ_DIR   <- file.path(PROJECT_DIR, "data",
                       "Circadian_gene_Nucleotide_Matrix")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

window_half <- 4


## =========================================================
## 2. Target SNPs (Fig4 G only)
## =========================================================

targets <- data.frame(
  Gene         = c("ATF5", "CARTPT"),
  Chr          = c("19", "5"),
  cds_coord    = c("49932628", "71719366"),
  stringsAsFactors = FALSE
)
targets$target_cds   <- paste0("CDS_", targets$cds_coord)
targets$target_label <- paste0("Chr", targets$Chr, ":", targets$cds_coord)
targets$out_stem     <- paste0(
  "Fig4_G_", targets$Gene, "_Chr", targets$Chr, "_", targets$cds_coord
)


## =========================================================
## 3. Colors
## =========================================================

base_colors <- c(
  "A" = "#476066",
  "T" = "#838469",
  "C" = "#d69e49",
  "G" = "#eadaa0",
  "N" = "#f9f6f2",
  "-" = "#ab5852"
)

group_colors <- c(
  "More than twice" = "#476066",
  "Once"            = "#b8cdab"
)


## =========================================================
## 4. Metadata (sleep COUNT)
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
## 5. Confirm SNP records + chromosome table
## =========================================================

snp_tab <- fread(SNP_FILE) |> as.data.frame()
chr_tab <- fread(CHR_FILE) |> as.data.frame()

for (k in seq_len(nrow(targets))) {
  hit <- snp_tab[
    snp_tab$Gene == targets$Gene[k] &
      snp_tab$Nucleotide_position == targets$target_cds[k], ]
  if (nrow(hit) == 0) {
    stop("SNP not found in SNP.tsv: ", targets$Gene[k], " ",
         targets$target_cds[k])
  }
  message(targets$Gene[k], " ", targets$target_label[k],
          " | adjP = ", signif(hit$adjP[1], 4))

  ## prefer chromosome from gene list if available
  chr_hit <- chr_tab$Chr[match(targets$Gene[k], chr_tab$Gene_symbol)]
  if (!is.na(chr_hit) && nzchar(as.character(chr_hit))) {
    targets$Chr[k] <- as.character(chr_hit)
    targets$target_label[k] <- paste0(
      "Chr", targets$Chr[k], ":", targets$cds_coord[k]
    )
    targets$out_stem[k] <- paste0(
      "Fig4_G_", targets$Gene[k], "_Chr", targets$Chr[k], "_",
      targets$cds_coord[k]
    )
  }
}


## =========================================================
## 6. Plot helper
## =========================================================

plot_freq_snp_alignment <- function(gene, target_cds, target_label) {

  seq_file <- file.path(SEQ_DIR, paste0(gene, ".tsv"))
  if (!file.exists(seq_file)) stop("Missing matrix: ", seq_file)

  aln_df <- fread(seq_file) |> as.data.frame()
  rownames(aln_df) <- tolower(aln_df[, 1])
  aln_df <- aln_df[, -1, drop = FALSE]
  aln_df <- aln_df[, grepl("V|CDS", colnames(aln_df)), drop = FALSE]

  if (!target_cds %in% colnames(aln_df)) {
    stop("Column not found in ", gene, " matrix: ", target_cds)
  }

  aln_df$Number_of_sleep_times_per_day <- meta$Number_of_sleep_times_per_day[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]
  aln_df$Sleep_Type <- meta$Type[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]
  aln_df$Species <- meta$Species_name_ensembl[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]

  aln_df <- aln_df |>
    filter(!is.na(Number_of_sleep_times_per_day),
           !is.na(Sleep_Type),
           Sleep_Type %in% names(group_colors))

  all_pos <- setdiff(
    colnames(aln_df),
    c("Species", "Number_of_sleep_times_per_day", "Sleep_Type")
  )
  j <- which(all_pos == target_cds)
  visible_pos <- all_pos[
    pmax(1, j - window_half):pmin(length(all_pos), j + window_half)
  ]

  ## display names: Chr*:pos for CDS columns
  rename_pos <- function(x) {
    ifelse(
      grepl("^CDS_", x),
      paste0("Chr", chr_tab$Chr[match(gene, chr_tab$Gene_symbol)],
             ":", sub("^CDS_", "", x)),
      x
    )
  }
  visible_lab <- rename_pos(visible_pos)
  names(visible_lab) <- visible_pos
  ## keep focal label explicit
  visible_lab[target_cds] <- target_label

  plot_df <- aln_df |>
    dplyr::select(Species, all_of(visible_pos),
                  Number_of_sleep_times_per_day, Sleep_Type) |>
    tidyr::pivot_longer(cols = all_of(visible_pos),
                        names_to = "Position",
                        values_to = "Base")

  plot_df$Position <- factor(plot_df$Position, levels = visible_pos)
  plot_df <- plot_df |>
    mutate(
      FillColor = ifelse(Position == target_cds, base_colors[Base], "white"),
      FontColor = ifelse(Position == target_cds, "white", "black")
    )

  plot_df_center <- plot_df |>
    filter(Position == target_cds) |>
    dplyr::select(Species, CenterBase = Base)

  species_order <- plot_df |>
    left_join(plot_df_center, by = "Species") |>
    distinct(Species, Sleep_Type, Number_of_sleep_times_per_day, CenterBase) |>
    arrange(Sleep_Type, CenterBase, Number_of_sleep_times_per_day) |>
    pull(Species)

  plot_df$Species <- factor(plot_df$Species, levels = species_order)

  sleep_tile <- aln_df |>
    dplyr::select(Species, Sleep_Type) |>
    distinct() |>
    mutate(Species = factor(Species, levels = species_order))

  n_sp <- length(species_order)

  ## single ggplot (avoids patchwork title/panel gap)
  trait_df <- sleep_tile |>
    mutate(
      x = 0,
      dot_col = unname(group_colors[as.character(Sleep_Type)])
    )

  ## numeric x for positions so trait dots sit just left of the window
  pos_num <- setNames(seq_along(visible_pos), visible_pos)
  plot_df$x <- pos_num[as.character(plot_df$Position)]

  p <- ggplot() +
    geom_point(data = trait_df,
               aes(x = x, y = Species),
               colour = trait_df$dot_col,
               size = 3.8) +
    geom_tile(data = plot_df,
              aes(x = x, y = Species, fill = FillColor),
              width = 0.55, height = 0.92, color = NA) +
    scale_fill_identity() +
    geom_text(data = plot_df,
              aes(x = x, y = Species, label = Base),
              colour = plot_df$FontColor,
              size = 3.8, fontface = "bold") +
    scale_x_continuous(
      breaks = unname(pos_num),
      labels = unname(visible_lab[names(pos_num)]),
      position = "top",
      limits = c(-0.7, length(visible_pos) + 0.4),
      expand = c(0, 0)
    ) +
    scale_y_discrete(limits = species_order, expand = expansion(add = 0.2)) +
    labs(title = paste(gene, "Gene", target_label, "Region"),
         x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5,
                                colour = "#476066",
                                margin = margin(t = 2, b = 6)),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0,
                                     face = "bold", size = 9),
      axis.text.y = element_text(hjust = 0, size = 9, colour = "grey20"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 8, r = 14, b = 4, l = 4)
    )

  attr(p, "n_species") <- n_sp
  p
}


## =========================================================
## 7. Build + save each target
## =========================================================

for (k in seq_len(nrow(targets))) {
  p <- plot_freq_snp_alignment(
    gene         = targets$Gene[k],
    target_cds   = targets$target_cds[k],
    target_label = targets$target_label[k]
  )

  n_sp <- attr(p, "n_species")
  if (is.null(n_sp)) n_sp <- 30
  ## compact row packing: ~0.28\" per species + title/axis chrome
  fig_h <- max(8, 0.28 * n_sp + 1.4)

  stem <- targets$out_stem[k]
  ggsave(file.path(OUT_DIR, paste0(stem, ".jpg")),
         p, width = 8, height = fig_h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")),
         p, width = 8, height = fig_h, bg = "white")

  message("Saved: Result/Fig4/", stem, ".jpg|.pdf  (height=",
          round(fig_h, 2), " in)")
}
