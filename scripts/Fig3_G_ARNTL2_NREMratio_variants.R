############################################################
# Fig3_G_ARNTL2_NREMratio_variants.R
#
# Figure 3G. ARNTL2 SNP Chr12:27401573 (CDS_27401573)
#   Nucleotide window alignment + logo/boxplot for NREM ratio
#   Adapted from Result2 SNP visualization / Result3 SNP panels,
#   using outputs related to Result3_nrem_ratio_variants.R
#
# Input:
#   - data/Result1/species_sleep_metadata.txt
#   - data/Circadian_gene_Nucleotide_Matrix/ARNTL2.tsv
#   - data/Result3/NREM_Ratio_SNP.tsv
#
# Output:
#   - Result/Fig3/Fig3_G_ARNTL2_Chr12_27401573_alignment.jpg|.pdf
#   - Result/Fig3/Fig3_G_ARNTL2_Chr12_27401573_boxplot_logo.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggseqlogo)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

META_FILE <- file.path(PROJECT_DIR, "data", "Result1",
                       "species_sleep_metadata.txt")
SEQ_FILE  <- file.path(PROJECT_DIR, "data",
                       "Circadian_gene_Nucleotide_Matrix",
                       "ARNTL2.tsv")
SNP_FILE  <- file.path(PROJECT_DIR, "data", "Result3",
                       "NREM_Ratio_SNP.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

gene        <- "ARNTL2"
target_cds  <- "CDS_27401573"
target_label <- "Chr12:27401573"
window_half <- 4


## =========================================================
## 2. Colors
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
  "High Ratio" = "#2c6e49",
  "Low Ratio"  = "#f5d7b0"
)


## =========================================================
## 3. Metadata (NREM ratio groups)
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
## 4. Confirm SNP record
## =========================================================

snp_tab <- fread(SNP_FILE) |> as.data.frame()
snp_row <- snp_tab |>
  filter(Gene == gene, Column_Index == target_cds)

if (nrow(snp_row) == 0) {
  stop("SNP not found in NREM_Ratio_SNP.tsv: ", target_cds)
}
message("SNP ", target_label, " | p.adj = ",
        signif(snp_row$p.adj[1], 4),
        " | method = ", snp_row$Method[1])


## =========================================================
## 5. Load ARNTL2 nucleotide matrix
## =========================================================

aln_df <- fread(SEQ_FILE) |> as.data.frame()
rownames(aln_df) <- aln_df[, 1]
aln_df <- aln_df[, -1, drop = FALSE]
aln_df <- aln_df[, grepl("V|CDS", colnames(aln_df)), drop = FALSE]

if (!target_cds %in% colnames(aln_df)) {
  stop("Column not found in ARNTL2 matrix: ", target_cds)
}

aln_df$NREM_Ratio <- meta$Percentage_of_NREM_time_per_day[
  match(rownames(aln_df), meta$Species_symbol_name_ensembl)
]
aln_df$Type <- meta$Type[
  match(rownames(aln_df), meta$Species_symbol_name_ensembl)
]
aln_df$Species <- meta$Species_name_ensembl[
  match(rownames(aln_df), meta$Species_symbol_name_ensembl)
]

aln_df <- aln_df |>
  filter(!is.na(NREM_Ratio),
         Type %in% c("High Ratio", "Low Ratio"))


## =========================================================
## 6. Window around focal SNP
## =========================================================

all_pos <- setdiff(colnames(aln_df),
                   c("Species", "NREM_Ratio", "Type"))
idx <- which(all_pos == target_cds)
window_pos <- all_pos[
  pmax(1, idx - window_half):pmin(length(all_pos), idx + window_half)
]

plot_df <- aln_df |>
  dplyr::select(Species, all_of(window_pos), NREM_Ratio, Type) |>
  tidyr::pivot_longer(cols = all_of(window_pos),
                      names_to = "Position",
                      values_to = "Base")

plot_df$Position <- factor(plot_df$Position, levels = window_pos)
plot_df <- plot_df |>
  mutate(
    FillColor = ifelse(Position == target_cds, base_colors[Base], "white"),
    FontColor = ifelse(Position == target_cds, "white", "black"),
    PosLabel  = ifelse(Position == target_cds, target_label, as.character(Position))
  )

## x-axis labels: show Chr label only at focal SNP
pos_lab <- setNames(
  ifelse(window_pos == target_cds, target_label,
         sub("^CDS_", "", window_pos)),
  window_pos
)

species_order <- plot_df |>
  filter(Position == target_cds) |>
  arrange(Type, Base, NREM_Ratio) |>
  pull(Species)

plot_df$Species <- factor(plot_df$Species,
                          levels = rev(unique(species_order)))


## =========================================================
## 7. Alignment plot
## =========================================================

p_snp <- ggplot(plot_df, aes(x = Position, y = Species)) +
  geom_tile(aes(fill = FillColor), width = 0.3) +
  geom_text(aes(label = Base, color = FontColor),
            size = 4, fontface = "bold") +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_discrete(position = "top", labels = pos_lab,
                   expand = expansion(mult = c(0.04, 0.10))) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 80, r = 24, b = 8, l = 4),
        axis.text.x.top = element_text(size = 12, angle = 45, hjust = 0,
                                       vjust = 0, colour = "grey15",
                                       margin = margin(b = 12)),
        axis.text.y = element_text(size = 10))

p_trait <- ggplot(
  distinct(plot_df, Species, Type),
  aes(x = 1, y = Species, color = Type)
) +
  geom_point(size = 4) +
  scale_color_manual(values = group_colors) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(t = 80, r = 0, b = 8, l = 4))

alignment_plot <- p_trait + p_snp +
  plot_layout(widths = c(0.3, 4)) +
  plot_annotation(
    title = paste0(gene, "  |  ", target_label),
    theme = theme(
      plot.title = element_text(face = "bold", colour = "#2c6e49",
                                size = 14,
                                margin = margin(t = 6, b = 18)),
      plot.margin = margin(t = 16, r = 16, b = 10, l = 10)
    )
  )


## =========================================================
## 8. Logo + boxplot at focal SNP
## =========================================================

box_df <- plot_df |>
  filter(Position == target_cds) |>
  dplyr::select(Species, Base, NREM_Ratio, Type) |>
  filter(!is.na(Base), Base %in% c("A", "T", "G", "C"))

p_box <- ggplot(box_df,
                aes(x = Base, y = NREM_Ratio, fill = Base)) +
  geom_boxplot(width = 0.5, outlier.shape = NA,
               color = "black", linewidth = 0.8) +
  geom_jitter(aes(color = Base), width = 0.15,
              size = 3, alpha = 0.6) +
  scale_fill_manual(values = base_colors) +
  scale_color_manual(values = base_colors) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Nucleotide type",
    y = "NREM ratio",
    title = paste(gene, "at", target_label)
  ) +
  theme(legend.position = "right")

custom_scheme <- make_col_scheme(
  chars = names(base_colors),
  cols  = base_colors
)

seq_mat <- aln_df[, target_cds, drop = FALSE]
p_logo <- ggseqlogo(seq_mat, method = "bits",
                    col_scheme = custom_scheme) +
  theme_void() +
  labs(title = "Sequence logo")

logo_box_plot <- p_logo / p_box +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(
    title = paste(gene, "functional effect at", target_label),
    theme = theme(plot.title = element_text(face = "bold",
                                            colour = "#2c6e49"))
  )


## =========================================================
## 9. Save
## =========================================================

stem <- paste0("Fig3_G_", gene, "_Chr12_27401573")

ggsave(file.path(OUT_DIR, paste0(stem, "_alignment.jpg")),
       alignment_plot, width = 10.5, height = 14, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(stem, "_alignment.pdf")),
       alignment_plot, width = 10.5, height = 14)

ggsave(file.path(OUT_DIR, paste0(stem, "_boxplot_logo.jpg")),
       logo_box_plot, width = 8, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(stem, "_boxplot_logo.pdf")),
       logo_box_plot, width = 8, height = 8)

message("Saved: Result/Fig3/", stem, "_alignment.jpg|.pdf")
message("Saved: Result/Fig3/", stem, "_boxplot_logo.jpg|.pdf")
