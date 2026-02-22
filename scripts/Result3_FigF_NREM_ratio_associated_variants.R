############################################################
# Result3_FigF_NREM_ratio_associated_variants.R
#
# SNP-centered nucleotide profile and phenotype association
# for NREM ratio–associated circadian genes.
#
# Input:
#   data/Metadata/species_sleep_metadata.txt
#   data/Result3/NREM_Ratio_SNP.tsv
#   data/Result3/NREM_key_12_optimal_Kruskal.tsv
#   data/NucleotideMatrix/*.tsv
#
# Output:
#   figures/Result3/SNP_Profile/
#
# Project root required:
#   Sleep_Evolution/
############################################################

## =========================================================
## 1. Libraries
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggseqlogo)
  library(patchwork)
})

## =========================================================
## 2. Paths
## =========================================================

PROJECT_DIR <- getwd()

META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")

SNP_FILE   <- file.path(PROJECT_DIR, "data", "Result3",
                        "NREM_Ratio_SNP.tsv")

PERM_FILE  <- file.path(PROJECT_DIR, "data", "Result3",
                        "NREM_ratio_Anova_Result.tsv")

SEQ_DIR    <- file.path(PROJECT_DIR, "data",
                        "Circadian_gene_Nucleotide_Matrix")

OUT_DIR    <- file.path(PROJECT_DIR,
                        "figures", "Result3",
                        "SNP_Profile")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

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

meta$Sleep_Type <- "Others"
meta$Sleep_Type[meta$Species_name_ensembl %in% top_species]    <- "High Ratio"
meta$Sleep_Type[meta$Species_name_ensembl %in% bottom_species] <- "Low Ratio"

meta <- meta |> filter(Sleep_Type != "Others")

## =========================================================
## 4. Significant genes
## =========================================================

perm_df <- fread(PERM_FILE)
sig_genes <- perm_df |> filter(adj.P < 0.05) |> pull(Gene)

## =========================================================
## 5. SNP data
## =========================================================

snp_df <- fread(SNP_FILE) |>
  filter(p.adj < 0.05,
         Gene %in% sig_genes)

## =========================================================
## 6. Color schemes
## =========================================================

base_colors <- c(
  "A"="#476066",
  "T"="#838469",
  "C"="#d69e49",
  "G"="#eadaa0",
  "N"="#f9f6f2",
  "-"="#ab5852"
)

custom_scheme <- make_col_scheme(
  chars = names(base_colors),
  cols  = base_colors
)

group_colors <- c(
  "High Ratio"="#58508d",
  "Low Ratio" ="#ffa600"
)

## =========================================================
## 7. Main loop
## =========================================================

for (gene in unique(snp_df$Gene)) {

  seq_file <- file.path(SEQ_DIR, paste0(gene, ".tsv"))
  if (!file.exists(seq_file)) next

  aln_df <- fread(seq_file) |> as.data.frame()
  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]

  aln_df <- aln_df[, grepl("V|CDS", colnames(aln_df))]

  aln_df$Species <- rownames(aln_df)

  aln_df <- aln_df |>
    left_join(meta,
              by = c("Species"="Species_name_ensembl")) |>
    filter(!is.na(Sleep_Type))

  snp_positions <- snp_df |>
    filter(Gene == gene) |>
    pull(Real_Position)

  for (pos in snp_positions) {

    if (pos > ncol(aln_df)) next

    visible_pos <- seq(max(1,pos-4),
                       min(ncol(aln_df)-3,pos+4))

    plot_df <- aln_df |>
      select(Species,
             Sleep_Type,
             Percentage_of_NREM_time_per_day,
             all_of(colnames(aln_df)[visible_pos])) |>
      pivot_longer(
        cols = 4:ncol(.),
        names_to = "Position",
        values_to = "Base"
      )

    ## ===============================
    ## SNP tile plot
    ## ===============================

    p_tile <- ggplot(plot_df,
                     aes(Position, Species)) +
      geom_tile(aes(fill=Base),
                color=NA) +
      scale_fill_manual(values=base_colors) +
      theme_minimal() +
      theme(
        axis.text.x=element_text(angle=45,hjust=0),
        panel.grid=element_blank()
      )

    ## ===============================
    ## Sequence logo
    ## ===============================

    center_bases <- aln_df[,pos]

    p_logo <- ggseqlogo(
      center_bases,
      method="probability",
      col_scheme=custom_scheme
    ) +
      theme_minimal()

    ## ===============================
    ## Boxplot
    ## ===============================

    box_df <- data.frame(
      Base = center_bases,
      NREM = aln_df$Percentage_of_NREM_time_per_day
    )

    p_box <- ggplot(box_df,
                    aes(Base,NREM,fill=Base)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.15,alpha=0.6) +
      scale_fill_manual(values=base_colors) +
      theme_minimal()

    final_plot <- (p_tile / (p_logo | p_box)) +
      plot_annotation(
        title=paste(gene,"Position",pos)
      )

    ggsave(file.path(OUT_DIR,
                     paste0(gene,"_pos",pos,".pdf")),
           final_plot,
           width=10,height=12)
  }
}
