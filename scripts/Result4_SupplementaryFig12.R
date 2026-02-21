# ==========================================================
# SNP Heatmap Visualization (Sleep Timing Analysis)
# GitHub-ready version (relative paths, clean structure)
# ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# -----------------------------
# 1. Load Data
# -----------------------------

metadata <- fread("data/Result1/species_sleep_metadata.txt") %>%
  as.data.frame()

metadata$Species_symbol_name_ensembl <-
  gsub(" ", "_", metadata$Species_symbol_name_ensembl)

metadata <- metadata[!is.na(metadata$Sleep_timing_per_day), ]
metadata$Type <- metadata$Sleep_timing_per_day

# SNP dataset (첫번째 그림 경로)
TS_site1 <- fread("data/Result4/Sleeptiming_Manhattan_Dataset.tsv") %>%
  as.data.frame()

TS_site1 <- TS_site1[TS_site1$P < 0.05, ]
TS_site1$Position <- sapply(strsplit(TS_site1$SNP, ":"), `[`, 2)

# Cochran result (두번째 그림 경로)
perm_test <- fread("data/Result4/Sleeptiming_Cochran_Result.tsv") %>%
  as.data.frame()

knee <- perm_test[perm_test$P.cochran < 0.05, ]
gene <- knee$gene.idx

# nucleotide matrix directory (세번째 그림 경로)
NUC_DIR <- "data/Circadian_gene_Nucleotide_Matrix/"
seqMatrix <- list.files(NUC_DIR, pattern = ".tsv", full.names = TRUE)

seqMatrix <- seqMatrix[
  tools::file_path_sans_ext(basename(seqMatrix)) %in% gene
]

# circadian gene chromosome info (네번째 그림 경로)
circadian_Chr <-
  fread("data/Result1/circadian_Gene_list.tsv") %>%
  as.data.frame()

# -----------------------------
# 2. Visualization
# -----------------------------

OUT_DIR <- "figures/SNP_heatmap/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

base_colors <- c(
  A="#476066", T="#838469",
  C="#d69e49", G="#eadaa0",
  N="#f9f6f2", "-"="#ab5852"
)

group_colors <- c(
  "Sleep at night"="#2c6e49",
  "Sleep at anytime"="grey85",
  "Sleep at daytime"="#d68c45"
)

for(i in seq_along(seqMatrix)) {

  gene_name <- tools::file_path_sans_ext(
    basename(seqMatrix[i])
  )

  aln_df <- fread(seqMatrix[i]) %>% as.data.frame()
  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]

  aln_df <- aln_df[, grepl("^CDS_", colnames(aln_df))]

  rownames(aln_df) <- str_to_title(rownames(aln_df))

  aln_df$Sleep_Timing <-
    metadata$Sleep_timing_per_day[
      match(rownames(aln_df),
            metadata$Species_symbol_name_ensembl)
    ]

  aln_df$Sleep_Type <-
    metadata$Type[
      match(rownames(aln_df),
            metadata$Species_symbol_name_ensembl)
    ]

  aln_df$Species <-
    metadata$Species_name_ensembl[
      match(rownames(aln_df),
            metadata$Species_symbol_name_ensembl)
    ]

  sleep_df <- aln_df %>%
    filter(!is.na(Sleep_Timing))

  snp_pos <- TS_site1$BP[
    TS_site1$Gene == gene_name
  ]

  if(length(snp_pos) == 0) next

  for(bp in snp_pos) {

    center_col <- paste0("CDS_", bp)
    if(!(center_col %in% colnames(sleep_df))) next

    center_index <- which(colnames(sleep_df) == center_col)

    idx_window <- pmax(1, center_index-4):
                  pmin(ncol(sleep_df)-3, center_index+4)

    visible_cols <- colnames(sleep_df)[idx_window]

    plot_df <- sleep_df %>%
      select(Species, all_of(visible_cols),
             Sleep_Timing, Sleep_Type) %>%
      pivot_longer(cols = all_of(visible_cols),
                   names_to="Position",
                   values_to="Base")

    plot_df$Position <- factor(plot_df$Position,
                               levels = visible_cols)

    plot_df <- plot_df %>%
      mutate(
        FillColor =
          ifelse(Position == center_col,
                 base_colors[Base], "white"),
        FontColor =
          ifelse(Position == center_col,
                 "white", "black")
      )

    p_snp <- ggplot(plot_df,
                    aes(Position, Species)) +
      geom_tile(aes(fill=FillColor),
                width=0.3) +
      scale_fill_identity() +
      geom_text(aes(label=Base,
                    color=FontColor),
                size=3.5,
                fontface="bold") +
      scale_color_identity() +
      theme_minimal() +
      theme(
        axis.text.x=element_text(angle=45,
                                 face="bold"),
        axis.title=element_blank(),
        panel.grid=element_blank()
      )

    p_trait <- sleep_df %>%
      select(Species, Sleep_Type) %>%
      ggplot(aes(1, Species,
                 color=Sleep_Type)) +
      geom_point(size=4) +
      scale_color_manual(values=group_colors) +
      theme_void() +
      theme(legend.position="none")

    final_plot <- p_trait + p_snp +
      plot_layout(widths=c(0.2,4)) +
      plot_annotation(
        title=paste(gene_name,
                    "CDS", bp)
      )

    ggsave(
      file.path(OUT_DIR,
                paste0(gene_name, "_", bp, ".pdf")),
      final_plot,
      width=8, height=9
    )
  }
}
