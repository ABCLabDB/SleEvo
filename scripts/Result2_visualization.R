############################################################
# Result2_visualization.R
#
# Visualization of evolutionary association between
# circadian genes and total sleep time across species.
#
# This script generates all main figures for Result2:
#
#  Figure 2A. Phylogenetic tree visualization
#    - Gene-wise dendrograms based on CDS similarity
#    - Branch colors indicate majority sleep phenotype
#
#  Figure 2B. Box / strip plot of sleep phenotypes
#    - Comparison of total sleep time between
#      evolutionarily matched Long_sleep vs Short_sleep groups
#
#  Figure 2E. Phylogenetic signal across genes
#    - Blomberg’s K and Moran’s I statistics
#
#  Figure 2F. SNP-centered nucleotide visualization
#    - Local nucleotide variation around sleep-associated SNPs
#
#  Figure 2H. Selection signature heatmap
#
# All paths are relative to the project root.
#
# Required working directory:
#   setwd("Sleep_Evolution")
############################################################


## =========================================================
## 1. Libraries
## =========================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

library(Biostrings)
library(ape)
library(dendextend)

library(ggplot2)
library(ggrepel)
library(ggdendro)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

library(pheatmap)
library(Cairo)
library(grid)
library(gridExtra)

## =========================================================
## 2. Project paths
## =========================================================

PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
STAT_FILE  <- file.path(PROJECT_DIR, "data", "Result2",
                        "TST_key_12_optimal_Kruskal.tsv")
SIGNAL_FILE <- file.path(PROJECT_DIR, "data", "Result2",
                         "Phylogenetic_Signal_Data_permutation.tsv")
SNP_SITE_FILE <- file.path(PROJECT_DIR, "data", "Result2",
                           "Total_Sleep_Time_SNP.tsv")
SEQ_MATRIX_DIR <- file.path(PROJECT_DIR, "data", "Result2",
                            "NucleotideMatrix")
HEATMAP_DATA_FILE <- file.path(PROJECT_DIR, "data", "Heatmap",
                               "Total_sleep_time.tsv")

FIG_DIR_TREE   <- file.path(PROJECT_DIR, "figures", "Result2", "PhyloTree")
FIG_DIR_BOX    <- file.path(PROJECT_DIR, "figures", "Result2", "Boxplot")
FIG_DIR_SIGNAL <- file.path(PROJECT_DIR, "figures", "Result2",
                            "Phylogenetic_Signal")
FIG_DIR_SNP    <- file.path(PROJECT_DIR, "figures", "Result2", "SNP")
FIG_DIR_HEAT   <- file.path(PROJECT_DIR, "figures", "Result2", "Heatmap")

dir.create(FIG_DIR_TREE,   recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_BOX,    recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_SIGNAL, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_SNP,    recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_HEAT, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Load metadata and define sleep categories
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
meta$Type[meta$Species_name_ensembl %in% bottom_species] <- "Short_sleep"
meta$Type[meta$Species_name_ensembl %in% top_species]    <- "Long_sleep"


## =========================================================
## 4. Load Result2 statistics
## =========================================================

stat <- fread(STAT_FILE) |> as.data.frame() |>
  filter(adj.P < 0.05)

stat$Cluster <- as.numeric(stat$Cluster)


## =========================================================
## 5. Load phylogenetic signal & SNP data
## =========================================================

signal_df <- fread(SIGNAL_FILE) |> as.data.frame()
signal_df[is.na(signal_df)] <- 0
colnames(signal_df)[c(2, 6)] <- c("Blomberg’s K", "Moran’s I")

signal_long <- signal_df |>
  dplyr::select(Gene, `Blomberg’s K`, `Moran’s I`) |>
  tidyr::pivot_longer(
    cols = c(`Blomberg’s K`, `Moran’s I`),
    names_to = "Signal",
    values_to = "Value"
  )

ts_sites <- fread(SNP_SITE_FILE) |> as.data.frame()


## =========================================================
## Figure 2A. Phylogenetic tree visualization
## =========================================================

sleep_colors <- c(
  "Long_sleep"  = "#194a7a",
  "Others"      = "grey70",
  "Short_sleep" = "#E7B800"
)

color_tree_by_sleep <- function(dend, sleep_info) {
  
  assign_color <- function(d) {
    
    if (is.leaf(d)) {
      st <- sleep_info[labels(d)]
      attr(d, "Type") <- st
      attr(d, "edgePar") <- list(col = sleep_colors[st], lwd = 3)
      return(d)
    }
    
    d[[1]] <- assign_color(d[[1]])
    d[[2]] <- assign_color(d[[2]])
    
    types <- na.omit(c(attr(d[[1]], "Type"), attr(d[[2]], "Type")))
    maj <- if (length(types) == 0) NA else names(sort(table(types), TRUE))[1]
    
    attr(d, "Type") <- maj
    attr(d, "edgePar") <- list(
      col = ifelse(is.na(maj), "grey70", sleep_colors[maj]),
      lwd = 3
    )
    d
  }
  
  assign_color(dend)
}

for (i in seq_len(nrow(stat))) {
  
  gene <- stat$Gene[i]
  k    <- stat$Cluster[i]
  if (k == 2) next
  
  fasta_path <- file.path(DATA_FASTA, paste0(gene, "_muscle.fasta"))
  if (!file.exists(fasta_path)) next
  
  cds <- readDNAStringSet(fasta_path)
  names(cds) <- sapply(strsplit(names(cds), ":"), `[`, 2)
  
  keep <- intersect(names(cds), meta$Species_symbol_name_ensembl)
  cds  <- cds[names(cds) %in% keep]
  
  names(cds) <- meta$Species_name_ensembl[
    match(names(cds), meta$Species_symbol_name_ensembl)
  ]
  
  use_meta <- meta[match(names(cds), meta$Species_name_ensembl), ]
  
  dna  <- as.DNAbin(cds)
  tree <- njs(dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE))
  dend <- as.dendrogram(hclust(as.dist(cophenetic.phylo(tree)), "average"))
  
  dend_colored <- color_tree_by_sleep(
    dend, setNames(use_meta$Type, use_meta$Species_name_ensembl)
  )
  
  leaf_types  <- use_meta$Type[match(labels(dend_colored),
                                     use_meta$Species_name_ensembl)]
  leaf_colors <- sleep_colors[leaf_types]
  leaf_colors[is.na(leaf_colors)] <- "grey70"
    
  dend_colored <- dend_colored %>%
    set("labels_cex", 0.8) %>%
    set("leaves_pch", 19) %>%        # ● solid circle
    set("leaves_col", leaf_colors) %>%
    set("leaves_cex", 1.5)
  
  pdf(file.path(FIG_DIR_TREE, paste0(gene, "_phylo_tree.pdf")), 12, 6)
  plot(dend_colored, main = gene)
  legend("topright", names(sleep_colors),
         fill = sleep_colors, title = "Sleep Type")
  dev.off()
}


## =========================================================
## Figure 2B. Box / strip plot of sleep phenotypes
## =========================================================

for (i in seq_len(nrow(stat))) {
  
  gene <- stat$Gene[i]
  k    <- stat$Cluster[i]
  if (k == 2) next
  
  fasta_path <- file.path(DATA_FASTA, paste0(gene, "_muscle.fasta"))
  if (!file.exists(fasta_path)) next
  
  cds <- readDNAStringSet(fasta_path)
  names(cds) <- sapply(strsplit(names(cds), ":"), `[`, 2)
  names(cds) <- str_to_title(names(cds))
  
  keep <- intersect(tolower(names(cds)), meta$Species_symbol_name_ensembl)
  cds  <- cds[tolower(names(cds)) %in% keep]
  
  names(cds) <- meta$Species_name_ensembl[
    match(tolower(names(cds)), meta$Species_symbol_name_ensembl)
  ]
  
  use_meta <- meta[match(names(cds), meta$Species_name_ensembl), ]
  
  dna <- as.DNAbin(cds)
  hc <- hclust(as.dist(cophenetic.phylo(
    njs(dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE))
  )), method = "average")
  
  clusters <- cutree(hc, k = k)
  
  class_assign <- rep(NA, length(clusters))
  names(class_assign) <- names(clusters)
  
  class_assign[names(clusters) %in% top_species]    <- "Long_sleep"
  class_assign[names(clusters) %in% bottom_species] <- "Short_sleep"
  
  sleep_df <- use_meta
  sleep_df$Infer_Class <- "Others"
  sleep_df$Infer_Class[sleep_df$Species_name_ensembl %in%
                         names(class_assign[class_assign == "Long_sleep"])] <- "Long_sleep"
  sleep_df$Infer_Class[sleep_df$Species_name_ensembl %in%
                         names(class_assign[class_assign == "Short_sleep"])] <- "Short_sleep"
  
  sleep_df$Real_Class <- "Others"
  sleep_df$Real_Class[sleep_df$Species_name_ensembl %in% top_species]    <- "Long_sleep"
  sleep_df$Real_Class[sleep_df$Species_name_ensembl %in% bottom_species] <- "Short_sleep"
  
  plot_df <- sleep_df |>
    filter(Infer_Class == Real_Class,
           Infer_Class %in% c("Long_sleep", "Short_sleep"))
  
  if (nrow(plot_df) < 4) next
  
  group_means <- plot_df |>
    group_by(Infer_Class) |>
    summarise(mean_sleep = mean(Total_sleep_time_per_day, na.rm = TRUE),
              .groups = "drop")
  
  overall_mean <- mean(plot_df$Total_sleep_time_per_day, na.rm = TRUE)
  
  # jitter 좌표
  plot_df$jitter_x <- as.numeric(factor(plot_df$Infer_Class)) +
    runif(nrow(plot_df), -0.25, 0.25)
  
  p <- ggplot() +
    
    # 박스플롯
    geom_boxplot(
      data = plot_df,
      aes(x = Infer_Class,
          y = Total_sleep_time_per_day,
          fill = Infer_Class),
      color = "black",
      linewidth = 0.9,
      width = 0.5,
      outlier.shape = NA
    ) +
    
    # 점
    geom_point(
      data = plot_df,
      aes(x = jitter_x,
          y = Total_sleep_time_per_day,
          color = Infer_Class),
      size = 5,
      alpha = 0.5
    ) +
    
    # 라벨
    geom_text_repel(
      data = plot_df,
      aes(x = jitter_x,
          y = Total_sleep_time_per_day,
          label = Species_name_ensembl),
      color = "black",
      size = 4,
      max.overlaps = 4,
      box.padding = 0.5,
      point.padding = 0.25,
      segment.color = "gray60",
      seed = 42
    ) +
    
    # Long_sleep mean
    geom_hline(
      yintercept = group_means$mean_sleep[group_means$Infer_Class == "Long_sleep"],
      color = sleep_colors["Long_sleep"],
      linetype = "dashed",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 2.3,
      y = group_means$mean_sleep[group_means$Infer_Class == "Long_sleep"] + 1,
      label = sprintf("Long_sleep mean: %.2f h",
                      group_means$mean_sleep[group_means$Infer_Class == "Long_sleep"]),
      color = sleep_colors["Long_sleep"],
      size = 4,
      fontface = "italic"
    ) +
    
    # Short_sleep mean
    geom_hline(
      yintercept = group_means$mean_sleep[group_means$Infer_Class == "Short_sleep"],
      color = sleep_colors["Short_sleep"],
      linetype = "dashed",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 0.7,
      y = group_means$mean_sleep[group_means$Infer_Class == "Short_sleep"] + 1,
      label = sprintf("Short_sleep mean: %.2f h",
                      group_means$mean_sleep[group_means$Infer_Class == "Short_sleep"]),
      color = sleep_colors["Short_sleep"],
      size = 4,
      fontface = "italic"
    ) +
    
    # 전체 평균
    geom_hline(
      yintercept = overall_mean,
      color = "gray40",
      linetype = "dashed",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 0.7,
      y = overall_mean + 1,
      label = sprintf("Overall mean: %.2f h", overall_mean),
      color = "gray30",
      size = 4,
      fontface = "italic"
    ) +
    
    scale_fill_manual(values = sleep_colors) +
    scale_color_manual(values = sleep_colors) +
    
    labs(
      x = "",
      y = "Sleep duration"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
      plot.title = element_text(hjust = 0.5)
    ) +
    
    coord_cartesian(ylim = c(0, 23.5))
  
  ggsave(file.path(FIG_DIR_BOX, paste0(gene, "_boxplot.pdf")),
         p, width = 10, height = 4)
}


## =========================================================
## Figure 2E. Phylogenetic signal across genes
## =========================================================

p_signal <- ggplot(signal_long,
                   aes(x = Gene, y = Value, color = Signal)) +
  geom_line(aes(group = Signal),
            linewidth = 1.2, alpha = 0.4, color = "grey60") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 1.5, color = "#95be8d") +
  scale_color_manual(values = c(
    "Blomberg’s K" = "#7593af",
    "Moran’s I"    = "#d69e49"
  )) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(y = "Phylogenetic signal score", x = NULL)

ggsave(file.path(FIG_DIR_SIGNAL,
                 "Figure2E_Phylogenetic_Signal.pdf"),
       p_signal, width = 6, height = 3.5)


## =========================================================
## Figure 2F. SNP-centered nucleotide visualization
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
  "Long_sleep"  = "#58508d",
  "Short_sleep" = "#ffa600"
)

for (gene in unique(ts_sites$Gene)) {

  seq_file <- file.path(SEQ_MATRIX_DIR, paste0(gene, ".tsv"))
  if (!file.exists(seq_file)) next

  aln_df <- fread(seq_file) |> as.data.frame()
  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]

  aln_df <- aln_df[, grepl("V|CDS", colnames(aln_df)), drop = FALSE]
  rownames(aln_df) <- str_to_title(rownames(aln_df))

  aln_df$Total_sleep <- meta$Total_sleep_time_per_day[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]
  aln_df$Sleep_Type <- meta$Type[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]
  aln_df$Species <- meta$Species_name_ensembl[
    match(rownames(aln_df), meta$Species_symbol_name_ensembl)
  ]

  aln_df <- aln_df |>
    filter(!is.na(Total_sleep),
           Sleep_Type %in% c("Long_sleep", "Short_sleep"))

  for (pos in ts_sites$Real_Position[ts_sites$Gene == gene]) {

    all_pos <- setdiff(colnames(aln_df),
                       c("Species", "Total_sleep", "Sleep_Type"))
    if (!pos %in% all_pos) next

    idx <- which(all_pos == pos)
    window_pos <- all_pos[pmax(1, idx - 4):pmin(length(all_pos), idx + 4)]

    plot_df <- aln_df |>
      select(Species, all_of(window_pos), Total_sleep, Sleep_Type) |>
      pivot_longer(cols = all_of(window_pos),
                   names_to = "Position",
                   values_to = "Base")

    plot_df$Position <- factor(plot_df$Position, levels = window_pos)

    plot_df <- plot_df |>
      mutate(
        FillColor = ifelse(Position == pos, base_colors[Base], "white"),
        FontColor = ifelse(Position == pos, "white", "black")
      )

    species_order <- plot_df |>
      filter(Position == pos) |>
      arrange(Sleep_Type, Base, Total_sleep) |>
      pull(Species)

    plot_df$Species <- factor(plot_df$Species,
                              levels = rev(unique(species_order)))

    p_snp <- ggplot(plot_df,
                    aes(x = Position, y = Species)) +
      geom_tile(aes(fill = FillColor), width = 0.3) +
      geom_text(aes(label = Base, color = FontColor),
                size = 4, fontface = "bold") +
      scale_fill_identity() +
      scale_color_identity() +
      scale_x_discrete(position = "top") +
      theme_minimal() +
      theme(axis.title = element_blank(),
            panel.grid = element_blank())

    p_trait <- ggplot(
      distinct(plot_df, Species, Sleep_Type),
      aes(x = 1, y = Species, color = Sleep_Type)
    ) +
      geom_point(size = 4) +
      scale_color_manual(values = group_colors) +
      theme_void() +
      theme(legend.position = "none")

    final_plot <- p_trait + p_snp +
      plot_layout(widths = c(0.2, 4)) +
      plot_annotation(
        title = paste(gene, "SNP", pos),
        theme = theme(plot.title = element_text(face = "bold"))
      )

    ggsave(file.path(FIG_DIR_SNP,
                     paste0(gene, "_", pos, "_SNP.pdf")),
           final_plot, width = 8, height = 12, dpi = 300)
  }
}


## =========================================================
## Figure 2H. Selection signature heatmap
## =========================================================

merge_signDF <- fread(HEATMAP_DATA_FILE) |> as.data.frame()

sig_genes <- stat$Gene
merge_signDF <- merge_signDF |>
  filter(Gene %in% sig_genes) |>
  arrange(desc(Cluster), Gene)

rownames(merge_signDF) <- merge_signDF$Gene

mat <- merge_signDF |>
  select(Gene, Ps, π, `Tajima's D`, `dN/dS`)
rownames(mat) <- mat$Gene
mat <- mat[, -1]


## ---------------------------------------------------------
## Heatmaps
## ---------------------------------------------------------

heatmap_Tajima <- pheatmap(
  t(mat[, "Tajima's D", drop = FALSE]),
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 22, cellheight = 20,
  color = colorRampPalette(c("#053061", "white", "#67001f"))(100),
  border_color = "gray80",
  main = "Tajima's D",
  show_colnames = TRUE, show_rownames = FALSE,
  silent = TRUE
)

heatmap_dNdS <- pheatmap(
  t(mat[, "dN/dS", drop = FALSE]),
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 22, cellheight = 20,
  color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
  border_color = "gray80",
  main = "dN/dS",
  show_colnames = TRUE, show_rownames = FALSE,
  silent = TRUE
)

heatmap_pi <- pheatmap(
  t(mat[, "π", drop = FALSE]),
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 22, cellheight = 20,
  color = colorRampPalette(c("white", "#E67E22", "#B03A2E"))(100),
  border_color = "gray80",
  main = "Nucleotide diversity (π)",
  show_colnames = TRUE, show_rownames = FALSE,
  silent = TRUE
)


## ---------------------------------------------------------
## Save Figure 2H
## ---------------------------------------------------------

CairoPDF(
  file.path(FIG_DIR_HEAT, "Figure2H_Selection_Heatmap.pdf"),
  width = 14, height = 8
)

grid.arrange(
  grobs = list(
    heatmap_Tajima[[4]],
    heatmap_dNdS[[4]],
    heatmap_pi[[4]]
  ),
  nrow = 3
)

dev.off()
