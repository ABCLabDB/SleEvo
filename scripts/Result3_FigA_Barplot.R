############################################################
# Result3 – NREM ratio evolutionary boxplot (simplified)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Biostrings)
library(ape)

PROJECT_DIR <- getwd()

FASTA_DIR <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE <- file.path(PROJECT_DIR, "data", "Metadata", "METADATA_fixedVersion.txt")
PERM_FILE <- file.path(PROJECT_DIR, "data", "Result3", "NREM_key_12_optimal_Kruskal.tsv")
OUT_DIR   <- file.path(PROJECT_DIR, "figures", "Result3", "Figure3_BD")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 1. Metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <-
  gsub(" ", "_", meta$Species_symbol_name_ensembl)

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

meta$Real_Class <- "Others"
meta$Real_Class[meta$Species_name_ensembl %in% top_species]    <- "High_Ratio"
meta$Real_Class[meta$Species_name_ensembl %in% bottom_species] <- "Low_Ratio"

## =========================================================
## 2. Significant genes
## =========================================================

sig_genes <- fread(PERM_FILE) |>
  filter(adj.P < 0.05)

sig_genes$Cluster <- as.numeric(sig_genes$Cluster)

## =========================================================
## 3. Main loop
## =========================================================

for (i in seq_len(nrow(sig_genes))) {
  
  gene <- sig_genes$Gene[i]
  k    <- sig_genes$Cluster[i]
  
  fasta_path <- file.path(FASTA_DIR,
                          paste0(gene, "_muscle.fasta"))
  
  if (!file.exists(fasta_path)) next
  
  cds <- readDNAStringSet(fasta_path)
  
  names(cds) <- sapply(strsplit(names(cds), ":"), `[`, 2)
  
  keep <- intersect(names(cds),
                    meta$Species_symbol_name_ensembl)
  
  cds <- cds[names(cds) %in% keep]
  
  names(cds) <- meta$Species_name_ensembl[
    match(names(cds),
          meta$Species_symbol_name_ensembl)
  ]
  
  use_meta <- meta[
    match(names(cds),
          meta$Species_name_ensembl), ]
  
  ## -------------------------
  ## phylogenetic clustering
  ## -------------------------
  
  dna <- as.DNAbin(cds)
  
  hc <- hclust(
    as.dist(cophenetic.phylo(
      njs(dist.dna(dna, as.matrix = TRUE,
                   pairwise.deletion = TRUE))
    )),
    method = "average"
  )
  
  clusters <- cutree(hc, k = k)
  
  cluster_df <- data.frame(
    Species = names(clusters),
    Cluster = clusters,
    Real_Class = use_meta$Real_Class
  )
  
  ## -------------------------
  ## majority class per cluster
  ## -------------------------
  
  majority_map <- cluster_df |>
    filter(Real_Class %in% c("High_Ratio", "Low_Ratio")) |>
    group_by(Cluster) |>
    count(Real_Class) |>
    slice_max(n, n = 1, with_ties = FALSE) |>
    select(Cluster, Majority = Real_Class)
  
  cluster_df <- cluster_df |>
    left_join(majority_map, by = "Cluster")
  
  cluster_df$Evolution_relation <-
    ifelse(cluster_df$Real_Class ==
             cluster_df$Majority,
           cluster_df$Real_Class,
           "no_relation")
  
  plot_df <- cluster_df |>
    filter(Evolution_relation %in%
             c("High_Ratio", "Low_Ratio")) |>
    left_join(use_meta,
              by = c("Species" = "Species_name_ensembl"))
  
  if (nrow(plot_df) < 4) next
  
  ## =========================================================
  ## 4. Plot
  ## =========================================================
  
  overall_mean <-
    mean(meta$Percentage_of_NREM_time_per_day)
  
  group_means <- plot_df |>
    group_by(Evolution_relation) |>
    summarise(mean_sleep =
                mean(Percentage_of_NREM_time_per_day))
  
  set.seed(43)
  
  plot_df$jitter_x <-
    as.numeric(as.factor(plot_df$Evolution_relation)) +
    runif(nrow(plot_df), -0.25, 0.25)
  
  p <- ggplot(plot_df) +
    geom_boxplot(
      aes(x = Evolution_relation,
          y = Percentage_of_NREM_time_per_day,
          fill = Evolution_relation),
      width = 0.5, outlier.shape = NA
    ) +
    geom_point(
      aes(x = jitter_x,
          y = Percentage_of_NREM_time_per_day,
          color = Evolution_relation),
      size = 4, alpha = 0.5
    ) +
    geom_text_repel(
      aes(x = jitter_x,
          y = Percentage_of_NREM_time_per_day,
          label = Species),
      size = 3
    ) +
    geom_hline(yintercept = overall_mean,
               linetype = "dashed",
               color = "grey40") +
    scale_fill_manual(values = c(
      "High_Ratio" = "#2c6e49",
      "Low_Ratio"  = "#f5d7b0"
    )) +
    scale_color_manual(values = c(
      "High_Ratio" = "#2c6e49",
      "Low_Ratio"  = "#f5d7b0"
    )) +
    theme_minimal() +
    labs(
      title = gene,
      x = "Sleep type",
      y = "NREM sleep (%)"
    )
  
  ggsave(file.path(OUT_DIR,
                   paste0(gene, "_Figure3BD_boxplot.pdf")),
         p, width = 8, height = 4)
  
  message("Saved: ", gene)
}
