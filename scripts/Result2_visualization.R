############################################################
# Result2: Visualization of evolutionary association
#
# This script contains visualization code for Result2,
# which links phylogenetic clustering of circadian genes
# with total sleep time across species.
#
# The script consists of two main parts:
#
# Part A. Phylogenetic tree visualization
#   - Gene-wise dendrograms based on CDS sequence similarity
#   - Branch colors represent majority sleep phenotype
#   - Used for Figure 2A and related supplementary figures
#
# Part B. Box/strip plot of sleep phenotypes
#   - Compares total sleep time between evolutionarily
#     matched Long-sleep vs Short-sleep groups
#   - Used for Figure 2B
#
############################################################


## =========================================================
## 1. Libraries
## =========================================================
library(data.table)
library(dplyr)
library(stringr)
library(Biostrings)
library(ape)
library(dendextend)
library(ggdendro)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)


## =========================================================
## 2. Project paths
## =========================================================
PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1", "species_sleep_metadata.txt")
STAT_FILE  <- file.path(PROJECT_DIR, "data", "Result2", "TST_key_12_optimal_Kruskal.tsv")

FIG_DIR_TREE <- file.path(PROJECT_DIR, "figures", "Result2", "PhyloTree")
FIG_DIR_BOX  <- file.path(PROJECT_DIR, "figures", "Result2", "Boxplot")

dir.create(FIG_DIR_TREE, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_BOX,  recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Load metadata and define sleep categories
## =========================================================
meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

meta <- meta |> filter(!is.na(Total_sleep_time_per_day)) |>
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
  filter(P_anova < 0.05)

stat$Cluster <- as.numeric(stat$Cluster)


## =========================================================
## Figure2A,C,D. Phylogenetic tree visualization
## =========================================================

sleep_colors <- c(
  "Long_sleep"  = "#194a7a",
  "Others"     = "grey70",
  "Short_sleep"= "#E7B800"
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
  cds <- cds[names(cds) %in% keep]

  names(cds) <- meta$Species_name_ensembl[
    match(names(cds), meta$Species_symbol_name_ensembl)
  ]

  use_meta <- meta[match(names(cds), meta$Species_name_ensembl),]

  dna <- as.DNAbin(cds)
  tree <- njs(dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE))

  dend <- as.dendrogram(hclust(as.dist(cophenetic.phylo(tree)), "average"))

  sleep_info <- setNames(use_meta$Type, use_meta$Species_name_ensembl)
  dend_colored <- color_tree_by_sleep(dend, sleep_info)

  pdf(file.path(FIG_DIR_TREE, paste0(gene, "_phylo_tree.pdf")), 12, 6)
  plot(dend_colored, main = gene)
  legend("topright", legend = names(sleep_colors),
         fill = sleep_colors, title = "Sleep Type", cex = 0.8)
  dev.off()
}


## =========================================================
## Figure 2B. Box / strip plot of sleep phenotypes (Figure 2B)
##
## This section compares total sleep time between
## evolutionarily matched Long_sleep and Short_sleep groups
## defined based on phylogenetic clustering.
## =========================================================

for (i in seq_len(nrow(stat))) {

  gene <- stat$Gene[i]
  k    <- stat$Cluster[i]
  if (k == 2) next

  message("Boxplot for gene: ", gene)

  fasta_path <- file.path(DATA_FASTA, paste0(gene, "_muscle.fasta"))
  if (!file.exists(fasta_path)) next

  ## ---------------------------------------
  ## Reconstruct clustering
  ## ---------------------------------------
  cds <- readDNAStringSet(fasta_path)
  names(cds) <- sapply(strsplit(names(cds), ":"), `[`, 2)
  names(cds) <- str_to_title(names(cds))

  keep <- intersect(names(cds), meta$Species_symbol_name_ensembl)
  cds  <- cds[names(cds) %in% keep]

  names(cds) <- meta$Species_name_ensembl[
    match(names(cds), meta$Species_symbol_name_ensembl)
  ]

  use_meta <- meta[match(names(cds), meta$Species_name_ensembl), ]
  rownames(use_meta) <- NULL

  dna <- as.DNAbin(cds)
  tree <- njs(dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE))
  hc <- hclust(as.dist(cophenetic.phylo(tree)), method = "average")
  clusters <- cutree(hc, k = k)

  ## ---------------------------------------
  ## Assign evolutionary classes
  ## ---------------------------------------
  class_assign <- rep(NA, length(clusters))
  names(class_assign) <- names(clusters)

  class_assign[names(clusters) %in% top_species]    <- "Long_sleep"
  class_assign[names(clusters) %in% bottom_species] <- "Short_sleep"

  df <- data.frame(
    Species = names(class_assign),
    Cluster = clusters,
    Class   = class_assign,
    stringsAsFactors = FALSE
  ) |> filter(!is.na(Class))

  ## ---------------------------------------
  ## Match inferred vs real sleep class
  ## ---------------------------------------
  sleep_df <- use_meta
  sleep_df$Infer_Class <- "Others"
  sleep_df$Infer_Class[sleep_df$Species_name_ensembl %in% df$Species[df$Class == "Long_sleep"]]  <- "Long_sleep"
  sleep_df$Infer_Class[sleep_df$Species_name_ensembl %in% df$Species[df$Class == "Short_sleep"]] <- "Short_sleep"

  sleep_df$Real_Class <- "Others"
  sleep_df$Real_Class[sleep_df$Species_name_ensembl %in% top_species]    <- "Long_sleep"
  sleep_df$Real_Class[sleep_df$Species_name_ensembl %in% bottom_species] <- "Short_sleep"

  sleep_df <- sleep_df |>
    mutate(Evolution_relation =
             if_else(Infer_Class == Real_Class, Infer_Class, "no_match"))

  plot_df <- sleep_df |>
    filter(Evolution_relation %in% c("Long_sleep", "Short_sleep"))

  if (nrow(plot_df) < 4) next

  ## ---------------------------------------
  ## Statistics for plotting
  ## ---------------------------------------
  overall_mean <- mean(meta$Total_sleep_time_per_day, na.rm = TRUE)

  group_means <- plot_df |>
    group_by(Evolution_relation) |>
    summarise(mean_sleep = mean(Total_sleep_time_per_day, na.rm = TRUE),
              .groups = "drop")

  set.seed(42)
  plot_df <- plot_df |>
    mutate(
      jitter_x = as.numeric(factor(Evolution_relation)) +
        runif(n(), -0.25, 0.25)
    )

  ## ---------------------------------------
  ## Plot
  ## ---------------------------------------
  p <- ggplot() +

    geom_boxplot(
      data = plot_df,
      aes(x = Evolution_relation, y = Total_sleep_time_per_day,
          fill = Evolution_relation),
      width = 0.5, outlier.shape = NA, color = "black"
    ) +

    geom_point(
      data = plot_df,
      aes(x = jitter_x, y = Total_sleep_time_per_day,
          color = Evolution_relation),
      size = 4, alpha = 0.6
    ) +

    geom_text_repel(
      data = plot_df,
      aes(x = jitter_x, y = Total_sleep_time_per_day,
          label = Species_name_ensembl),
      size = 3.5, max.overlaps = 4
    ) +

    geom_hline(
      yintercept = overall_mean,
      linetype = "dashed", color = "grey40"
    ) +

    scale_fill_manual(values = sleep_colors) +
    scale_color_manual(values = sleep_colors) +

    labs(
      title = gene,
      x = "Evolutionary sleep class",
      y = "Total sleep time (hours per day)"
    ) +

    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, 23.5))

  ggsave(
    file.path(FIG_DIR_BOX, paste0(gene, "_boxplot.pdf")),
    plot = p, width = 10, height = 4
  )


  ## =========================================================
## Figure 2E. Phylogenetic signal across genes
##
## This section visualizes phylogenetic signal statistics
## (Blomberg’s K and Moran’s I) for circadian genes.
##
## Input:
##  - Permutation-based phylogenetic signal results
##
## Output:
##  - Line / point plot summarizing signal strength per gene
##
## Used for:
##  - Figure 2E
## =========================================================

SIGNAL_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "Result2",
  "Phylogenetic_Signal_Data_permutation.tsv"
)

FIG_DIR_SIGNAL <- file.path(
  PROJECT_DIR,
  "figures",
  "Result2",
  "Phylogenetic_Signal"
)

dir.create(FIG_DIR_SIGNAL, recursive = TRUE, showWarnings = FALSE)


## -----------------------------------------
## Load and reshape data
## -----------------------------------------
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


## -----------------------------------------
## Plot
## -----------------------------------------
p_signal <- ggplot(signal_long,
                   aes(x = Gene, y = Value, color = Signal)) +

  geom_line(aes(group = Signal),
            linewidth = 1.2, alpha = 0.4, color = "grey60") +

  geom_point(size = 3, alpha = 0.9) +

  geom_smooth(method = "lm", se = TRUE,
              linewidth = 1.5, color = "#95be8d") +

  scale_color_manual(
    values = c(
      "Blomberg’s K" = "#7593af",
      "Moran’s I"    = "#d69e49"
    )
  ) +

  theme_bw(base_size = 12) +

  labs(
    x = NULL,
    y = "Phylogenetic signal score",
    color = "Phylogenetic signal"
  ) +

  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 10)
  )


ggsave(
  file.path(FIG_DIR_SIGNAL, "Figure2C_Phylogenetic_Signal.pdf"),
  plot = p_signal,
  width = 6,
  height = 3.5
)
