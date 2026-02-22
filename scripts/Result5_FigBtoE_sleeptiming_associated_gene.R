############################################################
# Result5_FigBCDE_sleeptiming_associated_gene.R
#
# Sleep frequency–associated genes (ATF5 and NFIL3 only)
#
# This script generates:
#   - Fig5B/E: Phylogenetic tree colored by sleep frequency
#   - Fig5C/D: Mosaic plot (cluster majority vs phenotype)
#
# Target genes:
#   ATF5
#   NFIL3
#
# Input:
#   data/Result1/species_sleep_metadata.txt
#   data/Fasta/<Gene>_muscle.fasta
#   data/Result5/Sleep_frequency_Cochran.tsv
#
# Output:
#   figures/Result5/Figure5B_Tree_<Gene>.pdf
#   figures/Result5/Figure5C_Mosaic_<Gene>.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ape)
  library(Biostrings)
  library(dendextend)
  library(ggplot2)
  library(ggmosaic)
})

# ----------------------------------------------------------
# 1. Load species metadata
# ----------------------------------------------------------

metadata <- fread("data/Result1/species_sleep_metadata.txt") %>%
  as.data.frame()

metadata$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_", metadata$Species_symbol_name_ensembl))

metadata <- metadata[
  !is.na(metadata$Number_of_sleep_times_per_day),
]

metadata$Group <- ifelse(
  metadata$Number_of_sleep_times_per_day == "Once",
  "Once",
  "More than twice"
)

# ----------------------------------------------------------
# 2. Load Cochran test results
# ----------------------------------------------------------

perm_test <- fread("data/Result5/Sleep_frequency_Cochran.tsv") %>%
  as.data.frame()

knee <- perm_test %>%
  filter(P_Value < 0.05)

# ---- Restrict to ATF5 and NFIL3 only ----
knee <- knee[knee$Gene %in% c("ATF5", "NFIL3"), ]
rownames(knee) <- NULL

# ----------------------------------------------------------
# 3. Directories
# ----------------------------------------------------------

FASTA_DIR <- "data/Fasta/"
OUT_DIR   <- "figures/Result5/"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

custom_colors <- c(
  "More than twice" = "#476066",
  "Once" = "#b8cdab"
)

# ----------------------------------------------------------
# 4. Main loop (ATF5 and NFIL3 only)
# ----------------------------------------------------------

for(i in seq_len(nrow(knee))){

  gene_name <- knee$Gene[i]
  k <- knee$Cluster[i]

  if(k == 2) next

  fasta_file <- file.path(
    FASTA_DIR,
    paste0(gene_name, "_muscle.fasta")
  )

  cds <- readDNAStringSet(fasta_file)

  cds@ranges@NAMES <-
    sapply(strsplit(names(cds), ":"), `[`, 2)

  valid_species <- intersect(
    cds@ranges@NAMES,
    metadata$Species_symbol_name_ensembl
  )

  cds <- cds[cds@ranges@NAMES %in% valid_species]

  cds@ranges@NAMES <-
    metadata$Species_name_ensembl[
      match(cds@ranges@NAMES,
            metadata$Species_symbol_name_ensembl)
    ]

  use_meta <- metadata[
    match(cds@ranges@NAMES,
          metadata$Species_name_ensembl),
  ]

  # ------------------------------------------------------
  # Build phylogenetic tree (Fig5B/E)
  # ------------------------------------------------------

  dna <- as.DNAbin(cds)

  dm <- dist.dna(
    dna,
    model = "T92",
    pairwise.deletion = TRUE
  )

  tree <- nj(dm)
  hc   <- hclust(as.dist(cophenetic(tree)), method = "average")
  dend <- as.dendrogram(tree)

  sleep_info <- setNames(
    use_meta$Group,
    use_meta$Species_name_ensembl
  )

  assign_branch_color <- function(d){

    if(is.leaf(d)){
      sp <- labels(d)
      attr(d, "edgePar") <-
        list(
          col = custom_colors[sleep_info[sp]],
          lwd = 3
        )
      return(d)
    }

    d[[1]] <- assign_branch_color(d[[1]])
    d[[2]] <- assign_branch_color(d[[2]])
    return(d)
  }

  dend_colored <- assign_branch_color(dend)

  labels_colors(dend_colored) <-
    custom_colors[sleep_info[labels(dend_colored)]]

  ggd <- dend_colored %>%
    set("labels_cex", 0.9) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 2)

  pdf(
    file.path(
      OUT_DIR,
      paste0("Figure5B_Tree_", gene_name, ".pdf")
    ),
    width = 8,
    height = 5
  )

  plot(ggd)

  legend(
    "topright",
    legend = names(custom_colors),
    fill = custom_colors,
    title = "Sleep frequency",
    cex = 0.8
  )

  dev.off()

  # ------------------------------------------------------
  # Mosaic plot (Fig5C/D)
  # ------------------------------------------------------

  clusters <- cutree(hc, k = k)

  df <- data.frame(
    Species = names(clusters),
    Cluster = clusters
  ) %>%
    left_join(
      use_meta[, c("Species_name_ensembl","Group")],
      by = c("Species"="Species_name_ensembl")
    )

  cluster_majority <- df %>%
    group_by(Cluster) %>%
    summarise(
      Majority = names(which.max(table(Group))),
      .groups = "drop"
    )

  df <- df %>%
    left_join(cluster_majority, by = "Cluster")

  df_tbl <- as.data.frame(
    table(df$Group, df$Majority)
  )

  df_tbl <- df_tbl %>%
    group_by(Var1) %>%
    mutate(
      prop  = Freq / sum(Freq),
      label = paste0(Freq, " (", round(prop*100,1), "%)")
    ) %>%
    ungroup()

  p_mosaic <- ggplot(df_tbl) +
    geom_mosaic(
      aes(
        x = product(Var1),
        fill = Var2,
        weight = Freq
      ),
      color = "black"
    ) +
    geom_mosaic_text(
      aes(
        x = product(Var1),
        fill = Var2,
        weight = Freq,
        label = label
      ),
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top") +
    labs(
      x = "",
      y = "",
      fill = "Cluster majority"
    )

  ggsave(
    file.path(
      OUT_DIR,
      paste0("Figure5C_Mosaic_", gene_name, ".pdf")
    ),
    p_mosaic,
    width = 7,
    height = 5
  )
}
