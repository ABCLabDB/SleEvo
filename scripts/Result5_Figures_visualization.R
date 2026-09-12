############################################################
# Result5_03_Figures.R
#
# Generate all Result5 figure panels:
#   - Fig5A: Cochran test lollipop plot
#   - Fig5B/E: Phylogenetic tree (ATF5, NFIL3)
#   - Fig5C/D: Mosaic plot (cluster vs phenotype)
#
# Input:
#   data/Result1/species_sleep_metadata.txt
#   data/Result5/Sleep_frequency_Fisher.tsv
#   data/Fasta/<Gene>_muscle.fasta
#
# Output:
#   figures/Result5/*.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ape)
  library(Biostrings)
  library(dendextend)
  library(ggmosaic)
})

OUT_DIR <- "figures/Result5/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# --------------------- Fig5A ------------------------------
# Fisher enrichment lollipop plot
############################################################

cochran_df <- fread(
  "data/Result5/Sleep_frequency_Fisher.tsv"
) %>% as.data.frame()

cochran_df$log_Pvalue <- -log10(cochran_df$P_Value)
threshold <- -log10(0.05)

cochran_df <- cochran_df %>%
  mutate(
    Significant =
      ifelse(log_Pvalue > threshold,
             "Sleep frequency related",
             "Not significant")
  ) %>%
  arrange(desc(log_Pvalue)) %>%
  slice(1:32)

p_lollipop <- ggplot(
  cochran_df,
  aes(x = reorder(Gene, -log_Pvalue),
      y = log_Pvalue,
      color = Significant)
) +
  geom_segment(
    aes(xend = reorder(Gene, -log_Pvalue),
        y = 0,
        yend = log_Pvalue),
    linewidth = 5
  ) +
  geom_point(size = 5) +
  geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    color = "red",
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Sleep frequency related" = "#ffa600",
      "Not significant" = "grey85"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x =
      element_text(angle = 90,
                   hjust = 1,
                   vjust = 1,
                   size = 12, face = "italic)
  ) +
  labs(x = NULL,
       y = expression(-log[10](P)),
       color = "")

ggsave(
  file.path(OUT_DIR, "Figure5A_Lollipop.pdf"),
  p_lollipop,
  width = 8,
  height = 5
)

############################################################
# ---------------- Fig5B–E --------------------------------
# Trees + Mosaic (ATF5, NFIL3)
############################################################

metadata <- fread(
  "data/Result1/species_sleep_metadata.txt"
) %>% as.data.frame()

metadata$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_",
               metadata$Species_symbol_name_ensembl))

metadata <- metadata[
  !is.na(metadata$Number_of_sleep_times_per_day),
]

metadata$Group <- ifelse(
  metadata$Number_of_sleep_times_per_day == "Once",
  "Once",
  "More than twice"
)

perm_test <- fread(
  "data/Result5/Sleep_frequency_Fisher.tsv"
) %>% as.data.frame()

knee <- perm_test %>%
  filter(P_Value < 0.05)

knee <- knee[knee$Gene %in% c("ATF5", "NFIL3"), ]
rownames(knee) <- NULL

custom_colors <- c(
  "More than twice" = "#476066",
  "Once" = "#b8cdab"
)

FASTA_DIR <- "data/Fasta/"

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

  # ---- Tree ----

  dna <- as.DNAbin(cds)
  dm  <- dist.dna(dna, model="T92", pairwise.deletion=TRUE)
  tree <- nj(dm)
  hc   <- hclust(as.dist(cophenetic(tree)), method="average")
  dend <- as.dendrogram(hc)

  sleep_info <- setNames(
    use_meta$Group,
    use_meta$Species_name_ensembl
  )

  assign_branch_color_by_mode <- function(d) {
    if (is.leaf(d)) {
      sleep_type <- sleep_info[labels(d)]
      attr(d, "Group") <- sleep_type
      if (!is.na(sleep_type)) {
        attr(d, "edgePar") <- list(col = custom_colors[sleep_type], lwd = 4)
      } else {
        attr(d, "edgePar") <- list(col = "gray70", lwd = 4)
      }
      return(d)
    }
    
    d[[1]] <- assign_branch_color_by_mode(d[[1]])
    d[[2]] <- assign_branch_color_by_mode(d[[2]])
    
    left <- attr(d[[1]], "Group")
    right <- attr(d[[2]], "Group")
    all_types <- na.omit(c(left, right))
    
    if (length(all_types) == 0) {
      majority_type <- NA
      branch_color <- "gray70"
    } else {
      majority_type <- names(sort(table(all_types), decreasing = TRUE))[1]
      branch_color <- custom_colors[majority_type]
    }
    
    attr(d, "Group") <- majority_type
    attr(d, "edgePar") <- list(col = branch_color, lwd = 4)
    
    return(d)
  }

  dend_colored <- assign_branch_color_by_mode(dend)
  
  label_colors <- sleep_info[labels(dend_colored)]  
  label_colors <- custom_colors[label_colors] 
  label_colors[is.na(label_colors)] <- "gray70"  
  
  labels_colors(dend_colored) <- label_colors

  ggd <- dend_colored %>%
    set("labels_cex", 1) %>%
    set("leaves_pch", 19) %>%
    set("leaves_col", meta_for_color$color) %>%
    set("leaves_cex", 2.5)
  
  pdf(
    file.path(
      OUT_DIR,
      paste0("Figure5B_Tree_", gene_name, ".pdf")
    ),
    width = 8,
    height = 5
  )

  plot(ggd)
  legend("topright",
         legend = names(custom_colors),
         fill = custom_colors,
         title = "Sleep frequency",
         cex = 0.8)
  dev.off()

  # ---- Mosaic ----

  clusters <- cutree(hc, k = k)
  
  df <- data.frame(
    Species = names(clusters),
    Cluster = clusters
  ) %>%
    left_join(
      use_meta[,c("Species_name_ensembl","Group")],
      by=c("Species"="Species_name_ensembl")
    )
  
  cluster_majority <- df %>%
    group_by(Cluster) %>%
    summarise(
      Majority = names(which.max(table(Group))),
      .groups="drop"
    )
  
  df <- df %>%
    left_join(cluster_majority, by="Cluster")
  
  df_tbl <- as.data.frame(
    table(df$Group, df$Majority)
  )
  
  df_tbl <- df_tbl %>%
    group_by(Var1) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ungroup() %>%
    mutate(label = paste0(Freq, " Species\n(", round(prop*100, 1), "%)"),
           MatchStatus = ifelse(Var1 == Var2, "Evolution\ncorrelated", "Other effects"))
  
  p <- ggplot(df_tbl) +
    geom_mosaic(
      aes(x = product(Var1),
          fill = Var2,
          weight = Freq),
      color = "black", linewidth = 0.3
    ) +
    scale_fill_manual(values = c("#476066", "#b8cdab")) +
    labs(x = "", y = "") +
    theme_minimal(base_size = 14)+
    theme(legend.position = "top")
  
  gb <- ggplot_build(p)
  panel_data <- gb$data[[1]]
  
  hex_to_var <- setNames(df_tbl$Var2, unique(panel_data$fill))
  label_df <- panel_data %>%
    mutate(
      Var1 = x__Var1,
      Var2 = x__fill__Var2,
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2
    ) %>%
    left_join(df_tbl, by = c("Var1", "Var2"))
  
  label_df <- label_df %>%
    mutate(highlight = ifelse(Var1 == Var2, "Match", "Mismatch"))
  
  p_mosaic <- p + geom_text(
    data = label_df,
    aes(x = x, y = y, label = label.y),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold",
    color = "black"
  )
  
  ggsave(
    file.path(
      OUT_DIR,
      paste0("Figure5C_Mosaic_", gene_name, ".pdf")
    ),
    p_mosaic,
    width=7,
    height=5
  )
