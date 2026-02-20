############################################################
# Result3_Figure3A_NREM_phylo_tree.R
#
# Visualization of phylogenetic trees for
# NREM ratio–associated circadian genes.
#
# Input:
#   data/Result3/NREM_key_12_optimal_Kruskal.tsv
#   data/Fasta/*.fasta
#   data/Metadata/species_sleep_metadata.txt
#
# Output:
#   figures/Result3/Figure3A/*.pdf
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
  library(Biostrings)
  library(ape)
  library(dendextend)
})

## =========================================================
## 2. Project paths (relative to root)
## =========================================================

PROJECT_DIR <- getwd()

META_FILE <- file.path(PROJECT_DIR,
                       "data", "Result1",
                       "species_sleep_metadata.txt")

PERM_FILE <- file.path(PROJECT_DIR,
                       "data", "Result3",
                       "NREM_key_12_optimal_Kruskal.tsv")

FASTA_DIR <- file.path(PROJECT_DIR,
                       "data", "Fasta")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures", "Result3",
                     "Figure3A")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Load metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_", meta$Species_symbol_name_ensembl))

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

custom_colors <- c(
  "High Ratio" = "#194a7a",
  "Others"     = "grey70",
  "Low Ratio"  = "#E7B800"
)

## =========================================================
## 4. Load significant genes
## =========================================================

perm_test <- fread(PERM_FILE) |> as.data.frame()

sig_genes <- perm_test |>
  filter(adj.P < 0.05)

sig_genes$Cluster <- as.numeric(sig_genes$Cluster)

if (nrow(sig_genes) == 0) {
  stop("No significant genes found (adj.P < 0.05)")
}

## =========================================================
## 5. Function: color branches by majority sleep type
## =========================================================

color_tree_by_sleep <- function(dend, sleep_info) {
  
  assign_color <- function(d) {
    
    if (is.leaf(d)) {
      st <- sleep_info[labels(d)]
      attr(d, "Type") <- st
      attr(d, "edgePar") <- list(
        col = ifelse(is.na(st), "grey70", custom_colors[st]),
        lwd = 3
      )
      return(d)
    }
    
    d[[1]] <- assign_color(d[[1]])
    d[[2]] <- assign_color(d[[2]])
    
    types <- na.omit(c(attr(d[[1]], "Type"),
                       attr(d[[2]], "Type")))
    
    maj <- if (length(types) == 0)
      NA else names(sort(table(types), TRUE))[1]
    
    attr(d, "Type") <- maj
    attr(d, "edgePar") <- list(
      col = ifelse(is.na(maj), "grey70", custom_colors[maj]),
      lwd = 3
    )
    
    d
  }
  
  assign_color(dend)
}

## =========================================================
## 6. Main loop
## =========================================================

for (i in seq_len(nrow(sig_genes))) {
  
  gene <- sig_genes$Gene[i]
  
  fasta_path <- file.path(FASTA_DIR,
                          paste0(gene, "_muscle.fasta"))
  
  if (!file.exists(fasta_path)) {
    message("Skipping: FASTA not found for ", gene)
    next
  }
  
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
  
  dna <- as.DNAbin(cds)
  
  dm <- dist.dna(dna,
                 as.matrix = TRUE,
                 pairwise.deletion = TRUE)
  
  tree <- njs(dm)
  
  hc <- hclust(as.dist(cophenetic.phylo(tree)),
               method = "average")
  
  dend <- as.dendrogram(hc)
  
  sleep_info <- setNames(use_meta$Type,
                         use_meta$Species_name_ensembl)
  
  dend_colored <- color_tree_by_sleep(dend,
                                      sleep_info)
  
  pdf(file.path(OUT_DIR,
                paste0(gene, "_Figure3A_phylo_tree.pdf")),
      width = 12, height = 6)
  
  plot(dend_colored,
       main = paste0(gene,
                     " (NREM ratio-associated)"))
  
  legend("topright",
         legend = names(custom_colors),
         fill = custom_colors,
         title = "NREM Ratio Type")
  
  dev.off()
  
  message("Saved: ", gene)
}

## =========================================================
## End of script
## =========================================================
