############################################################
# Result2: Phylogenetic tree visualization
#
# This script generates phylogenetic dendrograms for genes
# significantly associated with total sleep time.
#
# Input:
#  - MUSCLE-aligned CDS FASTA files
#  - Species-level sleep metadata
#  - Result2 association results (TST_key_12_optimal_Kruskal.tsv)
#
# Output:
#  - Gene-wise phylogenetic tree plots (PDF)
#
# Figure:
#  - Main Figure 2A / Supplementary figures
############################################################


## =========================================================
## 1. Libraries
## =========================================================
library(data.table)
library(dplyr)
library(Biostrings)
library(ape)
library(dendextend)
library(ggdendro)
library(RColorBrewer)


## =========================================================
## 2. Project paths
## =========================================================
PROJECT_DIR <- getwd()

DATA_FASTA <- file.path(PROJECT_DIR, "data", "Fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1", "species_sleep_metadata.txt")
STAT_FILE  <- file.path(PROJECT_DIR, "data", "Result2", "TST_key_12_optimal_Kruskal.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "figures", "Result2", "PhyloTree")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 3. Load metadata
## =========================================================
meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()

meta <- meta |> arrange(Total_sleep_time_per_day)
rownames(meta) <- NULL

## Define sleep categories (bottom/top 33%)
prop_rate <- 0.33

bottom_species <- meta$Species_name_ensembl[
  1:ceiling(prop_rate * nrow(meta))
]

top_species <- meta$Species_name_ensembl[
  (nrow(meta) - ceiling(prop_rate * nrow(meta)) + 1):nrow(meta)
]

meta$Type <- "Others"
meta$Type[meta$Species_name_ensembl %in% bottom_species] <- "Short Sleep"
meta$Type[meta$Species_name_ensembl %in% top_species]    <- "Long Sleep"


## =========================================================
## 4. Load Result2 statistics
## =========================================================
stat <- fread(STAT_FILE) |> as.data.frame()
stat <- stat |> filter(P_anova < 0.05)
rownames(stat) <- NULL


## =========================================================
## 5. Color scheme
## =========================================================
sleep_colors <- c(
  "Long Sleep"  = "#194a7a",
  "Others"     = "grey70",
  "Short Sleep"= "#E7B800"
)


## =========================================================
## 6. Helper function: color branches by majority phenotype
## =========================================================
color_tree_by_sleep <- function(dend, sleep_info) {

  assign_color <- function(d) {

    if (is.leaf(d)) {
      st <- sleep_info[labels(d)]
      attr(d, "Type") <- st
      attr(d, "edgePar") <- list(
        col = sleep_colors[st],
        lwd = 3
      )
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


## =========================================================
## 7. Main loop: tree plotting
## =========================================================
for (i in seq_len(nrow(stat))) {

  gene <- stat$Gene[i]
  k    <- stat$Cluster[i]

  if (k == 2) next  # skip trivial clustering

  message("Plotting gene: ", gene)

  fasta_path <- file.path(DATA_FASTA, paste0(gene, "_muscle.fasta"))
  if (!file.exists(fasta_path)) next

  cds <- readDNAStringSet(fasta_path)
  names(cds) <- sapply(strsplit(names(cds), ":"), `[`, 2)

  keep <- intersect(names(cds), meta$Species_symbol_name_ensembl)
  cds <- cds[names(cds) %in% keep]

  names(cds) <- meta$Species_name_ensembl[
    match(names(cds), meta$Species_symbol_name_ensembl)
  ]

  use_meta <- meta[
    match(names(cds), meta$Species_name_ensembl),
  ]
  rownames(use_meta) <- NULL

  ## Build tree
  dna <- as.DNAbin(cds)
  dm  <- dist.dna(dna, as.matrix = TRUE, pairwise.deletion = TRUE)
  tree <- njs(dm)

  hc <- hclust(as.dist(cophenetic.phylo(tree)), method = "average")
  dend <- as.dendrogram(hc)

  ## Prepare colors
  sleep_info <- setNames(use_meta$Type, use_meta$Species_name_ensembl)
  dend_colored <- color_tree_by_sleep(dend, sleep_info)

  labels_colors(dend_colored) <- sleep_colors[sleep_info[labels(dend_colored)]]

  ## Plot
  pdf(
    file.path(OUT_DIR, paste0(gene, "_phylo_tree.pdf")),
    width = 12, height = 6
  )

  par(mar = c(10, 4, 2, 2))
  plot(dend_colored, main = gene)

  legend(
    "topright",
    legend = names(sleep_colors),
    fill   = sleep_colors,
    title  = "Sleep Type",
    cex    = 0.8
  )

  dev.off()
}
