############################################################
# Result4_Fig4B_SleepTiming_Tree.R
#
# Phylogenetic tree colored by sleep timing category
#
# Input:
#   data/Result4/Sleeptiming_Cochran_Result.tsv
#   data/Result1/species_sleep_metadata.txt
#   data/Fasta/*.fasta
#
# Output:
#   figures/Result4/Figure4B/<GENE>_SleepTimingTree.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Biostrings)
  library(ape)
  library(dendextend)
  library(stringr)
})

PROJECT_DIR <- getwd()

## =========================================================
## 1. Paths
## =========================================================

PERM_FILE <- file.path(PROJECT_DIR,
                       "data","Result4",
                       "Sleeptiming_Cochran_Result.tsv")

META_FILE <- file.path(PROJECT_DIR,
                       "data","Result1",
                       "species_sleep_metadata.txt")

FASTA_DIR <- file.path(PROJECT_DIR,
                       "data","Fasta")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures","Result4",
                     "Figure4B")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. Load metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <-
  tolower(gsub(" ","_",
               meta$Species_symbol_name_ensembl))

meta <- meta %>%
  filter(!is.na(Sleep_Timing))

meta$Type <- case_when(
  Sleep_Timing == "Sleep at night"   ~ "Sleep at night",
  Sleep_Timing == "Sleep at daytime" ~ "Sleep at daytime",
  TRUE                               ~ "Sleep at anytime"
)

custom_colors <- c(
  "Sleep at night"   = "#2c6e49",
  "Sleep at anytime" = "grey85",
  "Sleep at daytime" = "#d68c45"
)

## =========================================================
## 3. Load significant genes (Cochran P < 0.05)
## =========================================================

perm_test <- fread(PERM_FILE) |> as.data.frame()

perm_test$gene.cluster.idx <-
  sapply(strsplit(perm_test$gene.cluster.idx,"_"),
         `[`,2)

perm_test$gene.cluster.idx <-
  as.numeric(perm_test$gene.cluster.idx)

sig_genes <- perm_test %>%
  filter(P.cochran < 0.05)

## =========================================================
## 4. Helper function: build colored dendrogram
## =========================================================

build_tree_plot <- function(gene, k){

  fasta_file <- file.path(FASTA_DIR,
                          paste0(gene,"_muscle.fasta"))

  if(!file.exists(fasta_file)) return(NULL)

  cds <- readDNAStringSet(fasta_file)

  names(cds) <-
    sapply(strsplit(names(cds),":"),
           `[`,2)

  names(cds) <- str_to_title(names(cds))

  keep <- intersect(names(cds),
                    meta$Species_symbol_name_ensembl)

  cds <- cds[names(cds) %in% keep]

  names(cds) <- meta$Species_name_ensembl[
    match(names(cds),
          meta$Species_symbol_name_ensembl)
  ]

  dna <- as.DNAbin(cds)

  dm <- dist.dna(dna,
                 model="T92",
                 pairwise.deletion=TRUE)

  tree <- njs(dm)

  hc <- hclust(as.dist(cophenetic(tree)),
               method="average")

  dend <- as.dendrogram(hc)

  ## sleep mapping
  sleep_info <- setNames(
    meta$Type,
    meta$Species_name_ensembl
  )

  color_branches <- function(d){
    if(is.leaf(d)){
      sp <- labels(d)
      st <- sleep_info[sp]
      attr(d,"edgePar") <-
        list(col=custom_colors[st],
             lwd=4)
      return(d)
    }

    d[[1]] <- color_branches(d[[1]])
    d[[2]] <- color_branches(d[[2]])

    types <- na.omit(c(
      attr(d[[1]],"edgePar")$col,
      attr(d[[2]],"edgePar")$col
    ))

    if(length(types)==0){
      col_branch <- "gray70"
    } else {
      col_branch <- names(sort(table(types),
                               decreasing=TRUE))[1]
    }

    attr(d,"edgePar") <-
      list(col=col_branch,
           lwd=4)

    return(d)
  }

  dend <- color_branches(dend)

  return(dend)
}

## =========================================================
## 5. Loop genes
## =========================================================

for(i in seq_len(nrow(sig_genes))){

  gene <- sig_genes$gene.idx[i]
  k    <- sig_genes$gene.cluster.idx[i]

  if(k <= 1) next

  dend <- build_tree_plot(gene,k)

  if(is.null(dend)) next

  pdf(file.path(OUT_DIR,
                paste0(gene,
                       "_SleepTimingTree.pdf")),
      width=14,height=5)

  par(mar=c(8,4,2,2))

  plot(dend)

  legend("topright",
         legend=names(custom_colors),
         fill=custom_colors,
         title="Sleep Timing",
         cex=0.8)

  dev.off()
}

message("Figure4A trees saved.")
