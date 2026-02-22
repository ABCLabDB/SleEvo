############################################################
# Result5_FigBE_sleeptiming_associated_gene.R
#
# Phylogenetic tree colored by sleep frequency phenotype
# (Sleep frequency–associated genes)
#
# Input:
#   data/Result1/species_sleep_metadata.txt
#   data/Fasta/*.fasta
#   data/Result5/Sleep_frequency_Cochran.tsv
#
# Output:
#   figures/Result5/Figure5B_Tree_<Gene>.pdf
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
})

# ----------------------------------------------------------
# 1. Load metadata  (첫번째 그림 경로)
# ----------------------------------------------------------

metadata <- fread(
  "data/Result1/species_sleep_metadata.txt"
) %>% as.data.frame()

metadata$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_",
               metadata$Species_symbol_name_ensembl))

metadata <- metadata[
  !is.na(metadata$Number_of_sleep_times_per_day), ]

metadata$Type <- ifelse(
  metadata$Number_of_sleep_times_per_day == "Once",
  "Once",
  "More than twice"
)

# ----------------------------------------------------------
# 2. Load Cochran result  (세번째 그림 경로)
# ----------------------------------------------------------

perm_test <- fread(
  "data/Result5/Sleep_frequency_Cochran.tsv"
) %>% as.data.frame()

knee <- perm_test %>%
  filter(P_Value < 0.05)

target_genes <- knee$Gene

# ----------------------------------------------------------
# 3. Fasta directory  (두번째 그림 경로)
# ----------------------------------------------------------

FASTA_DIR <- "data/Fasta/"
fasta_files <- list.files(
  FASTA_DIR,
  pattern = "_muscle.fasta",
  full.names = TRUE
)

fasta_files <- fasta_files[
  tools::file_path_sans_ext(
    gsub("_muscle", "",
         basename(fasta_files))
  ) %in% target_genes
]

# custom color
custom_colors <- c(
  "More than twice" = "#476066",
  "Once" = "#b8cdab"
)

OUT_DIR <- "figures/Result5/"
dir.create(OUT_DIR,
           recursive = TRUE,
           showWarnings = FALSE)

# ----------------------------------------------------------
# 4. Main loop
# ----------------------------------------------------------

for(file in fasta_files){

  gene_name <- gsub(
    "_muscle",
    "",
    tools::file_path_sans_ext(basename(file))
  )

  cds_muscle <- readDNAStringSet(file)

  cds_muscle@ranges@NAMES <-
    sapply(strsplit(cds_muscle@ranges@NAMES, ":"),
           `[`, 2)

  # intersect species
  valid_species <- intersect(
    cds_muscle@ranges@NAMES,
    metadata$Species_symbol_name_ensembl
  )

  cds_muscle <- cds_muscle[
    cds_muscle@ranges@NAMES %in% valid_species
  ]

  cds_muscle@ranges@NAMES <-
    metadata$Species_name_ensembl[
      match(cds_muscle@ranges@NAMES,
            metadata$Species_symbol_name_ensembl)
    ]

  use_meta <- metadata[
    match(cds_muscle@ranges@NAMES,
          metadata$Species_name_ensembl),
  ]

  # distance matrix
  dna_Muscle <- as.DNAbin(cds_muscle)
  dm <- dist.dna(dna_Muscle,
                 model="T92",
                 pairwise.deletion=TRUE)

  tree <- nj(dm)
  dend <- as.dendrogram(tree)

  # ------------------------------------------------------
  # Color mapping
  # ------------------------------------------------------

  sleep_info <- setNames(
    use_meta$Type,
    use_meta$Species_name_ensembl
  )

  assign_branch_color <- function(d){

    if(is.leaf(d)){
      sp <- labels(d)
      sleep_type <- sleep_info[sp]
      attr(d, "edgePar") <-
        list(col = custom_colors[sleep_type],
             lwd = 3)
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

  # ------------------------------------------------------
  # Save
  # ------------------------------------------------------

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
         title = "Number of Sleep",
         cex = 0.8)

  dev.off()
}
