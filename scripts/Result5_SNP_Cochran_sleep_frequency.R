############################################################
# Result5_SNP_Cochran_sleep_frequency.R
#
# Identification of sleep frequency–associated SNPs
# using Cochran–Armitage trend test
#
# Input:
#   data/Result1/species_sleep_metadata.txt
#   data/Result5/Sleep_frequency_Cochran.tsv
#   data/nucleotideDF/*.tsv
#
# Output:
#   figures/Result5/SNP_Significant.tsv
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(DescTools)
})

# ----------------------------------------------------------
# 1. Load species metadata
# ----------------------------------------------------------

species_metadata <- fread(
  "data/Result1/species_sleep_metadata.txt"
) %>% as.data.frame()

species_metadata$Species_symbol_name_ensembl <-
  gsub(" ", "_",
       species_metadata$Species_symbol_name_ensembl)

# ----------------------------------------------------------
# 2. Load significant genes (Cochran gene-level result)
# ----------------------------------------------------------

perm_test <- fread(
  "data/Result5/Sleep_frequency_Cochran.tsv"
) %>% as.data.frame()

perm_test <- perm_test %>%
  filter(P_Value < 0.05)

target_genes <- perm_test$Gene

# ----------------------------------------------------------
# 3. Load nucleotide matrices
# ----------------------------------------------------------

NUC_DIR <- "data/nucleotideDF/"

seq_files <- list.files(
  NUC_DIR,
  pattern = ".tsv$",
  full.names = TRUE
)

file_gene_names <- tools::file_path_sans_ext(
  basename(seq_files)
)

seq_files <- seq_files[
  file_gene_names %in% target_genes
]

# ----------------------------------------------------------
# 4. Cochran–Armitage test per gene
# ----------------------------------------------------------

TS_site <- list()

for(file in seq_files){

  gene_name <- tools::file_path_sans_ext(
    basename(file)
  )

  aln_df <- fread(file) %>% as.data.frame()

  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]

  # keep only CDS columns
  keep_cols <- grepl("CDS|V", colnames(aln_df))
  aln_df <- aln_df[, keep_cols]

  rownames(aln_df) <- str_to_title(rownames(aln_df))

  aln_df$Number_of_sleep_times_per_day <-
    species_metadata$Number_of_sleep_times_per_day[
      match(rownames(aln_df),
            species_metadata$Species_symbol_name_ensembl)
    ]

  aln_df <- aln_df[
    !is.na(aln_df$Number_of_sleep_times_per_day),
  ]

  # remove monomorphic sites
  unique_counts <- sapply(
    aln_df[,-ncol(aln_df)],
    function(x) length(unique(x))
  )

  polymorphic_cols <- which(unique_counts > 1)

  if(length(polymorphic_cols) == 0) next

  nuc_matrix <- aln_df[, c(polymorphic_cols, ncol(aln_df))]
  colnames(nuc_matrix)[ncol(nuc_matrix)] <-
    "Number_of_sleep_times_per_day"

  pvals <- sapply(
    colnames(nuc_matrix)[-ncol(nuc_matrix)],
    function(col){

      tab <- table(
        nuc_matrix[[col]],
        nuc_matrix$Number_of_sleep_times_per_day
      )

      CochranArmitageTest(tab)$p.value
    }
  )

  result_df <- data.frame(
    Gene = gene_name,
    Nucleotide_position = colnames(nuc_matrix)[-ncol(nuc_matrix)],
    Gene_length = ncol(aln_df) - 1,
    P = pvals
  )

  result_df$adjP <- p.adjust(result_df$P, method = "BH")

  TS_site[[gene_name]] <- result_df
}

TS_site <- bind_rows(TS_site)

# ----------------------------------------------------------
# 5. Significant SNPs
# ----------------------------------------------------------

TS_site_sig <- TS_site %>%
  filter(adjP < 0.05)

# ----------------------------------------------------------
# 6. Save result
# ----------------------------------------------------------

OUT_DIR <- "figures/Result5/"
dir.create(OUT_DIR,
           recursive = TRUE,
           showWarnings = FALSE)

fwrite(
  TS_site_sig,
  file.path(OUT_DIR, "SNP_Significant.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
