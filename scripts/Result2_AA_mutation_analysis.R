############################################################
# Result2_AA_mutation_analysis.R
#
# Amino-acid mutation analysis for sleep-associated SNPs
# across circadian genes.
#
# Output:
#   data/Result2/AminoAcid_Mutation_DF.tsv
############################################################


## =========================================================
## 1. Libraries
## =========================================================

library(data.table)
library(dplyr)
library(stringr)
library(Biostrings)


## =========================================================
## 2. Paths
## =========================================================

PROJECT_DIR <- getwd()

META_FILE <- file.path(PROJECT_DIR,
                       "data/Result1/species_sleep_metadata.txt")

PERM_FILE <- file.path(PROJECT_DIR,
                       "data/Result2/clustering.permutation200.Sleep_real.tsv")

SNP_FILE <- file.path(PROJECT_DIR,
                      "data/Result2/ALL_SNP.tsv")

SEQ_DIR <- file.path(PROJECT_DIR,
                     "data/Result2/NucleotideMatrix")

OUT_FILE <- file.path(PROJECT_DIR,
                      "data/Result2/AminoAcid_Mutation_DF.tsv")


## =========================================================
## 3. Load metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()
meta$Species_symbol_name_ensembl <- meta$Species_symbol_name_ensembl |>
  gsub(" ", "_", x = _) |>
  tolower()


## =========================================================
## 4. Significant genes & SNPs
## =========================================================

perm_test <- fread(PERM_FILE) |> as.data.frame()
sig_genes <- perm_test$gene.idx[perm_test$P.fisher < 0.05]

all_snp <- fread(SNP_FILE) |> as.data.frame()
all_snp <- all_snp |> filter(grepl("CDS", Column_Index))
sleep_genes <- intersect(unique(all_snp$Gene), sig_genes)


## =========================================================
## 5. Utility functions
## =========================================================

get_codons <- function(cds_row) {
  seq <- as.character(unlist(cds_row))
  sapply(seq(seq_len(length(seq) - 2), by = 3),
         function(i) paste(seq[i:(i + 2)], collapse = ""))
}

translate_codons <- function(codons) {
  as.character(translate(DNAStringSet(codons)))
}


get_aa_properties <- function() {
  data.frame(
    Amino_Acid = c("A","R","N","D","C","E","Q","G","H","I",
                   "L","K","M","F","P","S","T","W","Y","V","*"),
    Chemical_Property = c("Non-polar","Basic","Polar","Acidic","Polar",
                          "Acidic","Polar","Non-polar","Basic","Non-polar",
                          "Non-polar","Basic","Non-polar","Non-polar",
                          "Non-polar","Polar","Polar","Non-polar",
                          "Polar","Non-polar","Stop")
  )
}


## =========================================================
## 6. Main analysis
## =========================================================

aa_mutation_df <- list()

for (gene in sleep_genes) {

  seq_file <- file.path(SEQ_DIR, paste0(gene, ".tsv"))
  if (!file.exists(seq_file)) next

  aln <- fread(seq_file) |> as.data.frame()
  rownames(aln) <- aln[,1]
  aln <- aln[,-1]

  cds <- aln[, grepl("CDS", colnames(aln)), drop = FALSE]
  rownames(cds) <- str_to_title(rownames(cds))

  human_cds <- cds["Homo_sapiens", , drop = FALSE]
  orig_codons <- get_codons(human_cds)
  orig_aa <- translate_codons(orig_codons)

  snp_positions <- all_snp$Column_Index[all_snp$Gene == gene]

  for (pos in snp_positions) {

    if (!pos %in% colnames(cds)) next

    mutated_cds <- human_cds
    mutated_cds[1, pos] <- setdiff(unique(cds[,pos]), human_cds[1,pos])[1]

    mut_codons <- get_codons(mutated_cds)
    mut_aa <- translate_codons(mut_codons)

    diff_idx <- which(orig_codons != mut_codons)
    if (length(diff_idx) == 0) next

    aa_mutation_df[[length(aa_mutation_df) + 1]] <-
      data.frame(
        Gene = gene,
        Codon_Position = diff_idx,
        Original_Codon = orig_codons[diff_idx],
        Mutated_Codon = mut_codons[diff_idx],
        Original_Amino_Acid = orig_aa[diff_idx],
        Mutated_Amino_Acid = mut_aa[diff_idx],
        Mutation_Type = ifelse(orig_aa[diff_idx] == mut_aa[diff_idx],
                               "Synonymous", "Non-synonymous"),
        Full_Length = length(orig_aa)
      )
  }
}

aa_mutation_df <- bind_rows(aa_mutation_df)


## =========================================================
## 7. Annotate chemical properties
## =========================================================

aa_prop <- get_aa_properties()

aa_mutation_df <- aa_mutation_df |>
  left_join(aa_prop, by = c("Original_Amino_Acid" = "Amino_Acid")) |>
  rename(Original_Chemical_Property = Chemical_Property) |>
  left_join(aa_prop, by = c("Mutated_Amino_Acid" = "Amino_Acid")) |>
  rename(Mutated_Chemical_Property = Chemical_Property) |>
  mutate(Property_Change =
           ifelse(Original_Chemical_Property ==
                    Mutated_Chemical_Property,
                  "Unchanged", "Changed"))


## =========================================================
## 8. Save
## =========================================================

fwrite(aa_mutation_df, OUT_FILE, sep = "\t")
