############################################################
# Result5_AA_effect_from_SNP.R
#
# Evaluate amino acid–level effects of sleep frequency–
# associated SNPs (Result5)
#
# Steps:
#   1. Load significant SNP sites
#   2. Extract CDS alignment for each gene
#   3. Replace nucleotide at SNP position (human reference)
#   4. Translate codons
#   5. Classify synonymous vs non-synonymous
#   6. Annotate amino acid chemical property change
#
# Input:
#   data/Result1/species_sleep_metadata.txt
#   data/Result5/SNP_Significant.tsv
#   data/nucleotideDF/<Gene>.tsv
#
# Output:
#   figures/Result5/AminoAcid_Mutation_DF.tsv
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Biostrings)
})

# ----------------------------------------------------------
# 1. Load metadata
# ----------------------------------------------------------

metadata <- fread("data/Result1/species_sleep_metadata.txt") %>%
  as.data.frame()

metadata$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_",
               metadata$Species_symbol_name_ensembl))

metadata <- metadata[
  !is.na(metadata$Number_of_sleep_times_per_day),
]

# ----------------------------------------------------------
# 2. Load significant SNPs
# ----------------------------------------------------------

TS_site <- fread("data/Result5/SNP_Significant.tsv") %>%
  as.data.frame()

TS_site <- TS_site %>%
  filter(grepl("CDS", Nucleotide_position))

sleep_genes <- unique(TS_site$Gene)

# ----------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------

get_codons <- function(dna_vec){
  sapply(seq(1, length(dna_vec) - 2, by = 3),
         function(i) paste(dna_vec[i:(i+2)], collapse=""))
}

translate_codons <- function(codon_seq){
  as.character(
    Biostrings::translate(DNAStringSet(codon_seq))
  )
}

# ----------------------------------------------------------
# 4. Main analysis
# ----------------------------------------------------------

NUC_DIR <- "data/nucleotideDF/"

AA_results <- list()

for(g in sleep_genes){

  nuc_file <- file.path(NUC_DIR, paste0(g, ".tsv"))
  if(!file.exists(nuc_file)) next

  aln <- fread(nuc_file) %>% as.data.frame()
  rownames(aln) <- aln[,1]
  aln <- aln[,-1]

  cds_region <- aln[, grepl("CDS", colnames(aln))]

  if(!"homo_sapiens" %in% rownames(cds_region)) next

  human_CDS <- cds_region["homo_sapiens", , drop=FALSE]

  original_codons <- get_codons(as.character(human_CDS[1,]))
  original_aa     <- translate_codons(original_codons)

  snp_positions <- TS_site$Nucleotide_position[
    TS_site$Gene == g
  ]

  for(pos in snp_positions){

    if(!pos %in% colnames(cds_region)) next

    # dominant allele in each phenotype group
    snp_df <- cds_region[, pos, drop=FALSE]

    snp_df$Sleep_Group <-
      metadata$Number_of_sleep_times_per_day[
        match(rownames(snp_df),
              metadata$Species_symbol_name_ensembl)
      ]

    snp_df <- snp_df[!is.na(snp_df$Sleep_Group),]

    major_base <- snp_df %>%
      group_by(Sleep_Group) %>%
      summarise(
        Base = names(which.max(table(.data[[pos]]))),
        .groups="drop"
      )

    for(i in seq_len(nrow(major_base))){

      mutated <- human_CDS
      mutated[1, pos] <- major_base$Base[i]

      mutated_codons <- get_codons(as.character(mutated[1,]))
      mutated_aa     <- translate_codons(mutated_codons)

      diff_idx <- which(original_aa != mutated_aa)

      if(length(diff_idx) == 0) next

      AA_results[[length(AA_results)+1]] <- data.frame(
        Gene = g,
        Nucleotide_position = pos,
        Codon_Position = diff_idx,
        Original_AA = original_aa[diff_idx],
        Mutated_AA  = mutated_aa[diff_idx],
        Mutation_Type = "Non-synonymous"
      )
    }
  }
}

AA_results <- bind_rows(AA_results)

# ----------------------------------------------------------
# 5. Annotate chemical property change
# ----------------------------------------------------------

aa_properties <- data.frame(
  Amino_Acid = c("A","R","N","D","C","E","Q","G","H","I",
                 "L","K","M","F","P","S","T","W","Y","V","*"),
  Chemical_Property = c("Non-polar","Basic","Polar","Acidic","Polar",
                        "Acidic","Polar","Non-polar","Basic","Non-polar",
                        "Non-polar","Basic","Non-polar","Non-polar","Non-polar",
                        "Polar","Polar","Non-polar","Polar","Non-polar","Stop")
)

AA_results <- AA_results %>%
  left_join(aa_properties,
            by=c("Original_AA"="Amino_Acid")) %>%
  rename(Original_Property = Chemical_Property) %>%
  left_join(aa_properties,
            by=c("Mutated_AA"="Amino_Acid")) %>%
  rename(Mutated_Property = Chemical_Property) %>%
  mutate(Property_Change =
           ifelse(Original_Property ==
                    Mutated_Property,
                  "Unchanged","Changed"))

# ----------------------------------------------------------
# 6. Save
# ----------------------------------------------------------

OUT_DIR <- "figures/Result5/"
dir.create(OUT_DIR,
           recursive=TRUE,
           showWarnings=FALSE)

fwrite(
  AA_results,
  file.path(OUT_DIR,
            "AminoAcid_Mutation_DF.tsv"),
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)
