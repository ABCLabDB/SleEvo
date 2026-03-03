############################################################
# Script: Amino Acid Mutation Analysis for NREM-Associated SNPs
# Description:
# This script identifies non-synonymous mutations associated
# with NREM ratio variation and annotates amino acid
# physicochemical property changes.
# Output: Supplementary_table11.tsv
############################################################


############################
# 1. Load Required Libraries
############################

library(ape)
library(data.table)
library(Biostrings)
library(dplyr)
library(tidyr)


############################
# 2. Define File Paths (Modify if Needed)
############################

META_FILE <- "/disk4/bijsy/2.Sleep/2.Result/1.Result1/METADATA_fixedVersion.txt"
FASTA_DIR <- "/disk4/bijsy/Evolution/0.DATA/2.OUT/3.CDS/1.Muscle/"
NUCLEOTIDE_DIR <- "/disk4/bijsy/1.Important/2.Result/1.DATA/4.correlation_Result/nucleotideDF/"
PERM_FILE <- "/disk4/bijsy/2.Sleep/Figure/UPGMA_Bootstrapping/Evolution/Supplementary_table7.tsv"
SNP_FILE <- "/disk4/bijsy/2.Sleep/Figure/UPGMA_Bootstrapping/Evolution/Supplementary_table10.tsv"
OUTPUT_FILE <- "/disk4/bijsy/2.Sleep/Figure/UPGMA_Bootstrapping/Evolution/Supplementary_table11.tsv"


############################
# 3. Load Metadata
############################

metadata <- as.data.frame(fread(META_FILE))
metadata$Species_symbol_name_ensembl <- gsub(" ", "_",
                                            metadata$Species_symbol_name_ensembl)
metadata$Species_symbol_name_ensembl <- tolower(
                                            metadata$Species_symbol_name_ensembl)
meta <- metadata


############################
# 4. Identify Significant Genes
############################

perm_test <- as.data.frame(fread(PERM_FILE))
knee <- perm_test[which(perm_test$adj.P < 0.05),]

all_SNP <- as.data.frame(fread(SNP_FILE))
all_SNP <- all_SNP[which(all_SNP$p.adj < 0.05),]

sleep_gene <- unique(all_SNP$Gene)


############################
# 5. Index CDS Alignment Files
############################

seqMatrix <- list.files(NUCLEOTIDE_DIR, ".tsv")
seqMatrix <- sapply(strsplit(seqMatrix,".",fixed = TRUE),
                    function(x) x[1])
seqMatrix <- seqMatrix[which(seqMatrix %in% sleep_gene)]


############################
# 6. Helper Functions
############################

# Extract codons from human CDS
get_codons <- function(human_CDS) {
  seq <- as.character(unlist(human_CDS[1, ]))
  codons <- sapply(seq(1, length(seq) - 2, by = 3),
                   function(i) paste(seq[i:(i+2)], collapse = ""))
  return(codons)
}

# Translate codons to amino acids
translate_codons <- function(codon_seq) {
  return(translate(DNAStringSet(codon_seq)))
}


############################
# 7. Mutation Identification
############################

aminoAcid_Mutation_DF <- c()

for(i in 1:length(sleep_gene)){
  
  # Load nucleotide alignment
  a <- as.data.frame(
        fread(file.path(NUCLEOTIDE_DIR,
                        paste0(seqMatrix[i], ".tsv")))
      )
  
  rownames(a) <- a$V1
  a <- a[,-1]
  
  cds_region <- a[,which(grepl("CDS",colnames(a)))]
  
  # Extract human reference CDS
  human_CDS <- cds_region[
                  which(rownames(cds_region) == "homo_sapiens"),]
  
  original_codons <- get_codons(human_CDS)
  original_aa <- translate_codons(original_codons)
  
  # Add NREM ratio phenotype
  cds_region <- cbind(
    cds_region,
    meta$Percentage_of_NREM_time_per_day[
      match(rownames(cds_region),
            meta$Species_symbol_name_ensembl)
    ]
  )
  
  colnames(cds_region)[length(cds_region)] <- "NREM_ratio"
  cds_region <- cds_region[
                  -which(is.na(cds_region$NREM_ratio)),]
  
  # Extract SNP positions
  snp_cds_Region <- cds_region[
      ,which(colnames(cds_region) %in%
             c(all_SNP$Column_Index, "NREM_ratio"))
  ]
  
  nrem_ratio <- snp_cds_Region$NREM_ratio
  snp_cds_filtered  <- snp_cds_Region[,-ncol(snp_cds_Region),
                                      drop=FALSE]
  
  # Convert to long format
  long_df <- snp_cds_filtered %>%
    mutate(Species = rownames(snp_cds_filtered)) %>%
    pivot_longer(cols = -Species,
                 names_to = "Position",
                 values_to = "Nucleotide") %>%
    mutate(NREM_ratio = rep(nrem_ratio,
                            ncol(snp_cds_filtered)))
  
  # Calculate mean NREM ratio per nucleotide
  results_df <- long_df %>%
    filter(!is.na(Nucleotide)) %>%
    group_by(Position, Nucleotide) %>%
    summarise(Mean_NREM_Ratio =
                mean(NREM_ratio, na.rm = TRUE),
              .groups = "drop") %>%
    as.data.frame()
  
  # Remove ambiguous bases
  results_df <- results_df[
      !results_df$Nucleotide %in% c("-", "N"),]
  
  # Select max/min phenotype groups
  filtered_results <- results_df %>%
    group_by(Position) %>%
    filter(Mean_NREM_Ratio ==
             max(Mean_NREM_Ratio) |
           Mean_NREM_Ratio ==
             min(Mean_NREM_Ratio)) %>%
    mutate(Category =
             ifelse(Mean_NREM_Ratio ==
                      max(Mean_NREM_Ratio),
                    "long", "short")) %>%
    ungroup() %>%
    as.data.frame()
  
  
  ############################################
  # Mutation effect comparison
  ############################################
  
  analyze_category <- function(category_df){
    
    if(nrow(category_df) == 0) return(NULL)
    
    snp_position <- category_df$Position
    snp_new_base <- category_df$Nucleotide
    
    codon_positions_raw <-
      which(colnames(human_CDS) %in%
            snp_position) / 3
    
    category_df$Codon_Position <-
      floor(codon_positions_raw) +
      ifelse(codon_positions_raw %% 1 > 0, 1, 0)
    
    mutated_CDS <- human_CDS
    mutated_CDS[1, snp_position] <- snp_new_base
    
    mutated_codons <- get_codons(mutated_CDS)
    mutated_aa <- translate_codons(mutated_codons)
    
    mutation_results <- data.frame(
      Original_Codon = original_codons,
      Mutated_Codon = mutated_codons,
      Original_Amino_Acid = as.character(original_aa),
      Mutated_Amino_Acid = as.character(mutated_aa),
      Mutation_Type =
        ifelse(original_aa == mutated_aa,
               "Synonymous",
               "Non-synonymous")
    )
    
    filtered_mutations <- mutation_results[
        mutation_results$Original_Codon !=
        mutation_results$Mutated_Codon, ]
    
    filtered_mutations <- cbind(
        filtered_mutations,
        category_df[
          match(rownames(filtered_mutations),
                category_df$Codon_Position),])
    
    return(filtered_mutations)
  }
  
  short_df <- analyze_category(
                filtered_results[
                  filtered_results$Category=="short",])
  
  long_df <- analyze_category(
                filtered_results[
                  filtered_results$Category=="long",])
  
  merge_AA_mutation <- rbind(short_df, long_df)
  
  if (!is.null(merge_AA_mutation) &&
      nrow(merge_AA_mutation) > 0) {
    
    merge_AA_mutation$Gene <- sleep_gene[i]
    merge_AA_mutation$Full_Length <- length(original_aa)
  }
  
  aminoAcid_Mutation_DF <-
    rbind(aminoAcid_Mutation_DF,
          merge_AA_mutation)
}


############################
# 8. Annotate Amino Acid Properties
############################

get_amino_acid_properties <- function() {
  data.frame(
    Amino_Acid = c("A","R","N","D","C","E","Q","G","H","I",
                   "L","K","M","F","P","S","T","W","Y","V","*"),
    Chemical_Property =
      c("Non-polar","Basic","Polar","Acidic","Polar",
        "Acidic","Polar","Non-polar","Basic","Non-polar",
        "Non-polar","Basic","Non-polar","Non-polar","Non-polar",
        "Polar","Polar","Non-polar","Polar","Non-polar","Stop")
  )
}

analyze_mutation_effects <- function(mutation_df) {
  
  aa_properties <- get_amino_acid_properties()
  
  mutation_df <- left_join(
      mutation_df, aa_properties,
      by = c("Original_Amino_Acid"="Amino_Acid"))
  
  colnames(mutation_df)[
    colnames(mutation_df)=="Chemical_Property"] <-
    "Original_Chemical_Property"
  
  mutation_df <- left_join(
      mutation_df, aa_properties,
      by = c("Mutated_Amino_Acid"="Amino_Acid"))
  
  colnames(mutation_df)[
    colnames(mutation_df)=="Chemical_Property"] <-
    "Mutated_Chemical_Property"
  
  mutation_df <- mutation_df %>%
    mutate(Property_Change =
             ifelse(Original_Chemical_Property ==
                      Mutated_Chemical_Property,
                    "Unchanged","Changed"))
  
  return(mutation_df)
}


############################
# 9. Final Filtering and Export
############################

mutation_results <- analyze_mutation_effects(
                      aminoAcid_Mutation_DF)

mutation_results <- mutation_results[
  which(grepl("Non-synonymous",
              mutation_results$Mutation_Type)),]

fwrite(mutation_results,
       OUTPUT_FILE,
       sep = "\t",
       row.names = TRUE,
       col.names = TRUE,
       quote = FALSE)

############################################################
# End of Script
############################################################
