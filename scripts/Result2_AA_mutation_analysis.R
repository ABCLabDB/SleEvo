############################################################
# Result2_AA_mutation_analysis.R
#
# Amino-acid mutation analysis for sleep-associated SNPs
# across circadian genes.
#
# Output:
#   data/Result2/AminoAcid_Mutation_DF.tsv
############################################################


############################
# 1. Load Required Libraries
############################

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
                       "data/Result2/Total_sleep_time_Anova_Result.tsv")

SNP_FILE <- file.path(PROJECT_DIR,
                      "data/Result2/Total_Sleep_Time.tsv")

SEQ_DIR <- file.path(PROJECT_DIR,
                     "data/Result2/NucleotideMatrix")

OUT_FILE <- file.path(PROJECT_DIR,
                      "data/Result2/Total_sleep_time_AA_Mutation_Data.tsv")


############################
# 3. Load Metadata
############################

metadata <- as.data.frame(fread(META_FILE))

# Standardize species naming
metadata$Species_symbol_name_ensembl <- gsub(" ", "_", metadata$Species_symbol_name_ensembl)
metadata$Species_symbol_name_ensembl <- tolower(metadata$Species_symbol_name_ensembl)

meta <- metadata

############################
# 4. Load SNP List
############################

all_SNP <- as.data.frame(fread(SNP_FILE))
sleep_gene <- unique(all_SNP$Gene)

############################
# 5. Identify Corresponding CDS Alignment Files
############################

seqMatrix <- list.files(NUCLEOTIDE_DIR, ".tsv")
seqMatrix <- sapply(strsplit(seqMatrix, ".", fixed = TRUE), function(x) x[1])
seqMatrix <- seqMatrix[which(seqMatrix %in% sleep_gene)]

############################
# 6. Define Helper Functions
############################

# Extract codons from CDS sequence
get_codons <- function(human_CDS) {
  seq <- as.character(unlist(human_CDS[1, ]))
  codons <- sapply(seq(1, length(seq) - 2, by = 3),
                   function(i) paste(seq[i:(i+2)], collapse = ""))
  return(codons)
}

# Translate codon sequence into amino acids
translate_codons <- function(codon_seq) {
  return(translate(DNAStringSet(codon_seq)))
}

############################
# 7. Identify Amino Acid Mutations
############################

aminoAcid_Mutation_DF <- c()

for(i in 1:length(sleep_gene)){
  
  # Load nucleotide alignment matrix
  a <- as.data.frame(
    fread(paste0(NUCLEOTIDE_DIR, seqMatrix[i], ".tsv"))
  )
  
  rownames(a) <- a$V1
  a <- a[,-1]
  
  # Extract CDS region
  cds_region <- a[,which(grepl("CDS",colnames(a)))]
  
  # Extract human CDS
  human_CDS <- cds_region[which(rownames(cds_region) == "homo_sapiens"),]
  
  # Get original codons and amino acids
  original_codons <- get_codons(human_CDS)
  original_aa <- translate_codons(original_codons)
  
  # Attach sleep phenotype
  cds_region <- cbind(
    cds_region,
    meta$Total_sleep_time_per_day[
      match(rownames(cds_region),
            meta$Species_symbol_name_ensembl)
    ]
  )
  colnames(cds_region)[length(cds_region)] <- "Total_Sleep_Time"
  
  # Extract SNP positions for this gene
  snp_cds_Region <- cds_region[,which(colnames(cds_region) %in% c(all_SNP$Column_Index,"Total_Sleep_Time"))]
  
  ############################################
  # Summarize nucleotide effects
  ############################################
  
  nucleotide_cols <- colnames(snp_cds_Region)[grepl("CDS_", colnames(snp_cds_Region))]
  results_list <- list()
  
  for (nucleotide_Position in nucleotide_cols) {
    
    tmp <- snp_cds_Region[, c(nucleotide_Position, "Total_Sleep_Time")]
    
    summary_df <- tmp %>%
      group_by(!!sym(nucleotide_Position)) %>%
      summarise(
        Mean_TST = mean(Total_Sleep_Time, na.rm = TRUE),
        Sample_Count = n()
      ) %>%
      rename(Nucleotide = !!sym(nucleotide_Position)) %>%
      mutate(Position = nucleotide_Position)
    
    results_list[[nucleotide_Position]] <- summary_df
  }
  
  results_df <- as.data.frame(bind_rows(results_list))
  
  # Remove ambiguous or singleton categories
  results_df <- results_df[results_df$Nucleotide != "-", ]
  results_df <- results_df[results_df$Nucleotide != "N", ]
  results_df <- results_df[results_df$Sample_Count > 1, ]
  
  ############################################
  # Identify extreme sleep categories
  ############################################
  
  filtered_results <- results_df %>%
    group_by(Position) %>%
    filter(Mean_TST == max(Mean_TST) |
           Mean_TST == min(Mean_TST)) %>%
    mutate(Category = ifelse(
      Mean_TST == max(Mean_TST), "long", "short"
    )) %>%
    ungroup() %>%
    as.data.frame()
  
  
  ############################################
  # Compare amino acid mutations (short vs long)
  ############################################
  
  for(category_type in c("short", "long")){
    
    snp_df <- filtered_results[
      filtered_results$Category %in% category_type, ]
    
    snp_df <- snp_df[
      !duplicated(snp_df[, c("Position", "Mean_TST")]), ]
    
    snp_position <- snp_df$Position
    snp_new_base <- snp_df$Nucleotide
    
    codon_positions_raw <-
      which(colnames(human_CDS) %in% snp_position) / 3
    
    snp_df$Codon_Position <-
      floor(codon_positions_raw) +
      ifelse(codon_positions_raw %% 1 > 0, 1, 0)
    
    snp_df$Codon_In_Position <-
      ifelse(codon_positions_raw %% 1 == 0, 3,
             ifelse(codon_positions_raw %% 1 <= 0.34, 1, 2))
    
    mutated_CDS <- human_CDS
    mutated_CDS[1, snp_position] <- snp_new_base
    
    mutated_codons <- get_codons(mutated_CDS)
    mutated_aa <- translate_codons(mutated_codons)
    
    mutation_results <- data.frame(
      Original_Codon = original_codons,
      Mutated_Codon = mutated_codons,
      Original_Amino_Acid = as.character(original_aa),
      Mutated_Amino_Acid = as.character(mutated_aa),
      Mutation_Type = ifelse(
        original_aa == mutated_aa,
        "Synonymous", "Non-synonymous"
      )
    )
    
    filtered_mutations <-
      mutation_results[
        mutation_results$Original_Codon !=
          mutation_results$Mutated_Codon, ]
    
    filtered_mutations <-
      cbind(filtered_mutations,
            snp_df[match(
              rownames(filtered_mutations),
              snp_df$Codon_Position),])
    
    if(nrow(filtered_mutations) > 0){
      filtered_mutations$Gene <- sleep_gene[i]
      filtered_mutations$Full_Length <- length(original_aa)
      aminoAcid_Mutation_DF <-
        rbind(aminoAcid_Mutation_DF,
              filtered_mutations)
    }
  }
}


############################
# 8. Amino Acid Chemical Property Analysis
############################

get_amino_acid_properties <- function() {
  data.frame(
    Amino_Acid = c("A","R","N","D","C","E","Q","G","H","I",
                   "L","K","M","F","P","S","T","W","Y","V","*"),
    Chemical_Property = c("Non-polar","Basic","Polar","Acidic",
                          "Polar","Acidic","Polar","Non-polar",
                          "Basic","Non-polar","Non-polar",
                          "Basic","Non-polar","Non-polar",
                          "Non-polar","Polar","Polar",
                          "Non-polar","Polar","Non-polar","Stop")
  )
}

analyze_mutation_effects <- function(mutation_df) {
  
  aa_properties <- get_amino_acid_properties()
  
  mutation_df <- mutation_df %>%
    left_join(aa_properties,
              by = c("Original_Amino_Acid"="Amino_Acid")) %>%
    rename(Original_Chemical_Property = Chemical_Property)
  
  mutation_df <- mutation_df %>%
    left_join(aa_properties,
              by = c("Mutated_Amino_Acid"="Amino_Acid")) %>%
    rename(Mutated_Chemical_Property = Chemical_Property)
  
  mutation_df <- mutation_df %>%
    mutate(Property_Change =
             ifelse(Original_Chemical_Property ==
                      Mutated_Chemical_Property,
                    "Unchanged","Changed"))
  
  return(mutation_df)
}

############################
# 9. Filter Non-synonymous Mutations
############################

mutation_results <- analyze_mutation_effects(aminoAcid_Mutation_DF)
mutation_results <- mutation_results[
  grepl("Non-synonymous",
        mutation_results$Mutation_Type), ]


############################
# 10. Add Chromosome Annotation
############################

gene_Library <- as.data.frame(fread(GENE_LIBRARY_FILE))
gene_Library <- gene_Library[,c(1:4)]

mutation_results$Chr <-
  gene_Library$Chr[
    match(mutation_results$Gene,
          gene_Library$Gene_symbol)
  ]

mutation_results$Position <-
  sapply(strsplit(mutation_results$Position,
                  "CDS_", fixed = TRUE),
         function(x) x[2])


############################
# 11. Export Final Table
############################

mutation_results <- mutation_results[
  ,c(13,18,9,6,1:5,7,10,11,15:17)]

fwrite(mutation_results,
       OUTPUT_FILE,
       sep = "\t",
       col.names = TRUE,
       quote = FALSE)

############################################################
# End of Script
############################################################
