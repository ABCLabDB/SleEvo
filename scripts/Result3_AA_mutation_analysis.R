############################################################
# Script: Amino Acid Mutation Analysis and Visualization
# Description:
#   This script performs the following analyses:
#
#   1) Identifies non-synonymous mutations associated with
#      NREM ratio variation from significant SNP sites.
#   2) Compares original and mutated codons relative to the
#      human reference CDS.
#   3) Annotates amino acid physicochemical property changes.
#   4) Exports mutation table (Supplementary_table11.tsv).
#   5) Generates protein domain mutation visualization
#      for selected genes (Figure4F).
#
# Input:
#   - Species sleep metadata
#   - NREM ratio enrichment results
#   - Significant SNP table
#   - CDS nucleotide alignment matrices
#
# Output:
#   - Supplementary_table11.tsv
#   - Protein mutation domain plot (Figure4F.pdf)
#
# Notes:
#   - Only non-synonymous mutations are retained.
#   - Human CDS is used as the reference sequence.
#   - Physicochemical property changes are annotated.
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

META_FILE <- "/data/Result1/species_sleep_metadata.txt"
FASTA_DIR <- "/data/Fasta/"
NUCLEOTIDE_DIR <- "/data/Circadian_gene_Nucelotide_Matrix/"
PERM_FILE <- "/data/Result3/NREM_ratio_Anova_Result.tsv"
SNP_FILE <- "/data/Result3/NREM_ratio_SNP.tsv"
OUTPUT_FILE <- "/data/Result3/NREM_ratio_AminoAcid_mutation.tsv"


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
# 10. Protein Domain Mutation Plot (Figure4F)
############################################################

library(ggrepel)
library(ggplot2)
library(dplyr)

# Load generated mutation table
mutation_results_plot <- as.data.frame(
  fread(OUTPUT_FILE)
)

# Remove row index column if present
mutation_results_plot <- mutation_results_plot[,-c(1)]

# Clean Position format
mutation_results_plot$Position <-
  sapply(strsplit(mutation_results_plot$Position,"_"),
         function(x) x[2])

mutation_results_plot$Position <-
  paste0("Chr",
         mutation_results_plot$Chr,
         ":",
         mutation_results_plot$Position)

############################################################
# Filter Target Gene (ARNTL2)
############################################################

protein_Gene <- mutation_results_plot %>%
  filter(Gene == "ARNTL2") %>%
  mutate(
    mutation_label =
      paste0(Original_Amino_Acid,
             Codon_Position,
             Mutated_Amino_Acid,
             " (",Position,")"),
    
    y = ifelse(Category == "long", 1.1, 1),
    
    yend = ifelse(Category == "long",
                  runif(n(), min = 1.3, max = 1.38),
                  runif(n(), min = 0.65, max = 0.75)),
    
    label_y = yend
  )

full_length <- unique(protein_Gene$Full_Length)

############################################################
# Force specific mutation label to be visible
############################################################

protein_Gene$mutation_label[
  protein_Gene$Position != "Chr12:27401345"
] <- protein_Gene$mutation_label[
  protein_Gene$Position != "Chr12:27401345"
]

############################################################
# Domain annotation
############################################################

domain_df <- data.frame(
  name  = c("bHLH","PAS","PAS","PAC"),
  start = c(107,178,357,432),
  end   = c(160,250,427,475),
  fill  = c("#fbf2c4","#74a892","#74a892","#003f5c")
)

############################################################
# Plot
############################################################

p1 <- ggplot() +
  
  geom_rect(aes(xmin = 1,
                xmax = full_length,
                ymin = 1,
                ymax = 1.1),
            fill = "#ffffff",
            color = "black") +
  
  geom_rect(data = domain_df,
            aes(xmin = start,
                xmax = end,
                ymin = 1,
                ymax = 1.1,
                fill = name),
            color = "black",
            alpha = 0.8) +
  
  geom_segment(data = protein_Gene,
               aes(x = Codon_Position,
                   xend = Codon_Position,
                   y = y,
                   yend = yend),
               color = "black") +
  
  geom_point(data = protein_Gene,
             aes(x = Codon_Position,
                 y = yend,
                 color = Category),
             size = 3) +
  
  geom_text_repel(
    data = protein_Gene,
    aes(x = Codon_Position,
        y = yend,
        label = mutation_label),
    direction = "y",
    size = 3.5,
    fontface = "bold",
    box.padding = 0.2,
    point.padding = 0.2,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  scale_fill_manual(values = setNames(domain_df$fill,
                                      domain_df$name)) +
  
  scale_x_continuous(
    breaks = seq(0, full_length, by = 100),
    limits = c(0, full_length + 10)
  ) +
  
  coord_cartesian(ylim = c(0.6, 1.5)) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5,
                              face = "bold")
  ) +
  
  labs(
    x = "",
    title = "ARNTL2 Gene Mutation"
  )

ggsave(
  "/disk4/bijsy/2.Sleep/Figure/UPGMA_Bootstrapping/Evolution/Figure4F.pdf",
  plot = p1,
  width = 7,
  height = 3,
  units = "in",
  dpi = 300
)

############################################################
# End of Extended Script
############################################################
