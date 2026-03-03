############################################################
# Script: LD-like Co-evolution Heatmap for NREM-associated SNPs
# Description:
# This script visualizes LD-like correlation (mean r²)
# among non-synonymous SNPs associated with NREM ratio.
############################################################

############################
# 1. Load Required Libraries
############################

library(data.table)
library(dplyr)
library(stringr)
library(ggcorrplot)

############################
# 2. Define File Paths (Edit as Needed)
############################

NREM_GENE_FILE <- "/data/Result3/NREM_ratio_Anova_Result.tsv"
SNP_FILE       <- "/data/Result3/NREM_ratio_SNP.tsv"
MUT_FILE       <- "/data/Result3/NREM_ratio_AminoAcid_mutation.tsv"
META_FILE      <- "/data/Result1/species_sleep_metadata.txt"
CDS_DIR        <- "/data/Circadian_gene_Nucleotide_Matrix/"
CHR_FILE       <- "/data/Result1/circadian_Gene_list.tsv"

############################
# 3. Load Significant Genes
############################

nrem_gene <- fread(NREM_GENE_FILE) %>%
  filter(adj.P < 0.05)

all_SNP <- fread(SNP_FILE) %>%
  filter(Gene %in% nrem_gene$Gene,
         grepl("CDS", Column_Index),
         p.adj < 0.05)

sleep_gene <- unique(all_SNP$Gene)

############################
# 4. Load Metadata
############################

metadata <- fread(META_FILE)
metadata$Species_symbol_name_ensembl <- tolower(
  gsub(" ", "_", metadata$Species_symbol_name_ensembl)
)

############################
# 5. Example Gene Selection
############################

target_gene <- "ARNTL2"

############################
# 6. Load Alignment
############################

aln_df <- fread(file.path(CDS_DIR, paste0(target_gene, ".tsv")))
rownames(aln_df) <- aln_df[[1]]
aln_df <- aln_df[,-1]

############################
# 7. Load Mutation Data
############################

mutation_results <- fread(MUT_FILE)

mutation_results <- mutation_results %>%
  filter(Gene == target_gene,
         Mutation_Type == "Non-synonymous",
         Property_Change == "Changed")

############################
# 8. Extract Relevant SNP Positions
############################

positions <- unique(mutation_results$Position)

cor_DF <- aln_df[, c(positions, "Percentage_of_NREM_time_per_day"), with = FALSE]

############################
# 9. Binary Encoding (Reference = Human)
############################

ref_seq <- cor_DF["Human", 1:length(positions)]

binary_df <- as.data.frame(
  t(apply(cor_DF[, 1:length(positions)], 1,
          function(x) as.numeric(x != ref_seq)))
)

colnames(binary_df) <- positions

############################
# 10. Compute r² Matrix
############################

cor_matrix <- cor(binary_df,
                  use = "pairwise.complete.obs")^2

mean_r2 <- mean(cor_matrix[upper.tri(cor_matrix)],
                na.rm = TRUE)

############################
# 11. Plot Heatmap
############################

ggcorrplot(
  cor_matrix,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  colors = c("#d9f0a3", "#1b7837"),
  title = sprintf("LD-like co-evolution heatmap (mean r² = %.2f)",
                  mean_r2)
)

############################################################
# End of Script
############################################################
