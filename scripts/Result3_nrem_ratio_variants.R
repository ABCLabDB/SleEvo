############################################################
# Script: SNP Association Analysis for NREM Ratio
# Description:
# This script identifies nucleotide positions associated
# with NREM ratio across species using non-parametric and
# parametric tests depending on sample size.
# Output: NREM_Ratio_SNP.tsv
############################################################


############################
# 1. Load Required Libraries
############################

library(ape)
library(data.table)
library(ggplot2)
library(ggtree)
library(patchwork)
library(phyloseq)
library(Biostrings)
library(ade4)
library(adegenet)
library(pheatmap)
library(dplyr)
library(stringr)


############################
# 2. Define File Paths (Modify as Needed)
############################

BASE_DIR <- "/data/"
NUCLEOTIDE_DIR <- "/data/Circadian_gene_Nucelotide_Matrix/"
META_FILE <- "/data/Result1/species_sleep_metadata.txt"
PERM_FILE <- paste0(BASE_DIR, "Result3/NREM_ratio_Anova_Result.tsv")
GENE_LIBRARY_FILE <- paste0(BASE_DIR, "Result1/circadian_Gene_list.tsv")
OUTPUT_FILE <- paste0(BASE_DIR, "Result3/NREM_Ratio_SNP.tsv") #Supplementary Table 10


############################
# 3. Load Species Metadata
############################

species_Metadata <- as.data.frame(fread(META_FILE))
species_Metadata$Species_symbol_name_ensembl <- 
  gsub(" ", "_", species_Metadata$Species_symbol_name_ensembl)


############################
# 4. Load Significant Genes
############################

perm_test <- as.data.frame(fread(PERM_FILE))
perm_test <- perm_test[which(perm_test$adj.P < 0.05),]
rownames(perm_test) <- NULL


############################
# 5. Match Alignment Files
############################

seqMatrix <- list.files(NUCLEOTIDE_DIR, ".tsv")

gene_names <- sapply(strsplit(seqMatrix, ".", fixed = TRUE),
                     function(x) x[1])

seqMatrix <- seqMatrix[match(perm_test$Gene, gene_names)]
gene <- perm_test$Gene


############################
# 6. SNP Association Analysis
############################

TS_site <- c()

for(i in 1:length(gene)){
  
  # Load alignment file
  aln_df <- as.data.frame(
    fread(file.path(NUCLEOTIDE_DIR, seqMatrix[i]))
  )
  
  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]
  
  # Keep only variant and CDS columns
  remove_IDX <- which(!grepl("V|CDS", colnames(aln_df)))
  aln_df1 <- aln_df[,-remove_IDX]
  
  rownames(aln_df1) <- str_to_title(rownames(aln_df1))
  
  # Add NREM ratio
  aln_df1$NREM_Ratio <- 
    species_Metadata$Percentage_of_NREM_time_per_day[
      match(rownames(aln_df1),
            species_Metadata$Species_symbol_name_ensembl)
    ]
  
  sleep_df <- aln_df1
  
  # Remove missing phenotype
  if(length(which(is.na(sleep_df$NREM_Ratio))) > 0){
    sleep_df <- sleep_df[
      -which(is.na(sleep_df$NREM_Ratio)),]
  }
  
  # Count unique states per column
  unique_counts <- sapply(sleep_df,
                          function(x) length(unique(x)))
  
  one  <- which(unname(unique_counts)==1)
  two  <- which(unname(unique_counts)==2)
  more <- which(unname(unique_counts) %in% c(3,4,5,6))
  
  
  ############################################
  # 6A. Two-group comparison
  ############################################
  
  two_MA <- sleep_df[,two]
  two_MA <- cbind(two_MA, sleep_df$NREM_Ratio)
  
  num_columns <- ncol(two_MA)
  results_df <- data.frame(Column_Index = integer(0),
                           p_value = numeric(0))
  
  for (j in 1:(num_columns-1)) {
    
    a <- two_MA[,c(j,tail(j:ncol(two_MA), n = 1))]
    colnames(a)[1] <- "Nucleotide"
    colnames(a)[2] <- "NREM_Ratio"
    
    type1 <- names(table(a[,1]))[1]
    type2 <- names(table(a[,1]))[2]
    
    group1_df <- a$NREM_Ratio[
      which(grepl(type1,a$Nucleotide))]
    group2_df <- a$NREM_Ratio[
      which(grepl(type2,a$Nucleotide))]
    
    n1 <- length(group1_df)
    n2 <- length(group2_df)
    
    if (n1 < 3 | n2 < 3){
      
      results_df <- rbind(results_df,
        data.frame(Column_Index = colnames(two_MA)[j],
                   p_value = NA,
                   Method = "Too small sample size",
                   Real_Position = two[j]))
      
    } else if (n1 >= 30 & n2 >= 30){
      
      t_test_result <- t.test(group1_df, group2_df)
      
      results_df <- rbind(results_df,
        data.frame(Column_Index = colnames(two_MA)[j],
                   p_value = t_test_result$p.value,
                   Method = "T-test",
                   Real_Position = two[j]))
      
    } else {
      
      res <- wilcox.test(group1_df, group2_df)
      
      results_df <- rbind(results_df,
        data.frame(Column_Index = colnames(two_MA)[j],
                   p_value = res$p.value,
                   Method = "Wilcoxon",
                   Real_Position = two[j]))
    }
  }
  
  if(nrow(results_df) > 0){
    results_df$Gene <- gene[i]
    results_df$Phenotype <- "NREM_Ratio"
    results_df$Length <- ncol(sleep_df) - 1
    results_df$p.adj <- p.adjust(results_df$p_value,
                                 method = "BH")
  }
  
  
  ############################################
  # 6B. Multi-group comparison
  ############################################
  
  anova_MA <- sleep_df[,c(more,
                          tail(seq_along(sleep_df),1))]
  
  three_group_anova <- data.frame(
    Column_Index = integer(0),
    p_value = numeric(0))
  
  for(k in 1:(ncol(anova_MA)-1)){
    
    base_vec <- anova_MA[, k]
    tbl <- as.data.frame(table(base_vec))
    colnames(tbl) <- c("base", "freq")
    
    valid_group_count <- sum(tbl$freq >= 3)
    
    if (valid_group_count < 2) {
      
      three_group_anova <- rbind(three_group_anova,
        data.frame(Column_Index =
                     colnames(anova_MA)[k],
                   p_value = NA,
                   Method = "Too small sample size",
                   Real_Position = more[k]))
      
    } else {
      
      kru <- kruskal.test(
        NREM_Ratio ~ as.factor(anova_MA[,k]),
        data = anova_MA)
      
      three_group_anova <- rbind(three_group_anova,
        data.frame(Column_Index =
                     colnames(anova_MA)[k],
                   p_value = kru$p.value,
                   Method = "Kruskal-Wallis",
                   Real_Position = more[k]))
    }
  }
  
  if(nrow(three_group_anova) > 0){
    three_group_anova$Gene <- gene[i]
    three_group_anova$Phenotype <- "NREM_Ratio"
    three_group_anova$Length <- ncol(sleep_df) - 1
    three_group_anova$p.adj <-
      p.adjust(three_group_anova$p_value,
               method = "BH")
  }
  
  ############################################
  # 6C. Merge results
  ############################################
  
  if(nrow(results_df) == 0){
    merge_site <- three_group_anova
  } else{
    merge_site <- rbind(results_df,
                        three_group_anova)
  }
  
  rownames(merge_site) <- NULL
  TS_site <- rbind(TS_site, merge_site)
}


############################
# 7. Final Processing
############################

snp_Data <- TS_site[
  which(grepl("CDS_", TS_site$Column_Index)),]

rownames(snp_Data) <- NULL

gene_Library <- as.data.frame(fread(GENE_LIBRARY_FILE))
gene_Library <- gene_Library[,c(1:4)]

snp_Data$Chr <- gene_Library$Chr[
  match(snp_Data$Gene,
        gene_Library$Gene_symbol)
]

snp_Data <- snp_Data[,c(5,10,1,7,2,8)]


############################
# 8. Export Results
############################

fwrite(snp_Data,
       OUTPUT_FILE,
       sep = "\t",
       row.names = FALSE,
       col.names = TRUE,
       quote = FALSE)

############################################################
# End of Script
############################################################
