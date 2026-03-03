############################################################
# Script: Identification of Total Sleep Time–related SNPs
# Description:
# This script identifies nucleotide positions significantly
# associated with total sleep time across species using
# t-tests (2 groups) and ANOVA (≥3 groups).
# Significant SNP sites are exported as Supplementary Table 4.
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
# 2. Define File Paths
# (Modify these paths as needed)
############################

BASE_DIR <- "/data/"
NUCLEOTIDE_DIR <- paste0(BASE_DIR, "Circadian_gene_Nucleotide_Matrix/")
META_FILE <- paste0(BASE_DIR, "/Result1/species_sleep_metadata.txt")
PERM_FILE <- paste0(BASE_DIR, "/Result2/Total_sleep_time_Anova_Result.tsv")
GENE_LIBRARY_FILE <- paste0(BASE_DIR, "Result1/circadian_Gene_list.tsv")
OUTPUT_FILE <- paste0(BASE_DIR, "/Result2/Total_sleep_time_SNP.tsv")


############################
# 3. Load Species Metadata
############################

species_Metadata <- as.data.frame(fread(META_FILE))

# Replace spaces with underscores for consistency
species_Metadata$Species_symbol_name_ensembl <- 
  gsub(" ","_",species_Metadata$Species_symbol_name_ensembl)


############################
# 4. Load Significant Genes from Permutation Test
############################

perm_DF <- as.data.frame(fread(PERM_FILE))

# Select genes passing adjusted P-value threshold
perm_test <- perm_DF[which(perm_DF$adj.P < 0.05),]
rownames(perm_test) <- NULL


############################
# 5. Match Significant Genes to Alignment Files
############################

seqMatrix <- list.files(NUCLEOTIDE_DIR, ".tsv")
gene <- sapply(strsplit(seqMatrix,".", fixed = TRUE), function(x) x[1])

# Keep only files corresponding to significant genes
seqMatrix <- seqMatrix[match(perm_test$Gene, gene)]
gene <- perm_test$Gene


############################
# 6. Identify Total Sleep Time–Associated SNP Sites
############################

TS_site <- c()

for(i in 1:length(gene)){
  
  # Load alignment matrix for each gene
  aln_df <- as.data.frame(fread(file.path(NUCLEOTIDE_DIR, seqMatrix[i])))
  rownames(aln_df) <- aln_df[,1]
  aln_df <- aln_df[,-1]
  
  # Identify column types
  V <- which(grepl("V",colnames(aln_df)))
  CDS <- which(grepl("CDS",colnames(aln_df)))
  IDX <- which(grepl("Total",colnames(aln_df)))
  
  # Keep SNP + CDS + phenotype columns
  aln_df <- aln_df[,c(V,CDS,IDX)]
  sleep_df <- aln_df
  
  # Remove rows with missing sleep phenotype
  if(length(which(is.na(sleep_df$Total_sleep))) > 0){
    sleep_df <- sleep_df[-which(is.na(sleep_df$Total_sleep)),]
  }
  
  # Count number of unique states per column
  unique_counts <- sapply(sleep_df, function(x) length(unique(x)))
  
  one <- which(unname(unique_counts)==1)
  two <- which(unname(unique_counts)==2)
  more <- which(unname(unique_counts) %in% c(3,4,5,6))
  
  
  ############################################
  # 6A. Two-group comparison (t-test)
  ############################################
  
  two_MA <- sleep_df[,two]
  two_MA <- cbind(two_MA,sleep_df$Total_sleep)
  
  num_columns <- ncol(two_MA)
  results_df <- data.frame(Column_Index = integer(0), p_value = numeric(0))
  
  for (j in 1:(num_columns-1)) {
    
    a <- two_MA[,c(j,tail(j:ncol(two_MA), n = 1))]
    colnames(a)[1] <- "Nucleotide"
    colnames(a)[2] <- "Sleep_Time"
    
    type1 <- names(table(a[,1]))[1]
    type2 <- names(table(a[,1]))[2]
    
    group1_df <- a$Sleep_Time[which(grepl(type1,a$Nucleotide))]
    group2_df <- a$Sleep_Time[which(grepl(type2,a$Nucleotide))]
    
    n1 <- length(group1_df)
    n2 <- length(group2_df)
    
    # Perform t-test only if both groups have ≥2 samples
    if (n1 >= 2 && n2 >= 2) {
      t_test_result <- t.test(group1_df, group2_df, alternative = "two.sided")
      results_df <- rbind(results_df, 
                          data.frame(Column_Index = colnames(two_MA)[j], 
                                     p_value = t_test_result$p.value))
    }
  }
  
  if(nrow(results_df) == 0 ){
    results_df <- data.frame()
  }else{
    results_df$Gene <- gene[i]
    results_df$Phenotype <- "Total_sleep"
    results_df$Column_Index <- paste0("V",results_df$Column_Index)
    results_df$Length <- ncol(sleep_df) - 1
    results_df$p.adj <- p.adjust(results_df$p_value,method = "bonferroni")
  }
  
  
  ############################################
  # 6B. Multi-group comparison (ANOVA)
  ############################################
  
  anova_MA <- sleep_df[,c(more,tail(seq_along(sleep_df),1))]
  
  x <- c()
  y <- c()
  
  for(k in 1:ncol(anova_MA[,-c(tail(seq_along(anova_MA),1))])){
    
    anova <- aov(Total_sleep ~ as.factor(anova_MA[,k]), data = anova_MA)
    summary <- summary(anova)
    
    x <- append(x,summary[[1]]$`Pr(>F)`[1])
    y <- append(y,summary[[1]]$`F value`[1] )
  }
  
  three_group_anova <- data.frame(
    Column_Index = colnames(anova_MA[,-c(tail(seq_along(anova_MA),1))]),
    p_value = x,
    Gene= gene[i],
    Phenotype = "Total_sleep",
    Length = (ncol(sleep_df) - 1),
    F_value = y,
    p.adj = p.adjust(x,method = "bonferroni")
  )
  
  
  ############################################
  # 6C. Merge Significant Sites
  ############################################
  
  if(nrow(results_df) == 0 ){
    merge_site <- three_group_anova
  }else{
    results_df$F_value <- NA
    merge_site <- rbind(results_df,three_group_anova)
  }
  
  rownames(merge_site)=NULL
  TS_site <- rbind(TS_site, merge_site)
}


############################
# 7. Filter Significant SNP Sites
############################

per_Snp <- TS_site
per_Snp <- per_Snp[which(per_Snp$p.adj < 0.05),]

# Clean column names
per_Snp$Column_Index <- gsub("VCDS","CDS",per_Snp$Column_Index)
per_Snp <- per_Snp[-which(grepl("VV",per_Snp$Column_Index)),]
per_Snp <- per_Snp[-which(grepl("V",per_Snp$Column_Index)),]


############################
# 8. Add Chromosome Annotation
############################

gene_Library <- as.data.frame(fread(GENE_LIBRARY_FILE))
gene_Library <- gene_Library[,c(1:4)]

per_Snp$Chr <- gene_Library$Chr[
  match(per_Snp$Gene, gene_Library$Gene_symbol)
]

per_Snp <- per_Snp[,c(3,8,1,5,2,6)]


############################
# 9. Export Results
############################

fwrite(per_Snp, OUTPUT_FILE, 
       sep = "\t", 
       col.names = TRUE, 
       quote = FALSE)

############################################################
# End of Script
############################################################
