############################################################
# Directory settings
############################################################

DATA.DIR <- "/data/Fasta/"

############################################################
# Load required libraries
############################################################

library(ape)
library(hierfstat)
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
library(mclust)
library(factoextra)
library(stats)
library("ggpubr")
library(gridExtra)
library(cowplot)
library(dispRity)
library(fossil)
library(DescTools)


############################################################
# Load sequence alignment files
############################################################

setwd(DATA.DIR)
muscle <- list.files(DATA.DIR,".fasta")

# Extract gene names from FASTA file names
gene <- sapply(strsplit(muscle,"_"),function(x) x[1])


############################################################
# Function: making.sum.cluster
# Purpose:
# Resolve tie clusters and aggregate cluster counts
# for sleep phenotype groups.
############################################################

making.sum.cluster <- function(result_table.df,j){
  
  
  tie.idx <- 0
  
  ### Handle tie cases between clusters
  ## When the same counts appear in both phenotype rows
  if(sum(result_table.df[1,] == result_table.df[2,]) != 0){
    
    tie.df <- c()
    
    same.idx <- which(result_table.df[1,] == result_table.df[2,])
    
    same.cl <- names(result_table.df)[same.idx]
    
    tie.df <- result_table.df %>% dplyr::select(all_of(same.cl))
    
    
    if(sum(colSums(tie.df) != 0) != 0){
      
      idx <- which(colSums(tie.df) != 0)  %>% unlist()
      
      idx.names <- names(tie.df)[idx]
      
      
      tie.idx <- c()
      
      for(p in 1:length(idx)){
        print(j)
        
        tie.idx[p] <- unique((tie.df[,idx[p]])/2)
        result_table.df <- result_table.df %>% dplyr::select(-all_of(idx.names[p]))
        
      }
      
      tie.idx <- sum(tie.idx)
      
      cat(paste("tie idx is...", tie.idx), "\n")
      
    }
    
    
  }
  
  
  ############################################################
  # Cluster aggregation
  ############################################################
  
  cluster.sum <- colSums(result_table.df)
  
  ## Calculate phenotype ratios per cluster
  bottom.res <- result_table.df[rownames(result_table.df) == "less_sleep",]/cluster.sum
  top.res <- result_table.df[rownames(result_table.df) == "more_sleep",]/cluster.sum
  
  
  ### Replace NaN values with zero
  if(sum(is.na(bottom.res)) != 0){
    bottom.res[which(is.na(bottom.res))] <- 0
  } 
  
  if(sum(is.na(top.res)) != 0){
    top.res[which(is.na(top.res))] <- 0
  } 
  
  
  ### Combine ratio matrices
  result.oi <- rbind(bottom.res, top.res)
  
  ### Determine dominant phenotype per cluster
  max.idx <- apply(result.oi, 2,function(x){
    rownames(result.oi)[which.max(x)]
  })
  
  
  less.sleep.idx <- names(max.idx[max.idx == "less_sleep"])
  more.sleep.idx <- names(max.idx[max.idx == "more_sleep"])
  
  less.sleep.df <- result_table.df %>% dplyr::select(all_of(less.sleep.idx))
  more.sleep.df <- result_table.df %>% dplyr::select(all_of(more.sleep.idx))
  
  
  ### Final cluster aggregation
  cluster.sum.df  <- data.frame(Cluster.less = rowSums(less.sleep.df),
                                Cluster.more = rowSums(more.sleep.df))
  
  
  if(tie.idx != 0){
    cluster.sum.df <- cluster.sum.df + tie.idx
  }
  
  
  return(cluster.sum.df)
  
}


############################################################
# Function: calculate.purity
# Purpose:
# Compute purity index of phenotype assignment
# across clusters.
############################################################

calculate.purity <- function(cluster.sum.df){
  
  ### Compute cluster totals
  new.cluster.sum <- colSums(cluster.sum.df)
  
  ### Identify dominant phenotype per cluster
  new.test.df <- apply(cluster.sum.df, 2, max)
  
  new.result.oi <- new.test.df/new.cluster.sum
  
  if(sum(is.na(new.result.oi)) > 0){
    
    new.result.oi[is.na(new.result.oi)] <- 0
    
  }
  
  
  ### Calculate purity score
  purity.score <- sum(new.result.oi)/ncol(cluster.sum.df)
  
  return(purity.score)
  
}


############################################################
# Load metadata
############################################################

metadata <- as.data.frame(fread("/data/Result1/species_sleep_metadata.txt"))

# Standardize species naming
metadata$Species_symbol_name_ensembl <- gsub(" ", "_",metadata$Species_symbol_name_ensembl)

phenotype <- colnames(metadata)[c(8,11:13)]

meta <- metadata[order(metadata$Sleep_timing_per_day),
                 c(1,2,3,which(colnames(metadata)==phenotype[4]))]

meta <- meta[-which(is.na(meta$Sleep_timing_per_day)),]

rownames(meta) <-NULL


############################################################
# Define phenotype groups
############################################################

bottom_species_names <- meta$idx[which(meta$Sleep_timing_per_day == "Sleep at night")]

top_species_names <- meta$idx[which(meta$Sleep_timing_per_day == "Sleep at daytime")]


############################################################
# Initialize objects
############################################################

adjustedRandIndexDF <- c()
scoreDF <- c()

final.result <- c()

cluster.sum.df <- c()
purity.score <- c()


############################################################
# Main gene loop
############################################################

for(i in 1:length(gene)){
  
  
  cat(paste0("-------------", gene[i], "-------------"), "\n")
  
  setwd(DATA.DIR)
  
  ## Load aligned CDS sequence for each gene
  cds_muscle <- readDNAStringSet(muscle[i])
  
  cds_muscle@ranges@NAMES <- sapply(strsplit(cds_muscle@ranges@NAMES, ":"), function(x) x[2])
  
  
  ############################################################
  # Match species between alignment and metadata
  ############################################################
  
  name <- intersect(cds_muscle@ranges@NAMES, meta$Species_symbol_name_ensembl)
  
  ID <- which(meta$Species_symbol_name_ensembl %in% name)
  
  cds_muscle <- cds_muscle[which(cds_muscle@ranges@NAMES %in% name)]
  
  cds_muscle@ranges@NAMES <- meta$idx[match(cds_muscle@ranges@NAMES, meta$Species_symbol_name_ensembl)]
  
  use_meta <- meta[match(cds_muscle@ranges@NAMES, meta$idx),]
  
  rownames(use_meta) <- NULL
  
  
  ############################################################
  # Compute genetic distances (TN93/K80 default)
  ############################################################
  
  remove_ID <- which(is.na(use_meta[,1]))
  
  if(length(remove_ID) > 0) {
    cds_muscle <- cds_muscle[-remove_ID]
  }
  
  dna_Muscle <- as.DNAbin(cds_muscle)
  
  # Larger values indicate greater evolutionary distance
  dm <- dist.dna(dna_Muscle,as.matrix = T,pairwise.deletion = T)
  
  tree <- njs(dm)
  
  dm <- cophenetic.phylo(tree)
  
  
  ############################################################
  # Hierarchical clustering of species
  ############################################################
  
  average.dm <- hclust(as.dist(dm), method="average")
  
  scoreDF <- c()
  
  
  ############################################################
  # Evaluate cluster solutions (K = 2–30)
  ############################################################
  
  for(j in 2:30){
    
    cluster.sum.df <- c()
    purity.score <- c()
    
    clusters <- cutree(average.dm, j)
    
    clusters %>% table
    clusters %>% table %>% sum
    
    
    ############################################################
    # Assign phenotype labels to clusters
    ############################################################
    
    class_assignments <- rep(NA, length(clusters))
    
    names(class_assignments) <- names(clusters)
    
    class_assignments[names(clusters) %in% top_species_names] <- "more_sleep"
    
    class_assignments[names(clusters) %in% bottom_species_names] <- "less_sleep"
    
    
    df <- data.frame(Cluster = clusters[names(class_assignments)],
                     Class = class_assignments,
                     species_name = names(class_assignments),
                     row.names = NULL)
    
    
    ############################################################
    # Contingency table for phenotype vs cluster
    ############################################################
    
    result_table <- table(df$Class,df$Cluster)
    
    
    ############################################################
    # Identify clusters with equal counts across classes
    ############################################################
    
    similar_nonzero_columns <- apply(result_table, 2, function(column) {
      all(column != 0) && length(unique(column)) == 1
    })
    
    
    column_max_row <- apply(result_table, 2, which.max)
    
    df1 <- df[-which(is.na(df$Class)),]
    
    used <- df1[which(df1$Cluster %in% names(column_max_row)),]
    
    diff <- setdiff(1:j,as.numeric(names(column_max_row)))
    
    
    if(length(diff) ==0 ){
      
      used <- df1[which(df1$Cluster %in% names(column_max_row)),]
      
      Short_group1 <- used$species_name[which(used$Class == "Bottom")]
      
      Long_group1 <- used$species_name[which(used$Class == "Top")]
      
    }else{
      
      diff_species <- df1[which(df1$Cluster %in% c(diff)),]
      
      diff_species_L <- diff_species$species_name[which(diff_species$Class == "Top")]
      
      diff_species_S <- diff_species$species_name[which(diff_species$Class == "Bottom")]
      
      used <- df1[which(df1$Cluster %in% names(column_max_row)),]
      
      Short_group <- used$species_name[which(used$Class == "Bottom")]
      
      Long_group <- used$species_name[which(used$Class == "Top")]
      
      Short_group1 <- c(Short_group,diff_species_S)
      
      Long_group1 <- c(Long_group,diff_species_L)
      
    }
    
    
    result_table.df <- as.data.frame.matrix(result_table)
    
    colnames(result_table.df) <- paste0("cluster_", colnames(result_table.df))
    
    
    ############################################################
    # Aggregate clusters and compute statistics
    ############################################################
    
    cluster.sum.df <- making.sum.cluster(result_table.df, j)
    
    purity.score <- calculate.purity(cluster.sum.df)
    
    cochran <- CochranArmitageTest(result_table.df)
    
    estimate <- fisher.test(cluster.sum.df+1)
    
    
    ############################################################
    # Store clustering evaluation results
    ############################################################
    
    scoreDF <- rbind(scoreDF,
                     data.frame(Gene=gene[i],
                                Cluster= j,
                                purity_score = purity.score,
                                p_value=estimate$p.value,
                                P.cochran = cochran$p.value,
                                Estimate.cochran = cochran$statistic))
    
    
  }
  
  final.result <-rbind(final.result, scoreDF) 
  
}

############################################################
# Identify optimal cluster number (K) for each gene
############################################################

gene.ls <- unique(final.result$Gene)

length(gene.ls)

significant_gene <- c()
signifi.result <- c()

final.signifi.result <- c()

slopDF <- c()
result.ls <- c()

for(i in 1:length(gene.ls)){
  
  p2 <- c()
  
  geneDF <- c()
  significant_drop <- c()
  
  ############################################################
  # Subset clustering results for a single gene
  ############################################################
  
  geneDF <- final.result[which(final.result$Gene == gene.ls[i]),]
  
  
  ############################################################
  # Calculate purity score change across cluster numbers
  ############################################################
  
  geneDF$DiffChange = c(NA, diff(geneDF$purity_score))
  
  drop <- which(geneDF$DiffChange <= 0)
  
  
  ############################################################
  # Robust linear model to estimate purity trend
  ############################################################
  
  rlm.result <- rlm(purity_score ~ Cluster, geneDF)
  
  
  ############################################################
  # Determine optimal K
  ############################################################
  
  if(rlm.result$coefficients[2] < 0.005){
    
    residual.max.idx <- min(drop) - 1
    
  } else {
    
    residual.max.idx <- which.max(rlm.result$residuals[1:13])
    
  }
  
  
  ############################################################
  # Store optimal cluster result
  ############################################################
  
  signifi.result <- rbind(signifi.result,
                          geneDF[(residual.max.idx),])
  
}


############################################################
# Save optimal cluster results
############################################################

fwrite(signifi.result,
       "/data/Result4/Sleep_Timing_optimal_K.tsv",
       sep = "\t",
       row.names = F,
       col.names = T)
