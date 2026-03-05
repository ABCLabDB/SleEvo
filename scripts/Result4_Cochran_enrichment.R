library(data.table)
library(dplyr)
library(tidyverse)
library(ape)
library(hierfstat)
library(Biostrings)
library(ade4)
library(adegenet)
library(ggtree)
library(patchwork)
library(phyloseq)

library(foreach)
library(doParallel)
library(phytools)
library(dendextend)
library(plotly)
library(ggfortify)
library(DescTools)
library(umap)

############################################################
# Function: making.sum.cluster
# Purpose:
# Resolve cluster ties and compute aggregated cluster counts
# for "less_sleep" and "more_sleep" groups.
############################################################

making.sum.cluster <- function(result_table.df,j){
  
  tie.idx <- 0
  
  ### Handle tie cases
  ## When cluster counts are identical across rows
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
  
  ## Calculate cluster sums
  cluster.sum <- colSums(result_table.df)
  
  ## Compute cluster ratios
  bottom.res <- result_table.df[rownames(result_table.df) == "less_sleep",]/cluster.sum
  top.res <- result_table.df[rownames(result_table.df) == "more_sleep",]/cluster.sum
  
  
  ### Replace NaN values with 0
  if(sum(is.na(bottom.res)) != 0){
    bottom.res[which(is.na(bottom.res))] <- 0
  } 
  
  if(sum(is.na(top.res)) != 0){
    top.res[which(is.na(top.res))] <- 0
  } 
  
  ### Combine ratio results
  result.oi <- rbind(bottom.res, top.res)
  
  ### Identify dominant cluster label
  max.idx <- apply(result.oi, 2,function(x){
    rownames(result.oi)[which.max(x)]
  })
  
  less.sleep.idx <- names(max.idx[max.idx == "less_sleep"])
  more.sleep.idx <- names(max.idx[max.idx == "more_sleep"])
  
  less.sleep.df <- result_table.df %>% dplyr::select(all_of(less.sleep.idx))
  more.sleep.df <- result_table.df %>% dplyr::select(all_of(more.sleep.idx))
  
  
  ### Final aggregated cluster counts
  cluster.sum.df  <- data.frame(Cluster.less = rowSums(less.sleep.df),
                                Cluster.more = rowSums(more.sleep.df))
  
  if(tie.idx != 0){
    cluster.sum.df <- cluster.sum.df + tie.idx
  }
  
  return(cluster.sum.df)
  
}


############################################################
# Set up parallel computing environment
############################################################

cores=detectCores()
cores <- 48

cl <- makeCluster(cores[1]-1) # avoid overloading the system
registerDoParallel(cl)


############################################################
# Define directories
############################################################

DATA.DIR <- "/data/Cicadian_gene_Nucleotide_Matrix/"


############################################################
# Metadata preprocessing
############################################################

meta <- as.data.frame(fread("/data/Result1/species_sleep_metadata.txt"))
phenotype <- colnames(meta)[c(8,11:13)]
meta <- meta[order(meta$Sleep_timing_per_day),c(1,2,3,which(colnames(meta)==phenotype[4]))]

if(length(which(is.na(meta$Sleep_timing_per_day))) > 0 ){
  meta <- meta[-which(is.na(meta$Sleep_timing_per_day)),]
}

rownames(meta) <-NULL

meta$Sleep_label <- NA

meta$Sleep_label[which(meta$Sleep_timing == "Sleep at daytime")] <- "more_sleep"
meta$Sleep_label[which(meta$Sleep_timing == "Sleep at night")] <- "less_sleep"

colnames(meta)[1] <- "idx"


############################################################
# Load optimal cluster information
############################################################
max.purity.df <- as.data.frame(fread("/data/Result4/Sleeptiming_optimal_K.tsv"))
colnames(max.purity.df)[1:4] <- c("Gene", "best.Cluster","purity.score","p_value")
max.purity.df$best.Cluster <- paste0("cluster_",max.purity.df$best.Cluster)


############################################################
# Initialize objects
############################################################

result.df <- c()

gene.idx <- c()
gene.cluster.idx <- c()
original.purity <- c()

Gene.cluster <- c()
best.Gene.cluster <- c()
best.cluster.df <- c()
permu.df <- c()

nsim <- 1000000
res <- c()
perm_p_values <- c()
P_odds <- c()

new.permu.df <- c()
result_table <- c()
result_table.df <- c()
cluster.sum <- c()
test.df <- c()
result.oi <- c()

p_value <- c()


############################################################
# Parallel permutation analysis across genes
############################################################

results <- foreach(k = 1:nrow(max.purity.df), .packages = c("data.table", "dplyr"), .combine = rbind) %dopar% {
  
  ## Load required libraries inside each worker
  library(ape)
  library(hierfstat)
  library(Biostrings)
  library(ade4)
  library(phytools)
  library(dendextend)
  library(plotly)
  library(ggfortify)
  library(DescTools)
  library(umap)
  
  gene.idx <- c()
  gene.cluster.idx <- c()
  
  
  Gene.cluster <- c()
  best.Gene.cluster <- c()
  best.cluster.df <- c()
  permu.df <- c()
  
  
  ## Retrieve gene-specific parameters
  gene.idx <- max.purity.df$Gene[k]
  gene.cluster.idx <- max.purity.df$best.Cluster[k]
  original.cochran <- max.purity.df$Estimate.cochran[k]
  
  cat(paste0("-------------", gene.idx, "-------------"), "\n")
  
  ## Load aligned CDS sequences
  setwd(DATA.DIR)
  cds_muscle <- readDNAStringSet(paste0(gene.idx,"_muscle.fasta"))
  cds_muscle@ranges@NAMES <- sapply(strsplit(cds_muscle@ranges@NAMES, ":"), function(x) x[2])
  
  ## Match species present in metadata
  meta$Species_symbol_name_ensembl <- gsub(" ","_",tolower(meta$Species_symbol_name_ensembl))
  name <- intersect(cds_muscle@ranges@NAMES, meta$Species_symbol_name_ensembl)
  ID <- which(meta$Species_symbol_name_ensembl %in% name)
  cds_muscle <- cds_muscle[which(cds_muscle@ranges@NAMES %in% name)]
  
  cds_muscle@ranges@NAMES <- meta$idx[match(cds_muscle@ranges@NAMES, meta$Species_symbol_name_ensembl)]
  use_meta <- meta[match(cds_muscle@ranges@NAMES, meta$idx),]
  rownames(use_meta) <- NULL
  
  remove_ID <- which(is.na(use_meta[,1]))
  
  if(length(remove_ID) > 0) {
    cds_muscle <- cds_muscle[-remove_ID]
  }
  
  ## Convert to DNAbin format
  dna_Muscle <- as.DNAbin(cds_muscle)
  
  ## Compute genetic distance and build phylogenetic tree
  dm <- dist.dna(dna_Muscle,as.matrix = T,pairwise.deletion = T)
  tree <- njs(dm)
  
  ## Compute cophenetic distance matrix
  dm <- cophenetic.phylo(tree)
  
  ## Hierarchical clustering
  average.dm <- hclust(as.dist(dm), method="average")
  
  K_group <- sapply(strsplit(gene.cluster.idx, "_"), function(x) x[2])
  clusters <- cutree(average.dm, K_group)
  
  clusters <- as.data.frame(clusters)
  clusters$idx <- rownames(clusters)
  
  colnames(clusters) <- c("clusterInfo","idx")
  clusters <- clusters[,c(2,1)]
  
  best.cluster.df <- merge(clusters, meta, by = "idx", all = FALSE, sort = FALSE)
  
  permu.df <- best.cluster.df %>% dplyr::select(idx, clusterInfo, Sleep_label)
  
  ############################################################
  # Permutation test
  ############################################################
  set.seed(319)
  Co_estimate <- numeric(nsim)
  Co_estimate_merge <- numeric(nsim)
  
  for(i in 1:nsim){
    
    perm <- sample(c(unique(permu.df$clusterInfo)), nrow(permu.df), replace = TRUE)
    
    new.permu.df <- transform(permu.df, New.clusters = perm)
    
    result_table <- table(new.permu.df$Sleep_label, new.permu.df$New.clusters)
    result_table.df <- as.data.frame.matrix(result_table)
    
    cochran_test <- CochranArmitageTest(result_table)
    
    ### Aggregate clusters
    result_table.df <- making.sum.cluster(result_table.df, as.numeric(gsub("cluster_","",gene.cluster.idx)))
    colnames(result_table.df) <- paste0("cluster_", colnames(result_table.df))
    
    cochran_merge_test <- CochranArmitageTest(result_table.df)
    
    Co_estimate[i] <- cochran_test$statistic
    Co_estimate_merge[i] <- cochran_merge_test$statistic
  }
  
  ############################################################
  # Compute permutation p-value
  ############################################################
  P.cochran <- (sum(abs(Co_estimate) >= abs(original.cochran)) + 1) / nsim
  P.cochran.merge <- (sum(abs(Co_estimate_merge) >= abs(original.cochran)) + 1) / nsim
  
  
  data.frame(gene.idx, gene.cluster.idx, P.cochran,P.cochran.merge)
}

############################################################
# Stop parallel cluster
############################################################

stopCluster(cl)

############################################################
# Save final results
############################################################

result.df <- results

fwrite(result.df,
       "/data/Result4/Sleeptiming_Cochran_Result.tsv",
       sep = "\t",
       row.names = F,
       col.names = T)
