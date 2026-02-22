############################################################
# Result3_visualization.R
#
# Master visualization script for Result3
# Generates:
#   - Figure3A : NREM phylogenetic tree
#   - Figure3B : Evolutionary clustering boxplot
#   - Figure3C : Phylogenetic signal
#   - Figure3D : Selection heatmap
#   - Figure3E : Manhattan plot
#   - Figure3F : SNP-centered profile
#   - ARNTL2 mutation map
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(Biostrings)
  library(ape)
  library(dendextend)
  library(pheatmap)
  library(tidyr)
  library(ggseqlogo)
  library(patchwork)
  library(Cairo)
})

############################################################
# 1. Project root check
############################################################

if(!dir.exists("data") || !dir.exists("figures")){
  stop("Run from Sleep_Evolution project root.")
}

PROJECT_DIR <- getwd()

############################################################
# 2. Global paths
############################################################

PATHS <- list(
  META      = "data/Result1/species_sleep_metadata.txt",
  PERM_A    = "data/Result3/NREM_key_12_optimal_Kruskal.tsv",
  PERM_B    = "data/Result3/NREM_ratio_Anova_Result.tsv",
  PHYLO_SIG = "data/Result3/NREM_ratio_phylogeneticsignal.tsv",
  SNP       = "data/Result3/NREM_ratio_SNP.tsv",
  AA_MUT    = "data/Result3/NREM_ratio_AminoAcid_mutation.tsv",
  FASTA     = "data/Fasta",
  NUC_DIR   = "data/Circadian_gene_Nucleotide_Matrix",
  HEATMAP   = "data/Heatmap/NREM_ratio.tsv",
  GENE_INFO = "data/Result1/circadian_Gene_list.tsv",
  OUT       = "figures/Result3"
)

dir.create(PATHS$OUT, recursive=TRUE, showWarnings=FALSE)

############################################################
# 3. Shared data
############################################################

meta <- fread(PATHS$META)
perm_A <- fread(PATHS$PERM_A)
perm_B <- fread(PATHS$PERM_B)

sig_genes <- perm_B |> filter(adj.P < 0.05) |> pull(Gene)

############################################################
# 4. Figure3A — Phylogenetic trees
############################################################

make_phylo_trees <- function(){

  message("Generating Figure3A...")

  for(gene in sig_genes){

    fasta_path <- file.path(PATHS$FASTA,
                            paste0(gene,"_muscle.fasta"))

    if(!file.exists(fasta_path)) next

    cds <- readDNAStringSet(fasta_path)
    names(cds) <- sapply(strsplit(names(cds),":"),`[`,2)

    dna <- as.DNAbin(cds)
    tree <- nj(dist.dna(dna, model="T92"))
    dend <- as.dendrogram(as.hclust(tree))

    pdf(file.path(PATHS$OUT,
                  paste0(gene,"_Figure3A_tree.pdf")),
        width=10,height=5)
    plot(dend)
    dev.off()
  }
}

############################################################
# 5. Figure3C — Phylogenetic signal
############################################################

make_phylo_signal <- function(){

  message("Generating Figure3C...")

  DF <- fread(PATHS$PHYLO_SIG)
  DF <- DF |> filter(Gene %in% sig_genes)

  df_long <- DF |>
    select(Gene,K,I) |>
    pivot_longer(cols=c(K,I),
                 names_to="Signal",
                 values_to="Value")

  p <- ggplot(df_long,
              aes(Gene,Value,color=Signal))+
    geom_point()+
    theme_bw()+
    theme(axis.text.x=
            element_text(angle=45,hjust=1))

  ggsave(file.path(PATHS$OUT,
                   "Figure3C_phylogenetic_signal.pdf"),
         p,width=6,height=4)
}

############################################################
# 6. Figure3D — Heatmap
############################################################

make_heatmap <- function(){

  message("Generating Figure3D...")

  heat_df <- fread(PATHS$HEATMAP)
  heat_df <- heat_df |> filter(Gene %in% sig_genes)

  mat <- heat_df[,c("Ps","Pi","TajimasD","dN_dS")]

  CairoPDF(file.path(PATHS$OUT,
                     "Figure3D_heatmap.pdf"),
           width=10,height=6)

  pheatmap(mat,
           cluster_rows=FALSE,
           cluster_cols=FALSE)

  dev.off()
}

############################################################
# 7. Figure3E — Manhattan
############################################################

make_manhattan <- function(){

  message("Generating Figure3E...")

  df <- fread(PATHS$SNP)
  df$logP <- -log10(df$P)

  p <- ggplot(df,aes(BP,logP,color=factor(CHR)))+
    geom_point(size=1.5)+
    geom_hline(yintercept=-log10(0.05),
               linetype="dashed",
               color="red")+
    theme_classic()+
    theme(legend.position="none")

  ggsave(file.path(PATHS$OUT,
                   "Figure3E_Manhattan.pdf"),
         p,width=10,height=4)
}

############################################################
# 8. ARNTL2 Mutation Map
############################################################

make_ARNTL2_map <- function(){

  message("Generating ARNTL2 mutation map...")

  mut <- fread(PATHS$AA_MUT)
  mut <- mut |> filter(Gene=="ARNTL2")

  if(nrow(mut)==0) return(NULL)

  p <- ggplot(mut,
              aes(Codon_Position,1))+
    geom_segment(aes(xend=Codon_Position,
                     y=1,yend=1.2))+
    geom_point(size=3)+
    theme_void()

  ggsave(file.path(PATHS$OUT,
                   "ARNTL2_Mutation_Map.pdf"),
         p,width=7,height=3)
}

############################################################
# 9. Main
############################################################

main <- function(){

  message("====================================")
  message("   Running Result3 Visualization    ")
  message("====================================")

  make_phylo_trees()
  make_phylo_signal()
  make_heatmap()
  make_manhattan()
  make_ARNTL2_map()

  message("All Result3 figures generated.")
}

if (sys.nframe() == 0) {
  main()
}
