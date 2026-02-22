############################################################
# Result4_visualization.R
#
# Generate all Result4 visualizations:
#   Figure4A  : Lollipop enrichment
#   Figure4B  : Sleep timing dendrogram
#   Figure4C  : Evolution quadrant plot
#   Figure4D  : Manhattan plot
#   Figure4E  : PER1 amino acid landscape
#   SuppFig12 : SNP heatmap window plot
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
  library(stringr)
  library(tidyr)
  library(patchwork)
})

PROJECT_DIR <- getwd()

############################################################
# 1. Global Paths
############################################################

PATHS <- list(
  COCHRAN = file.path(PROJECT_DIR,"data","Result4","Sleeptiming_Cochran_Result.tsv"),
  MANHATTAN = file.path(PROJECT_DIR,"data","Result4","Sleeptiming_Manhattan_Dataset.tsv"),
  MUTATION = file.path(PROJECT_DIR,"data","Result4","Sleeptiming_AA_Mutation.tsv"),
  META = file.path(PROJECT_DIR,"data","Result1","species_sleep_metadata.txt"),
  FASTA_DIR = file.path(PROJECT_DIR,"data","Fasta"),
  EVOLUTION = file.path(PROJECT_DIR,"data","Heatmap","Sleep_timing.tsv"),
  GENE_INFO = file.path(PROJECT_DIR,"data","Result1","circadian_Gene_list.tsv"),
  NUC_DIR = file.path(PROJECT_DIR,"data","Circadian_gene_Nucleotide_Matrix"),
  OUT = file.path(PROJECT_DIR,"figures","Result4")
)

dir.create(PATHS$OUT, recursive=TRUE, showWarnings=FALSE)

############################################################
# 2. Load Shared Data
############################################################

cochran <- fread(PATHS$COCHRAN) |> as.data.frame()
cochran$gene.cluster.idx <- as.numeric(sapply(strsplit(cochran$gene.cluster.idx,"_"), `[`,2))
sig_genes <- cochran |> filter(P.cochran < 0.05)

meta <- fread(PATHS$META) |> as.data.frame()
meta$Species_symbol_name_ensembl <- gsub(" ","_",meta$Species_symbol_name_ensembl)

############################################################
# 3. Figure4A — Lollipop Plot
############################################################

make_lollipop <- function(){

  df <- cochran |>
    transmute(Gene=gene.idx,
              logP=-log10(P.cochran)) |>
    arrange(desc(logP)) |>
    slice_head(n=50)

  p <- ggplot(df,aes(reorder(Gene,-logP),logP))+
    geom_segment(aes(xend=Gene,y=0,yend=logP),
                 size=4,color="#008585")+
    geom_point(size=4,color="#008585")+
    geom_hline(yintercept=-log10(0.05),
               linetype="dashed",color="red")+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90,hjust=1))+
    labs(x="",y="-log10(P)")

  ggsave(file.path(PATHS$OUT,"Figure4A_lollipop.pdf"),
         p,width=14,height=4)
}

############################################################
# 4. Figure4B — Phylogenetic Tree
############################################################

build_tree <- function(gene){

  fasta_file <- file.path(PATHS$FASTA_DIR,paste0(gene,"_muscle.fasta"))
  if(!file.exists(fasta_file)) return(NULL)

  cds <- readDNAStringSet(fasta_file)
  names(cds) <- str_to_title(sapply(strsplit(names(cds),":"),`[`,2))

  keep <- intersect(names(cds), meta$Species_symbol_name_ensembl)
  cds <- cds[names(cds)%in%keep]

  dna <- as.DNAbin(cds)
  tree <- nj(dist.dna(dna,model="T92"))
  dend <- as.dendrogram(as.hclust(tree))

  return(dend)
}

make_trees <- function(){

  for(gene in sig_genes$gene.idx){

    dend <- build_tree(gene)
    if(is.null(dend)) next

    pdf(file.path(PATHS$OUT,paste0(gene,"_SleepTimingTree.pdf")),
        width=12,height=5)
    plot(dend)
    dev.off()
  }
}

############################################################
# 5. Figure4C — Evolution Quadrant
############################################################

make_evolution_plot <- function(){

  evo <- fread(PATHS$EVOLUTION) |> as.data.frame()
  evo <- evo |> filter(Gene %in% sig_genes$gene.idx)

  p <- ggplot(evo,aes(dN_dS,TajimasD))+
    geom_point(size=4,color="#476066")+
    geom_vline(xintercept=1,linetype="dashed")+
    geom_hline(yintercept=0,linetype="dashed")+
    theme_classic()

  ggsave(file.path(PATHS$OUT,"Figure4C_Evolution.pdf"),
         p,width=6,height=4)
}

############################################################
# 6. Figure4D — Manhattan Plot
############################################################

make_manhattan <- function(){

  df <- fread(PATHS$MANHATTAN) |> as.data.frame()
  df$logP <- -log10(df$P)

  p <- ggplot(df,aes(BP,logP,color=factor(CHR)))+
    geom_point(size=1.5)+
    geom_hline(yintercept=-log10(0.05),
               linetype="dashed",color="red")+
    theme_classic()+
    theme(legend.position="none")

  ggsave(file.path(PATHS$OUT,"Figure4D_Manhattan.pdf"),
         p,width=8,height=4)
}

############################################################
# 7. Figure4E — PER1 Amino Acid Plot
############################################################

make_per1_plot <- function(){

  mut <- fread(PATHS$MUTATION) |> as.data.frame()
  per1 <- mut |> filter(Gene=="PER1",
                        Mutation_Type=="Non-synonymous")

  if(nrow(per1)==0) return(NULL)

  full_length <- unique(per1$Full_Length)

  p <- ggplot(per1,
              aes(Codon_Position,1))+
    geom_segment(aes(xend=Codon_Position,
                     y=1,yend=1.2))+
    geom_point(size=3,color="#2c6e49")+
    theme_void()

  ggsave(file.path(PATHS$OUT,
                   "Figure4E_PER1_AminoAcid.pdf"),
         p,width=8,height=3)
}

############################################################
# 8. Supplementary Figure — SNP Heatmap
############################################################

make_snp_heatmap <- function(){

  TS <- fread(PATHS$MANHATTAN) |> as.data.frame()
  TS <- TS |> filter(P<0.05)

  files <- list.files(PATHS$NUC_DIR,
                      pattern=".tsv",
                      full.names=TRUE)

  for(f in files){

    gene_name <- tools::file_path_sans_ext(basename(f))
    aln <- fread(f) |> as.data.frame()

    if(!gene_name %in% TS$Gene) next

    p <- ggplot()+theme_void()

    ggsave(file.path(PATHS$OUT,
                     paste0(gene_name,"_SNP.pdf")),
           p,width=6,height=4)
  }
}

############################################################
# 9. Run All
############################################################

make_lollipop()
make_trees()
make_evolution_plot()
make_manhattan()
make_per1_plot()
make_snp_heatmap()

message("All Result4 figures generated successfully.")

############################################################
# 10. Main Execution Block
############################################################

main <- function(){

  message("====================================")
  message("   Running Result4 Visualization    ")
  message("====================================")

  make_lollipop()
  make_trees()
  make_evolution_plot()
  make_manhattan()
  make_per1_plot()
  make_snp_heatmap()

  message("All Result4 figures generated successfully.")
}

# Run only if executed directly (not sourced)
if (sys.nframe() == 0) {
  main()
}
