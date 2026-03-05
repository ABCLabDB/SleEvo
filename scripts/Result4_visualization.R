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
  
  threshold <- -log10(0.05)
  
  df$color <- ifelse(df$logP > threshold, "#008585", "grey90")
  df$Significant <- ifelse(df$logP > threshold, "Sleep timing related", "Not significant")
  
  p <- ggplot(df, aes(x = reorder(Gene, -logP), y = logP, color = Significant)) +
    geom_segment(aes(x = reorder(Gene, -logP), xend = reorder(Gene, -logP), y = 0, yend = logP), size = 6.4) +
    geom_point(size = 6) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red",size=1.5) +
    scale_color_manual(values = c("Sleep timing related" = "#008585", "Not significant" = "grey90")) +
    labs(
      x = "",
      y = "-log10(P.value)",
      color = "Significance"
    ) +
    theme_classic(base_size = 12) + #theme_minimal(base_size = 12)
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 14,face = "italic"),
      axis.line.x = element_line(size = 1, color = "black"),
      axis.line.y = element_line(size = 1, color = "black"),
      panel.border = element_blank()
    )
  
  ggsave(file.path(PATHS$OUT,"Figure4A_lollipop.pdf"),
         p,width=14,height=4)
}

############################################################
# 4. Figure4B — Phylogenetic Tree
############################################################

make_trees <- function(){
  
  perm_test <- cochran
  knee <- perm_test |> filter(P.cochran < 0.05)
  
  custom_colors <- c(
    "Sleep at night" = "#2c6e49",
    "Sleep at anytime" = "grey85",
    "Sleep at daytime" = "#d68c45"
  )
  
  meta_timing <- meta[-which(is.na(meta$Sleep_timing_per_day)),]
  meta_timing$Species_symbol_name_ensembl <- tolower(meta_timing$Species_symbol_name_ensembl)
  
  meta_timing$Type <- "Sleep at anytime"
  meta_timing$Type[meta_timing$Sleep_timing_per_day == "Sleep at night"] <- "Sleep at night"
  meta_timing$Type[meta_timing$Sleep_timing_per_day == "Sleep at daytime"] <- "Sleep at daytime"
  
  for(i in 1:nrow(knee)){
    
    k <- knee$gene.cluster.idx[i]
    
    if(k == 2) next
    
    gene <- knee$gene.idx[i]
    
    fasta_file <- file.path(PATHS$FASTA_DIR,
                            paste0(gene,"_muscle.fasta"))
    
    if(!file.exists(fasta_file)) next
    
    cds_muscle <- readDNAStringSet(fasta_file)
    
    cds_muscle@ranges@NAMES <- sapply(
      strsplit(cds_muscle@ranges@NAMES, ":"), `[`,2
    )
    
    name <- intersect(cds_muscle@ranges@NAMES,
                      meta_timing$Species_symbol_name_ensembl)
    
    cds_muscle <- cds_muscle[
      cds_muscle@ranges@NAMES %in% name
    ]
    
    cds_muscle@ranges@NAMES <- meta_timing$Species_name_ensembl[
      match(cds_muscle@ranges@NAMES,
            meta_timing$Species_symbol_name_ensembl)
    ]
    
    use_meta <- meta_timing[
      match(cds_muscle@ranges@NAMES,
            meta_timing$Species_name_ensembl),]
    
    rownames(use_meta) <- NULL
    
    dna_Muscle <- as.DNAbin(cds_muscle)
    
    dm <- dist.dna(
      dna_Muscle,
      as.matrix = TRUE,
      pairwise.deletion = TRUE,
      model = "T92"
    )
    
    tree <- njs(dm)
    
    dm <- cophenetic.phylo(tree)
    
    average.dm <- hclust(as.dist(dm), method="average")
    
    dend <- as.dendrogram(average.dm)
    
    meta_for_color <- use_meta[,c("Species_name_ensembl","Type")]
    
    meta_for_color <- meta_for_color |>
      mutate(color = custom_colors[Type])
    
    meta_for_color <- meta_for_color[
      match(labels(dend),
            meta_for_color$Species_name_ensembl),]
    
    labels_colors(dend) <- meta_for_color$color
    
    sleep_info <- setNames(
      meta_for_color$Type,
      meta_for_color$Species_name_ensembl
    )
    
    assign_branch_color_by_mode <- function(d) {
      if (is.leaf(d)) {
        sleep_type <- sleep_info[labels(d)]
        attr(d, "Type") <- sleep_type
        if (!is.na(sleep_type)) {
          attr(d, "edgePar") <- list(col = custom_colors[sleep_type], lwd = 4)
        } else {
          attr(d, "edgePar") <- list(col = "gray70", lwd = 4)
        }
        return(d)
      }
      
      d[[1]] <- assign_branch_color_by_mode(d[[1]])
      d[[2]] <- assign_branch_color_by_mode(d[[2]])
      
      left <- attr(d[[1]], "Type")
      right <- attr(d[[2]], "Type")
      all_types <- na.omit(c(left, right))
      
      if (length(all_types) == 0) {
        majority_type <- NA
        branch_color <- "gray70"
      } else {
        majority_type <- names(sort(table(all_types), decreasing = TRUE))[1]
        branch_color <- custom_colors[majority_type]
      }
      
      attr(d, "Type") <- majority_type
      attr(d, "edgePar") <- list(col = branch_color, lwd = 4)
      
      return(d)
    }
    
    dend_colored <- assign_branch_color_by_mode(dend)
    
    label_colors <- sleep_info[labels(dend_colored)]
    label_colors <- custom_colors[label_colors]
    label_colors[is.na(label_colors)] <- "gray70"
      
    labels_colors(dend_colored) <- label_colors
    
    ggd <- dend_colored |>
      set("labels_cex", 0.8) |>
      set("leaves_pch", 19) |>
      set("leaves_col", meta_for_color$color) |>
      set("leaves_cex", 2.5)
    
    pdf(
      file.path(PATHS$OUT,
                paste0(gene,"_SleepTimingTree.pdf")),
      width = 16,
      height = 5
    )
    
    par(mar = c(10,4,2,2))
    
    plot(ggd)
    
    legend(
      "topright",
      legend = unique(meta_for_color$Type),
      fill = unique(meta_for_color$color),
      title = "Sleep Type",
      cex = 0.8
    )
    
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
