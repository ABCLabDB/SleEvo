############################################################
# Script: Result3_visualization.R
# Description:
#   A: Phylogenetic dendrogram (NREM ratio)
#   B: Phylogenetic signal visualization
#   C: Evolutionary pressure scatter plot
#   D: Manhattan plot of NREM-associated SNPs
# Logic unchanged. Data loading centralized.
############################################################

############################
# 1. Load Required Libraries
############################

library(data.table)
library(dplyr)
library(ape)
library(Biostrings)
library(dendextend)
library(ggdendro)
library(ggh4x)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)
library(grid)
library(ggrepel)
library(qqman)
library(scales)

############################
# 2. Define File Paths (Modify as Needed)
############################

BASE_DIR        <- "/disk4/bijsy/2.Sleep/"
META_FILE       <- paste0(BASE_DIR, "2.Result/1.Result1/METADATA_fixedVersion.txt")
FASTA_DIR       <- "/disk4/bijsy/Evolution/0.DATA/2.OUT/3.CDS/1.Muscle/"
NWK_DIR         <- "/disk4/bijsy/MEGA/MEGA/"
PERM_FILE       <- paste0(BASE_DIR, "Anova_enrichment_test/NREM_key_12_optimal_Kruskal.tsv")
SUPP8_FILE      <- paste0(BASE_DIR, "Figure/UPGMA_Bootstrapping/Evolution/Supplementary_table8.tsv")
REM_FILE        <- paste0(BASE_DIR, "2.Result/2.Result2/5.Heatmap/Supplemantary_Table/REM.tsv")
CIRCADIAN_FILE  <- paste0(BASE_DIR, "Figure/Result1/supplementaryS2.tsv")
SNP_FILE        <- paste0(BASE_DIR, "Figure/Result6/NREM_Ratio_SNP.tsv")

OUTPUT_DIR      <- paste0(BASE_DIR, "Figure/Result3/Figure3_A/1.Dend_version/New_Version/")

############################
# 3. Load Metadata
############################

metadata <- fread(META_FILE)

metadata$Species_symbol_name_ensembl <-
  tolower(gsub(" ", "_", metadata$Species_symbol_name_ensembl))

meta <- metadata

# Sort species by NREM ratio
meta <- meta[order(meta$Percentage_of_NREM_time_per_day), ]

# Remove missing phenotype
meta <- meta[!is.na(meta$Percentage_of_NREM_time_per_day), ]

rownames(meta) <- NULL

############################################################
# =====================  Result3 A  ========================
#        Phylogenetic Dendrogram by NREM Ratio
############################################################

prop_rate <- 0.33

bottom_species_names <- meta[1:ceiling(prop_rate * nrow(meta)),]$Species_name_ensembl
top_species_names <- meta[(nrow(meta) - ceiling(prop_rate * nrow(meta)) + 1):nrow(meta)]$Species_name_ensembl

meta$Type <- "Others"
meta$Type[meta$Species_name_ensembl %in% bottom_species_names] <- "Low Ratio"
meta$Type[meta$Species_name_ensembl %in% top_species_names]    <- "High Ratio"

custom_colors <- c(
  "High Ratio" = "#194a7a",
  "Others"     = "grey",
  "Low Ratio"  = "#E7B800"
)

############################
# 5. Load Significant Genes
############################

perm_test <- fread(PERM_FILE)
knee <- perm_test[perm_test$adj.P < 0.05, ]
knee$Cluster <- as.numeric(knee$Cluster)

############################
# 6. Filter Genes with Available Trees
############################

fasta_file <- list.files(FASTA_DIR, ".fasta")
gene_names <- sapply(strsplit(fasta_file,"_"),
                     function(x) x[1])

############################
# 7. Dendrogram Plot per Gene
############################

for(i in 1:length(knee$Gene)){
  
  target_gene <- knee$Gene[i]
  
  if(!(target_gene %in% gene_names)){
    next
  }
  
  ############################
  # Load Alignment
  ############################
  
  cds_muscle <- readDNAStringSet(
    paste0(FASTA_DIR, target_gene, "_muscle.fasta")
  )
  
  cds_muscle@ranges@NAMES <-
    sapply(strsplit(cds_muscle@ranges@NAMES, ":"),
           function(x) x[2])
  
  valid_species <- intersect(
    cds_muscle@ranges@NAMES,
    meta$Species_symbol_name_ensembl
  )
  
  cds_muscle <- cds_muscle[
    cds_muscle@ranges@NAMES %in% valid_species
  ]
  
  cds_muscle@ranges@NAMES <-
    meta$Species_name_ensembl[
      match(cds_muscle@ranges@NAMES,
            meta$Species_symbol_name_ensembl)
    ]
  
  use_meta <- meta[
    match(cds_muscle@ranges@NAMES,
          meta$Species_name_ensembl),
  ]
  
  ############################
  # Construct Distance Matrix
  ############################
  
  dna_bin <- as.DNAbin(cds_muscle)
  dm <- dist.dna(dna_bin,
                 as.matrix = TRUE,
                 pairwise.deletion = TRUE)
  
  tree <- njs(dm)
  
  cophenetic_dm <- cophenetic.phylo(tree)
  
  average_clust <-
    hclust(as.dist(cophenetic_dm),
           method = "average")
  
  dend <- as.dendrogram(average_clust)
  
  ############################
  # Assign Leaf Colors
  ############################
  
  meta_for_color <- use_meta[, c("Species_name_ensembl", "Type")]
  
  sleep_info <-
    setNames(meta_for_color$Type,
             meta_for_color$Species_name_ensembl)
  
  assign_branch_color_by_mode <- function(d){
    
    if(is.leaf(d)){
      sleep_type <- sleep_info[labels(d)]
      attr(d,"Type") <- sleep_type
      attr(d,"edgePar") <- list(
        col = custom_colors[sleep_type],
        lwd = 3
      )
      return(d)
    }
    
    d[[1]] <- assign_branch_color_by_mode(d[[1]])
    d[[2]] <- assign_branch_color_by_mode(d[[2]])
    
    left  <- attr(d[[1]], "Type")
    right <- attr(d[[2]], "Type")
    
    majority_type <-
      names(sort(table(c(left,right)),
                 decreasing=TRUE))[1]
    
    attr(d,"Type") <- majority_type
    attr(d,"edgePar") <- list(
      col = custom_colors[majority_type],
      lwd = 3
    )
    
    return(d)
  }
  
  dend_colored <- assign_branch_color_by_mode(dend)
  
  labels_colors(dend_colored) <-
    custom_colors[sleep_info[labels(dend_colored)]]
  
  ############################
  # Plot and Save
  ############################
  
  pdf(paste0(OUTPUT_DIR,
             target_gene,
             "_nobarplot_Change.pdf"),
      width = 12, height = 6)
  
  par(mar = c(10,4,2,2))
  plot(dend_colored)
  
  legend("topright",
         legend = names(custom_colors),
         fill = custom_colors,
         title = "Sleep Type",
         cex = 0.8)
  
  dev.off()
}


############################################################
# =====================  Result3 B  ========================
#        Phylogenetic Signal Visualization
############################################################
supp8 <- fread(SUPP8_FILE)
DF <- supp8
DF[is.na(DF)] <- 0

colnames(DF)[c(2,6)] <- c("Blomberg’s K", "Moran’s I")

top30_genes <- sort(unique(DF$Gene))[1:20]
DF <- DF %>% filter(Gene %in% top30_genes)

df_long <- DF %>%
  dplyr::select(Gene, `Blomberg’s K`, `Moran’s I`) %>%
  pivot_longer(
    cols = c(`Blomberg’s K`, `Moran’s I`),
    names_to = "Signal",
    values_to = "Value"
  )

Phylogenetic_Signal_barplot <-
  ggplot(df_long, aes(x = Gene, y = Value, color = Signal)) +
  geom_line(aes(group = Signal), size = 2, alpha = 0.4, color = "grey") +
  geom_point(size = 4.5, alpha = 0.9) +
  geom_smooth(se = TRUE, method = "lm",
              size = 1.2, color = "#95be8d") +
  scale_color_manual(values = c(
    "Blomberg’s K" = "#7593af",
    "Moran’s I"    = "#730220"
  )) +
  theme_bw() +
  labs(
    x = "",
    y = "Phylogenetic signal score",
    color = "Phylogenetic signal"
  ) +
  coord_flip() +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )

ggsave(paste0(OUTPUT_SIGNAL,"Figure4B.jpg"),
       plot = Phylogenetic_Signal_barplot,
       width = 13, height = 3.5,
       units = "in", dpi = 300)

ggsave(paste0(OUTPUT_SIGNAL,"Figure4B.pdf"),
       plot = Phylogenetic_Signal_barplot,
       width = 3.5, height = 7)


############################################################
# =====================  Result3 C  ========================
#   Evolutionary Pressure (dN/dS vs Tajima's D)
############################################################
evolution_score <- fread(REM_FILE)

sleep_candidate <- supp7
sleep_candidate <- sleep_candidate[sleep_candidate$adj.P < 0.05,]
sleep_gene <- sleep_candidate$Gene

evolution_score_DF <- evolution_score
evolution_score_DF <-
  evolution_score_DF[evolution_score_DF$Gene %in% sleep_gene,]
rownames(evolution_score_DF) <- NULL

evolution_score_DF1 <- evolution_score_DF %>%
  mutate(
    dN_dS_Group = ifelse(dN_dS > 1, "dN/dS > 1", "dN/dS <= 1"),
    TajimasD_Group = ifelse(TajimasD > 0,
                            "Tajima's D > 0",
                            "Tajima's D <= 0")
  )

highlight_genes <- c("ARNTL2", "CPT1A", "ATF4", "PPARA", "CRY1")

p1 <- ggplot(evolution_score_DF1,
             aes(x = dN_dS, y = TajimasD)) +
  
  geom_rect(aes(xmin = -Inf, xmax = 1, ymin = 0, ymax = Inf),
            fill = NA, color = "#ab5852", size=1.5) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = Inf),
            fill = NA, color = "#7593af", size=1.5) +
  geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 0),
            fill = NA, color = "#d69e49", size=1.5) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = -Inf, ymax = 0),
            fill = NA, color = "#eadaa0", size=1.5) +
  
  geom_point(
    data = evolution_score_DF1[!evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(color = interaction(dN_dS_Group, TajimasD_Group)),
    size = 6.5, alpha = 0.5
  ) +
  
  geom_point(
    data = evolution_score_DF1[evolution_score_DF1$Gene %in% highlight_genes, ],
    aes(color = interaction(dN_dS_Group, TajimasD_Group)),
    size = 6.5, alpha = 1
  ) +
  
  geom_text_repel(
    data = evolution_score_DF1,
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf
  ) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray",
             size=1) +
  
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "gray",
             size=1) +
  
  scale_color_manual(
    values = c(
      "dN/dS > 1.Tajima's D > 0"  = "#7593af",
      "dN/dS > 1.Tajima's D <= 0" = "#eadaa0",
      "dN/dS <= 1.Tajima's D > 0" = "#ab5852",
      "dN/dS <= 1.Tajima's D <= 0"= "#d69e49"
    ),
    name = "Type"
  ) +
  
  labs(
    x = "dN/dS",
    y = "Tajima's D"
  ) +
  
  theme_minimal()

ggsave(paste0(OUTPUT_SIGNAL,"Figure4C.jpg"),
       plot = p1,
       width = 6, height = 3.5,
       units = "in", dpi = 300)

ggsave(paste0(OUTPUT_SIGNAL,"Figure4C.pdf"),
       plot = p1,
       width = 10, height = 6)


############################################################
# =====================  Result3 D  ========================
#                    Manhattan Plot
############################################################

circadian_Chr   <- fread(CIRCADIAN_FILE)
TS_site1_raw    <- fread(SNP_FILE)

TS_site1 <- as.data.frame(TS_site1_raw)

TS_site1$Chr <- circadian_Chr$Chr[
  match(TS_site1$Gene, circadian_Chr$Gene_symbol)
]

IDX <- which(grepl("CDS",TS_site1$Column_Index))
TS_site1 <- TS_site1[IDX,]
rownames(TS_site1) <- NULL

TS_site1$Column_Index <- sapply(
  strsplit(TS_site1$Column_Index, "CDS_"),
  function(x) x[2]
)

TS_site1$Column_Index <-
  paste0("Chr",TS_site1$Chr,":",TS_site1$Column_Index)

TS_site1 <- TS_site1[-which(is.na(TS_site1$p.adj)),]
TS_site1 <- TS_site1[,c(1,10,4,8,5)]
colnames(TS_site1) <- c("SNP", "CHR", "BP","P","Gene")

TS_site1$CHR[TS_site1$CHR %in% "X"] <- 23
TS_site1$CHR <- as.numeric(TS_site1$CHR)

TS_site1 <- as.data.table(TS_site1)
TS_site1[, FDR := P]
TS_site1[, CHR := as.numeric(CHR)]
TS_site1[, Label := ifelse(-log10(P) > 3, SNP, NA_character_)]

don <- TS_site1 %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(TS_site1, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(
    BPcum = BP + tot,
    is_annotate = ifelse(!is.na(Label), "yes", "no")
  )

axisdf <- don %>%
  group_by(CHR) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2)

chr_colors <- rep(c("#1e3d58","#43b0f1"), 12)
sig_threshold <- -log10(0.05)

P1 <- ggplot(don, aes(x = BPcum, y = -log10(P))) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(-log10(don$P)) + 1)
  ) +
  geom_hline(
    yintercept = sig_threshold,
    color = "red",
    linetype = "dashed",
    size = 1.5
  ) +
  geom_text_repel(
    data = subset(don,
                  -log10(P) > sig_threshold &
                    Gene %in% c("ARNTL2","PPARA")),
    aes(label = SNP),
    size = 4,
    max.overlaps = 10
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](P))
  )

ggsave(paste0(OUTPUT_SIGNAL,"Figure4E.pdf"),
       P1, width = 14, height = 4, dpi = 300)

ggsave(paste0(OUTPUT_SIGNAL,"Figure4E.jpg"),
       P1, width = 10, height = 4, dpi = 300)

############################################################
# End of Script
############################################################
