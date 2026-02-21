############################################################
# Result4_Figure_SleepTiming.R
#
# A) Lollipop plot of Fisher enrichment (-log10 P)
# B) Sleep timing concordance matrix based on phylo clustering
#
# Input:
#   data/Result4/Supplementary_Table_S23.tsv
#   data/Fasta/*.fasta
#   data/Metadata/species_sleep_metadata.txt
#
# Output:
#   figures/Result4/Figure4A_lollipop.pdf
#   figures/Result4/Figure4A/<GENE>.pdf
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
})

PROJECT_DIR <- getwd()

## =========================================================
## 1. Paths
## =========================================================

PERM_FILE <- file.path(PROJECT_DIR,
                       "data","Result4",
                       "Sleeptiming_Cochran_Result.tsv")

META_FILE <- file.path(PROJECT_DIR,
                       "data","Result1",
                       "species_sleep_metadata.txt")

FASTA_DIR <- file.path(PROJECT_DIR,
                       "data","Fasta")

OUT_DIR_A <- file.path(PROJECT_DIR,
                       "figures","Result4A_1")

OUT_DIR_B <- file.path(PROJECT_DIR, "figures","Figure4A_2")

dir.create(OUT_DIR_A, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_B, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. Load enrichment data (Figure4A_1)
## =========================================================

perm_test <- fread(PERM_FILE) |> as.data.frame()

perm_test$gene.cluster.idx <-
  sapply(strsplit(perm_test$gene.cluster.idx,"_"),
         `[`,2)

perm_test$gene.cluster.idx <- as.numeric(perm_test$gene.cluster.idx)

Data <- perm_test %>%
  transmute(
    Gene,
    Cluster = gene.cluster.idx,
    P_Value = P.cochran,
    logP = -log10(P_Value)
  )

threshold <- -log10(0.05)

top50 <- Data %>%
  arrange(desc(logP)) %>%
  slice_head(n=50)

p_lollipop <- ggplot(top50,
                     aes(x=reorder(Gene,-logP),
                         y=logP)) +
  geom_segment(aes(xend=Gene,
                   y=0, yend=logP),
               size=5,
               color="#008585") +
  geom_point(size=5,color="#008585") +
  geom_hline(yintercept=threshold,
             linetype="dashed",
             color="red",
             linewidth=1.2) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_text(angle=90,
                             hjust=1,
                             size=12)
  ) +
  labs(x="",
       y="-log10(P value)")

ggsave(file.path(OUT_DIR_A,
                 "Figure4A_lollipop.pdf"),
       p_lollipop,
       width=14,height=4)

## =========================================================
## 3. Load species metadata
## =========================================================

meta <- fread(META_FILE) |> as.data.frame()

meta$Species_symbol_name_ensembl <-
  gsub(" ","_",
       meta$Species_symbol_name_ensembl)

meta <- meta %>%
  filter(!is.na(Sleep_timing_per_day))

meta$type <- case_when(
  Sleep_timing_per_day == "Sleep at night"   ~ "night",
  Sleep_timing_per_day == "Sleep at daytime" ~ "daytime",
  TRUE                                       ~ "anytime"
)

## =========================================================
## 4. Significant genes
## =========================================================

sig_genes <- Data %>%
  filter(P_Value < 0.05)

## =========================================================
## 5. Helper: cluster assignment
## =========================================================

assign_sleep_cluster <- function(gene, k){

  fasta_file <- file.path(FASTA_DIR,
                          paste0(gene,"_muscle.fasta"))

  if(!file.exists(fasta_file)) return(NULL)

  cds <- readDNAStringSet(fasta_file)

  names(cds) <-
    sapply(strsplit(names(cds),":"),
           `[`,2)

  names(cds) <- str_to_title(names(cds))

  keep <- intersect(names(cds),
                    meta$Species_symbol_name_ensembl)

  cds <- cds[names(cds) %in% keep]

  names(cds) <- meta$Species_name_ensembl[
    match(names(cds),
          meta$Species_symbol_name_ensembl)
  ]

  dna <- as.DNAbin(cds)

  dm <- dist.dna(dna,
                 model="T92",
                 pairwise.deletion=TRUE)

  tree <- njs(dm)

  hc <- hclust(as.dist(cophenetic(tree)),
               method="average")

  clusters <- cutree(hc, k=k)

  df <- data.frame(
    Species = names(clusters),
    Cluster = clusters
  )

  df$Real <- meta$Sleep_timing_per_day[
    match(df$Species,
          meta$Species_name_ensembl)
  ]

  return(df)
}

## =========================================================
## 6. Loop genes (Figure4B_2)
## =========================================================

for(i in seq_len(nrow(sig_genes))){

  gene <- sig_genes$Gene[i]
  k    <- sig_genes$Cluster[i]

  if(k <= 1) next

  df <- assign_sleep_cluster(gene,k)

  if(is.null(df)) next

  df <- df %>%
    mutate(match = Real %in%
             c("Sleep at daytime",
               "Sleep at night"))

  p <- ggplot(df,
              aes(x=factor(Cluster),
                  y=match,
                  color=Real)) +
    geom_jitter(width=0.2,
                height=0.1,
                size=4) +
    scale_color_manual(
      values=c(
        "Sleep at night"="#2c6e49",
        "Sleep at daytime"="#d68c45",
        "Sleep at anytime"="grey80"
      )
    ) +
    theme_minimal() +
    labs(title=gene,
         x="Cluster",
         y="Sleep Timing Related")

  ggsave(file.path(OUT_DIR_B,
                   paste0(gene,".pdf")),
         p,
         width=8,height=4)
}

message("Result4 figures saved.")
