############################################################
# Result4_FigE_PER1.R
#
# Amino acid mutation landscape for PER1
# Sleep timing–associated non-synonymous mutations
#
# Input:
#   data/Result4/Sleeptiming_AA_Mutation.tsv
#   data/Result4/Sleeptiming_Cochran_Result.tsv
#   data/Result1/circadian_Gene_list.tsv
#
# Output:
#   figures/Result4/Figure4E_PER1_AminoAcid.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# ----------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------

# Amino acid mutation table (그림 경로 반영)
mutation_results <- fread(
  "data/Result4/Sleeptiming_AA_Mutation.tsv"
) %>% as.data.frame()

# Cochran significant genes
perm_test <- fread(
  "data/Result4/Sleeptiming_Cochran_Result.tsv"
) %>% as.data.frame()

timing_gene <- perm_test %>%
  filter(P.cochran < 0.05)

# Chromosome info
circadian_Chr <- fread(
  "data/Result1/circadian_Gene_list.tsv"
) %>% as.data.frame()

# ----------------------------------------------------------
# 2. Preprocessing
# ----------------------------------------------------------

mutation_results <- mutation_results[,-c(1:2)]

mutation_results$Position <-
  sapply(strsplit(mutation_results$Position,"_"),
         `[`, 2)

mutation_results$Chr <-
  circadian_Chr$Chr[
    match(mutation_results$Gene,
          circadian_Chr$Gene_symbol)
  ]

mutation_results$Position <-
  paste0("Chr",
         mutation_results$Chr,
         ":",
         mutation_results$Position)

# timing gene + non-synonymous only
mutation_results <- mutation_results %>%
  filter(Gene %in% timing_gene$gene.idx,
         Mutation_Type == "Non-synonymous")

# ----------------------------------------------------------
# 3. PER1 extraction
# ----------------------------------------------------------

protein_Gene <- mutation_results %>%
  filter(Gene == "PER1") %>%
  mutate(
    mutation_label =
      paste0(Original_Amino_Acid,
             Codon_Position,
             Mutated_Amino_Acid,
             " (", Position, ")"),
    y =
      ifelse(Sleep_Timing ==
               "Sleep at night", 1.1, 1),
    yend =
      ifelse(Sleep_Timing ==
               "Sleep at night",
             runif(n(), 1.3, 1.38),
             runif(n(), 0.65, 0.75))
  )

full_length <- unique(protein_Gene$Full_Length)

# ----------------------------------------------------------
# 4. Domain annotation (PER1-specific)
# ----------------------------------------------------------

domain_df <- data.frame(
  Domain = c(
    "CRY binding",
    "Nuclear export signal",
    "PAS", "PAS",
    "PAC",
    "Nuclear export signal",
    "Nuclear localization signal",
    "Nuclear export signal",
    "LXXLL"
  ),
  start = c(
    1149, 138,
    208, 348,
    422,
    489, 827,
    982, 1043
  ),
  end = c(
    1290, 147,
    275, 414,
    465,
    498, 843,
    989, 1047
  ),
  fill = rep(c("#fbf2c4","#74a892"),
             length.out = 9)
)

# ----------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------

p1 <- ggplot() +
  geom_rect(
    aes(xmin=1,
        xmax=full_length,
        ymin=1,
        ymax=1.1),
    fill="white",
    color="black"
  ) +

  geom_rect(
    data=domain_df,
    aes(xmin=start,
        xmax=end,
        ymin=1,
        ymax=1.1,
        fill=Domain),
    color="black",
    alpha=0.8
  ) +

  geom_segment(
    data=protein_Gene,
    aes(x=Codon_Position,
        xend=Codon_Position,
        y=y,
        yend=yend),
    color="black"
  ) +

  geom_point(
    data=protein_Gene,
    aes(x=Codon_Position,
        y=yend,
        color=Sleep_Timing),
    size=3
  ) +

  geom_text_repel(
    data=protein_Gene,
    aes(x=Codon_Position,
        y=yend,
        label=mutation_label),
    size=3.2,
    fontface="bold",
    segment.color=NA,
    max.overlaps=10
  ) +

  scale_fill_manual(
    values=setNames(domain_df$fill,
                    domain_df$Domain)
  ) +

  scale_x_continuous(
    breaks=seq(0, full_length, 50),
    limits=c(0, full_length+10)
  ) +

  coord_cartesian(ylim=c(0.6,1.5)) +

  theme_minimal() +
  theme(
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid=element_blank(),
    plot.title=element_text(
      hjust=0.5,
      face="bold"
    )
  ) +

  labs(
    title="PER1 Amino Acid Mutations",
    x=NULL
  )

# ----------------------------------------------------------
# 6. Save
# ----------------------------------------------------------

OUT_DIR <- "figures/Result4/"
dir.create(OUT_DIR,
           recursive=TRUE,
           showWarnings=FALSE)

ggsave(
  file.path(
    OUT_DIR,
    "Figure4E_PER1_AminoAcid.pdf"
  ),
  p1,
  width=8,
  height=3
)
