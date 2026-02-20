############################################################
# Result6_Figure_ARNTL2_Mutation_Map.R
#
# Protein domain mutation map for NREM ratio–associated SNPs
#
# Input:
#   data/Result1/circadian_Gene_list.tsv
#   data/Result3/NREM_ratio_AminoAcid_mutation.tsv
#
# Output:
#   figures/Result3/ARNTL2_Mutation_Map.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

## =========================================================
## 1. Libraries
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

## =========================================================
## 2. Paths
## =========================================================

PROJECT_DIR <- getwd()

CHR_FILE <- file.path(PROJECT_DIR,
                      "data", "Result1",
                      "circadian_Gene_list.tsv")

MUT_FILE <- file.path(PROJECT_DIR,
                      "data", "Result3",
                      "NREM_ratio_AminoAcid_mutation.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures", "Result3")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 3. Load data
## =========================================================

circadian_chr <- fread(CHR_FILE)
mutation_df   <- fread(MUT_FILE)

mutation_df <- mutation_df[,-c(1:2)]

mutation_df$Position <- sapply(
  strsplit(mutation_df$Position, "_"),
  `[`, 2
)

mutation_df$Chr <- circadian_chr$Chr[
  match(mutation_df$Gene, circadian_chr$Gene_symbol)
]

mutation_df$Position <- paste0(
  "Chr", mutation_df$Chr, ":", mutation_df$Position
)

## Keep only non-synonymous
mutation_df <- mutation_df %>%
  filter(Mutation_Type == "Non-synonymous")

## =========================================================
## 4. Select gene of interest
## =========================================================

target_gene <- "ARNTL2"

protein_df <- mutation_df %>%
  filter(Gene == target_gene) %>%
  mutate(
    mutation_label = paste0(
      Original_Amino_Acid,
      Codon_Position,
      Mutated_Amino_Acid,
      " (", Position, ")"
    ),
    y = ifelse(Category == "long", 1.1, 1),
    yend = ifelse(
      Category == "long",
      runif(n(), 1.3, 1.38),
      runif(n(), 0.65, 0.75)
    )
  )

if (nrow(protein_df) == 0) {
  stop("No mutations found for ", target_gene)
}

full_length <- unique(protein_df$Full_Length)

## =========================================================
## 5. Domain annotation (edit if needed)
## =========================================================

domain_df <- data.frame(
  name  = c("bHLH","PAS","PAS","PAC"),
  start = c(107,178,357,432),
  end   = c(160,250,427,475),
  fill  = c("#fbf2c4","#74a892","#74a892","#003f5c")
)

## =========================================================
## 6. Plot
## =========================================================

p <- ggplot() +
  # Full protein backbone
  geom_rect(aes(xmin=1, xmax=full_length,
                ymin=1, ymax=1.1),
            fill="white", color="black") +

  # Domains
  geom_rect(data=domain_df,
            aes(xmin=start, xmax=end,
                ymin=1, ymax=1.1,
                fill=name),
            color="black", alpha=0.8) +

  # Mutation stems
  geom_segment(data=protein_df,
               aes(x=Codon_Position,
                   xend=Codon_Position,
                   y=y, yend=yend),
               color="black") +

  # Mutation points
  geom_point(data=protein_df,
             aes(x=Codon_Position,
                 y=yend,
                 color=Category),
             size=3) +

  # Labels
  geom_text_repel(
    data=protein_df,
    aes(x=Codon_Position,
        y=yend,
        label=mutation_label),
    direction="y",
    size=3.5,
    fontface="bold",
    box.padding=0.2,
    point.padding=0.2,
    segment.color=NA
  ) +

  scale_fill_manual(values=setNames(domain_df$fill,
                                    domain_df$name)) +

  scale_x_continuous(
    breaks=seq(0, full_length, by=100),
    limits=c(0, full_length+10)
  ) +

  coord_cartesian(ylim=c(0.6,1.5)) +

  theme_minimal() +
  theme(
    axis.text.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    axis.title.x = element_text(face="bold"),
    plot.title   = element_text(hjust=0.5,
                                face="bold")
  ) +

  labs(
    x="",
    title=paste(target_gene,
                "Mutation Map")
  )

## =========================================================
## 7. Save
## =========================================================

ggsave(file.path(OUT_DIR,
                 paste0(target_gene,
                        "_Mutation_Map.pdf")),
       p,
       width=7, height=3)

ggsave(file.path(OUT_DIR,
                 paste0(target_gene,
                        "_Mutation_Map.jpg")),
       p,
       width=7, height=3,
       dpi=300)

message("Mutation map saved.")
