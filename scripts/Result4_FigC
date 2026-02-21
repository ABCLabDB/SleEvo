############################################################
# Result4_Fig4C.R
#
# Evolutionary regime plot:
#   dN/dS vs Tajima's D for sleep timing–associated genes
#
# Input:
#   data/Result4/Sleeptiming_Cochran_Result.tsv
#   data/Heatmap/Sleep_timing.tsv
#
# Output:
#   figures/Result4/Figure4C_Evolution_Quadrant.pdf
#
# Project root required:
#   Sleep_Evolution/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

PROJECT_DIR <- getwd()

## =========================================================
## 1. Paths
## =========================================================

COCHRAN_FILE <- file.path(PROJECT_DIR,
                          "data","Result4",
                          "Sleeptiming_Cochran_Result.tsv")

EVOLUTION_FILE <- file.path(PROJECT_DIR,
                            "data","Heatmap",
                            "Sleep_timing.tsv")

OUT_DIR <- file.path(PROJECT_DIR,
                     "figures","Result4")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. Load sleep timing–associated genes
## =========================================================

sleep_candidate <- fread(COCHRAN_FILE) |> as.data.frame()

sleep_candidate <- sleep_candidate %>%
  filter(P.cochran < 0.05)

sleep_genes <- sleep_candidate$gene.idx

## =========================================================
## 3. Load evolutionary scores
## =========================================================

evolution_df <- fread(EVOLUTION_FILE) |> as.data.frame()

evolution_df <- evolution_df %>%
  filter(Gene %in% sleep_genes)

## =========================================================
## 4. Categorize quadrants
## =========================================================

evolution_df <- evolution_df %>%
  mutate(
    dN_dS_Group = ifelse(dN_dS > 1,
                         "dN/dS > 1",
                         "dN/dS <= 1"),
    TajimasD_Group = ifelse(TajimasD > 0,
                            "Tajima's D > 0",
                            "Tajima's D <= 0")
  )

## =========================================================
## 5. Plot
## =========================================================

p <- ggplot(evolution_df,
            aes(x = dN_dS,
                y = TajimasD)) +

  ## Quadrant borders
  geom_rect(aes(xmin=-Inf, xmax=1,
                ymin=0, ymax=Inf),
            fill=NA, color="#ab5852", linewidth=1.2) +
  geom_rect(aes(xmin=1, xmax=Inf,
                ymin=0, ymax=Inf),
            fill=NA, color="#7593af", linewidth=1.2) +
  geom_rect(aes(xmin=-Inf, xmax=1,
                ymin=-Inf, ymax=0),
            fill=NA, color="#d69e49", linewidth=1.2) +
  geom_rect(aes(xmin=1, xmax=Inf,
                ymin=-Inf, ymax=0),
            fill=NA, color="#eadaa0", linewidth=1.2) +

  geom_point(aes(color=interaction(dN_dS_Group,
                                   TajimasD_Group)),
             size=5) +

  geom_text(aes(label=Gene),
            size=4,
            vjust=-0.5,
            check_overlap=TRUE) +

  geom_hline(yintercept=0,
             linetype="dashed",
             color="gray40") +
  geom_vline(xintercept=1,
             linetype="dashed",
             color="gray40") +

  scale_color_manual(
    values=c(
      "dN/dS > 1.Tajima's D > 0"   ="#7593af",
      "dN/dS > 1.Tajima's D <= 0" ="#eadaa0",
      "dN/dS <= 1.Tajima's D > 0" ="#ab5852",
      "dN/dS <= 1.Tajima's D <= 0"="#d69e49"
    ),
    name="Evolutionary regime"
  ) +

  theme_minimal(base_size=14) +
  theme(
    legend.position="right",
    axis.title=element_text(size=14)
  ) +

  labs(
    x="dN/dS",
    y="Tajima's D",
    title=NULL
  )

## =========================================================
## 6. Save
## =========================================================

ggsave(file.path(OUT_DIR,
                 "Figure4C_Evolution_Quadrant.pdf"),
       p,
       width=7,height=4)

ggsave(file.path(OUT_DIR,
                 "Figure4C_Evolution_Quadrant.jpg"),
       p,
       width=7,height=4,
       dpi=300)

message("Figure4C saved.")
