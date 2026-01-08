############################################################
# Result2_AA_PML_visualization.R
#
# Representative amino-acid mutation visualization for PML
# (Figure 2I / 2J)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


## =========================================================
## Paths
## =========================================================

PROJECT_DIR <- getwd()

AA_FILE <- file.path(PROJECT_DIR,
                     "data/Result2/Total_sleep_time_AA_Mutation_Data.tsv")

FIG_DIR <- file.path(PROJECT_DIR,
                     "figures/Result2/AminoAcid")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## Load data
## =========================================================

mutation_results <- fread(AA_FILE) |> as.data.frame()

plot_df <- mutation_results |>
  filter(Gene == "PML") |>
  mutate(
    label = paste0(Original_Amino_Acid,
                   Codon_Position,
                   Mutated_Amino_Acid),
    y = 1,
    yend = runif(n(), 1.2, 1.4)
  )

full_length <- unique(plot_df$Full_Length)

domain_df <- data.frame(
  start = c(476, 556),
  end = c(490, 562)
)


## =========================================================
## Plot
## =========================================================

p <- ggplot() +
  geom_rect(aes(xmin = 1, xmax = full_length,
                ymin = 0.9, ymax = 1.1),
            fill = "white", color = "black") +
  geom_rect(data = domain_df,
            aes(xmin = start, xmax = end,
                ymin = 0.9, ymax = 1.1),
            fill = "#003f5c", alpha = 0.8) +
  geom_segment(data = plot_df,
               aes(x = Codon_Position,
                   xend = Codon_Position,
                   y = y, yend = yend)) +
  geom_point(data = plot_df,
             aes(x = Codon_Position, y = yend),
             size = 3) +
  geom_text_repel(data = plot_df,
                  aes(x = Codon_Position,
                      y = yend,
                      label = label),
                  size = 3) +
  theme_minimal() +
  labs(title = "PML amino-acid mutations",
       x = "Codon position",
       y = NULL)

ggsave(file.path(FIG_DIR, "PML_Mutation.pdf"),
       p, width = 7, height = 3, dpi = 300)
