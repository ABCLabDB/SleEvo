############################################################
# Fig4_E_PER1_AminoAcid_mutation.R
#
# Figure 4E. PER1 protein lollipop mutation plot
#   - UniProt features: PAS (x2), PAC, NES (x3), NLS, LXXLL, CRY binding
#   - deterministic tier stacking (reproducible)
#   - duplicate positions collapsed; point size = count
#   - Daytime (orange) above backbone; Night (green) below
#
# Input:
#   - data/Result4/Sleeptiming_AA_Mutation.tsv
#   - data/Result4/Sleeptiming_Cochran_Result.tsv
#
# Output:
#   - Result/Fig4/Fig4_E_PER1_AminoAcid_mutation.jpg
#   - Result/Fig4/Fig4_E_PER1_AminoAcid_mutation.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

MUTATION_FILE <- file.path(
  PROJECT_DIR, "data", "Result4",
  "Sleeptiming_AA_Mutation.tsv"
)
COCHRAN_FILE <- file.path(
  PROJECT_DIR, "data", "Result4",
  "Sleeptiming_Cochran_Result.tsv"
)

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Data
## =========================================================

mutation_results <- fread(MUTATION_FILE) |> as.data.frame()

## drop unnamed index column if present
if (names(mutation_results)[1] %in% c("V1", "")) {
  mutation_results <- mutation_results[, -1, drop = FALSE]
}

## genomic SNP id: Chr<Chromosome>:<CDS coord>
mutation_results$Position <- paste0(
  "Chr", mutation_results$Chromosome, ":",
  sub("^CDS_", "", mutation_results$Position)
)

timing_gene <- fread(COCHRAN_FILE) |> as.data.frame()
timing_gene <- timing_gene[timing_gene$P.cochran < 0.05, ]

mutation_results <- mutation_results[
  mutation_results$Gene %in% timing_gene$gene.idx, ]
mutation_results <- mutation_results[
  mutation_results$Mutation_Type == "Non-synonymous", ]


## =========================================================
## 3. Settings
## =========================================================

target_gene <- "PER1"
full_length <- 1290   # UniProt O15534

## Fig4 E: night = green (below), daytime = orange (above)
group_colors <- c(
  "Sleep at night"   = "#2c6e49",
  "Sleep at daytime" = "#d68c45"
)

feature_df <- data.frame(
  Feature = c("PAS", "PAS", "PAC",
              "Nuclear export signal", "Nuclear export signal",
              "Nuclear export signal",
              "Nuclear localization signal", "LXXLL",
              "CRY binding domain"),
  start   = c(208, 348, 422, 138, 489, 982, 827, 1043, 1149),
  end     = c(275, 414, 465, 147, 498, 989, 843, 1047, 1290)
)

feature_levels <- c(
  "PAS", "PAC", "Nuclear export signal",
  "Nuclear localization signal", "LXXLL", "CRY binding domain"
)
feature_colors <- c(
  "PAS"                         = "#3d6f9e",
  "PAC"                         = "#8fb8d8",
  "Nuclear export signal"       = "#c0603f",
  "Nuclear localization signal" = "#e0a45e",
  "LXXLL"                       = "#7b5aa6",
  "CRY binding domain"          = "#74a892"
)
feature_df$Feature <- factor(feature_df$Feature, levels = feature_levels)

## tiny motifs -> enforce minimum drawn width
min_w <- full_length * 0.006
feature_df <- feature_df %>%
  mutate(
    mid  = (start + end) / 2,
    w    = pmax(end - start, min_w),
    xmin = mid - w / 2,
    xmax = mid + w / 2
  )


## =========================================================
## 4. Collapse duplicates, assign tiers
## =========================================================

pg <- mutation_results %>%
  filter(Gene == target_gene) %>%
  mutate(
    Group = ifelse(Sleep_Timing == "Sleep at night",
                   "Sleep at night", "Sleep at daytime"),
    mutation_label = paste0(Original_Amino_Acid, Codon_Position,
                            Mutated_Amino_Acid)
  ) %>%
  group_by(Codon_Position, Group, mutation_label, Position) %>%
  summarise(n = sum(Count, na.rm = TRUE), .groups = "drop")

if (is.na(full_length) || length(full_length) == 0) {
  full_length <- as.numeric(
    unique(mutation_results$Full_Length[
      mutation_results$Gene == target_gene])
  )[1]
}

min_gap <- full_length * 0.045
n_tier  <- 3

assign_tier <- function(pos) {
  ord  <- order(pos)
  tier <- integer(length(pos))
  last <- rep(-Inf, n_tier)
  for (k in ord) {
    t <- 1
    while (t < n_tier && pos[k] - last[t] < min_gap) t <- t + 1
    tier[k] <- t
    last[t] <- pos[k]
  }
  tier
}

pg <- pg %>%
  group_by(Group) %>%
  mutate(tier = assign_tier(Codon_Position)) %>%
  ungroup()

## geometry: daytime UP (orange), night DOWN (green) — Fig4 E
up        <- pg$Group == "Sleep at daytime"
step      <- 0.115
pg$y_base <- ifelse(up, 1.02, 0.98)
pg$yend   <- ifelse(up, 1.02 + step * pg$tier, 0.98 - step * pg$tier)

ymax <- max(pg$yend) + 0.20
ymin <- min(pg$yend) - 0.20

## highlighted mutations (Fig4 E labels)
highlight_snps <- c("Chr17:8144782", "Chr17:8142834")  # T1153A, A962P
pg$is_hl <- pg$Position %in% highlight_snps
pg$label_show <- ifelse(
  pg$is_hl,
  paste0(pg$mutation_label, " (", pg$Position, ")"),
  pg$mutation_label
)


## =========================================================
## 5. Plot
## =========================================================

p1 <- ggplot() +

  ## protein backbone
  geom_rect(aes(xmin = 1, xmax = full_length, ymin = 0.98, ymax = 1.02),
            fill = "#eceef0") +

  ## UniProt features
  geom_rect(data = feature_df,
            aes(xmin = xmin, xmax = xmax, ymin = 0.973, ymax = 1.027,
                fill = Feature)) +
  scale_fill_manual(values = feature_colors, name = "Domain", drop = FALSE) +

  ## stems
  geom_segment(data = subset(pg, !is_hl),
               aes(x = Codon_Position, xend = Codon_Position,
                   y = y_base, yend = yend),
               colour = "#b8bcc0", linewidth = 0.3) +
  geom_segment(data = subset(pg, is_hl),
               aes(x = Codon_Position, xend = Codon_Position,
                   y = y_base, yend = yend),
               colour = "#6b7280", linewidth = 0.55) +

  ## heads
  geom_point(data = pg, aes(x = Codon_Position, y = yend, size = n),
             colour = "white") +
  geom_point(data = pg, aes(x = Codon_Position, y = yend,
                            colour = Group, size = n * 0.72)) +
  geom_point(data = subset(pg, is_hl),
             aes(x = Codon_Position, y = yend, size = n * 0.72),
             shape = 21, fill = NA, colour = "#1f2937", stroke = 0.7) +
  scale_colour_manual(values = group_colors, name = "Group") +
  scale_size_continuous(range = c(1.8, 5), guide = "none") +

  ## labels
  geom_text_repel(
    data = subset(pg, !is_hl),
    aes(x = Codon_Position, y = yend, label = label_show),
    size = 2.7, colour = "#33383d",
    direction = "y",
    nudge_y = ifelse(subset(pg, !is_hl)$Group == "Sleep at daytime",
                     0.05, -0.05),
    segment.colour = NA, box.padding = 0.18,
    max.overlaps = Inf, seed = 1
  ) +
  geom_text_repel(
    data = subset(pg, is_hl),
    aes(x = Codon_Position, y = yend, label = label_show),
    size = 3.1, fontface = "bold", colour = "#1f2937",
    direction = "y",
    nudge_y = ifelse(subset(pg, is_hl)$Group == "Sleep at daytime",
                     0.13, -0.13),
    segment.colour = "#6b7280", segment.size = 0.3,
    min.segment.length = 0, box.padding = 0.4,
    max.overlaps = Inf, seed = 1
  ) +

  scale_x_continuous(breaks = seq(0, full_length, by = 200),
                     limits = c(-full_length * 0.01, full_length * 1.02),
                     expand = expansion(mult = 0.005)) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(x = "Amino acid position", y = NULL,
       title = "Non-synonymous mutation in PER1") +

  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.line.y      = element_blank(),
    axis.line.x      = element_line(colour = "#c9cdd1", linewidth = 0.35),
    axis.ticks.x     = element_line(colour = "#c9cdd1", linewidth = 0.35),
    axis.text.x      = element_text(size = 7.5, colour = "#70767c"),
    axis.title.x     = element_text(size = 9, colour = "#4b5157",
                                    margin = margin(t = 6)),
    plot.title       = element_text(size = 12, face = "bold",
                                    colour = "#1f2937"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 7.5),
    legend.key.size  = unit(0.35, "cm"),
    legend.box       = "vertical",
    legend.margin    = margin(t = -2),
    plot.margin      = margin(6, 12, 4, 6)
  ) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         fill   = guide_legend(order = 2, nrow = 2))


## =========================================================
## 6. Save
## =========================================================

ggsave(file.path(OUT_DIR, "Fig4_E_PER1_AminoAcid_mutation.pdf"),
       p1, width = 200, height = 95, units = "mm")
ggsave(file.path(OUT_DIR, "Fig4_E_PER1_AminoAcid_mutation.jpg"),
       p1, width = 200, height = 95, units = "mm", dpi = 600, bg = "white")

message("PER1 mutations plotted: ", nrow(pg),
        " | highlighted: ", sum(pg$is_hl))
message("Saved: Result/Fig4/Fig4_E_PER1_AminoAcid_mutation.jpg|.pdf")
