############################################################
# Fig4_H_NFIL3_sleepfrequency_indel.R
#
# Figure 4H. NFIL3 indel schematic (Panel H style)
#   Grey gene bar + teal (shorter / human ref) vs red (longer / indel)
#   allele boxes, clustered by genomic proximity.
#
# Input:
#   - data/Result5/All_Indel_Split.tsv
#   - data/Result1/circadian_Gene_list.tsv
#   - data/Circadian_gene_Nucleotide_Matrix/NFIL3.tsv  (preferred)
#     OR data/nucleotideDF/NFIL3.tsv
#     OR fallback: data/Result3/NREM_ratio_Manhattan.tsv (NFIL3 BP map)
#
# Output:
#   - Result/Fig4/Fig4_H_NFIL3_sleepfrequency_indel.jpg
#   - Result/Fig4/Fig4_H_NFIL3_sleepfrequency_indel.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()
gene <- "NFIL3"

INDEL_FILE <- file.path(PROJECT_DIR, "data", "Result5",
                        "All_Indel_Split.tsv")
CHR_FILE   <- file.path(PROJECT_DIR, "data", "Result1",
                        "circadian_Gene_list.tsv")
SEQ_CANDIDATES <- c(
  file.path(PROJECT_DIR, "data", "Circadian_gene_Nucleotide_Matrix",
            paste0(gene, ".tsv")),
  file.path(PROJECT_DIR, "data", "nucleotideDF",
            paste0(gene, ".tsv"))
)
NREM_FILE <- file.path(PROJECT_DIR, "data", "Result3",
                       "NREM_ratio_Manhattan.tsv")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Alignment POS -> genomic label (TS_site2)
## =========================================================

circadian_Chr <- fread(CHR_FILE) |> as.data.frame()
chr_num <- as.character(
  circadian_Chr$Chr[match(gene, circadian_Chr$Gene_symbol)]
)
if (is.na(chr_num) || !nzchar(chr_num)) chr_num <- "9"

seq_file <- SEQ_CANDIDATES[file.exists(SEQ_CANDIDATES)][1]

if (!is.na(seq_file)) {
  message("POS map from nucleotide matrix: ", seq_file)
  aln_df <- fread(seq_file) |> as.data.frame()
  rownames(aln_df) <- aln_df[, 1]
  aln_df <- aln_df[, -1, drop = FALSE]
  aln_df1 <- aln_df[, grepl("V|CDS", colnames(aln_df)), drop = FALSE]

  TS_site2 <- data.frame(
    V1 = seq_len(ncol(aln_df1)),
    V2 = colnames(aln_df1),
    V3 = gene,
    stringsAsFactors = FALSE
  )
  TS_site2$Chr <- chr_num
  TS_site2$V2 <- ifelse(
    grepl("CDS", TS_site2$V2),
    paste0("Chr", TS_site2$Chr, ":", TS_site2$V2),
    TS_site2$V2
  )
  TS_site2$V2 <- ifelse(
    grepl("CDS_", TS_site2$V2),
    gsub("CDS_", "", TS_site2$V2),
    TS_site2$V2
  )
} else {
  ## Fallback: NREM Manhattan BP (= alignment index) -> Chr9:genomic
  message("NFIL3 nucleotide matrix missing; ",
          "building POS map from NREM_ratio_Manhattan.tsv")
  nrem <- fread(NREM_FILE) |> as.data.frame()
  nrem <- nrem[nrem$Gene == gene, ]
  if (nrow(nrem) == 0) {
    stop("No NFIL3 rows in NREM_ratio_Manhattan.tsv and no matrix file")
  }
  bp  <- as.integer(nrem$BP)
  gpos <- as.integer(sub(".*:", "", nrem$SNP))
  ord <- order(bp)
  bp <- bp[ord]; gpos <- gpos[ord]

  ## fill every alignment index between min/max by step-wise
  ## genomic advance when BP increases with gpos
  max_pos <- max(bp, na.rm = TRUE)
  lab <- rep(NA_character_, max_pos)
  for (i in seq_along(bp)) {
    lab[bp[i]] <- paste0("Chr", chr_num, ":", gpos[i])
  }
  ## interpolate / carry within segments where gpos rises with bp
  for (i in seq_len(length(bp) - 1)) {
    b0 <- bp[i]; b1 <- bp[i + 1]
    g0 <- gpos[i]; g1 <- gpos[i + 1]
    if (b1 <= b0) next
    if (g1 >= g0 && (b1 - b0) == (g1 - g0)) {
      ## contiguous CDS stretch
      for (b in b0:b1) {
        lab[b] <- paste0("Chr", chr_num, ":", g0 + (b - b0))
      }
    } else {
      ## gap-rich stretch: keep known ends; mark interiors as V#
      for (b in (b0 + 1):(b1 - 1)) {
        if (is.na(lab[b])) lab[b] <- paste0("V", b)
      }
    }
  }
  for (b in which(is.na(lab))) lab[b] <- paste0("V", b)

  TS_site2 <- data.frame(
    V1 = seq_along(lab),
    V2 = lab,
    V3 = gene,
    Chr = chr_num,
    stringsAsFactors = FALSE
  )
}


## =========================================================
## 3. Load indel data (non-SNV)
## =========================================================

indel_DF <- fread(INDEL_FILE) |> as.data.frame()
indel_DF <- indel_DF[indel_DF$Variant_Type != "SNV", ]
target_gene_Indel <- indel_DF[indel_DF$Gene %in% gene, ]

if (nrow(target_gene_Indel) == 0) {
  stop("No non-SNV indels found for ", gene)
}


## =========================================================
## 4. Highest-AF ALT per POS
## =========================================================

tg <- target_gene_Indel %>%
  group_by(POS) %>%
  slice_max(AF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(POS)

tg$gpos_label <- TS_site2$V2[match(tg$POS, TS_site2$V1)]

## if still V# / missing, snap to nearest Chr*:coord label by POS
chr_map <- TS_site2 %>%
  filter(grepl("^Chr", V2)) %>%
  transmute(POS = V1, lab = V2, g = as.numeric(sub(".*:", "", V2)))

snap_label <- function(pos) {
  if (nrow(chr_map) == 0) return(paste0("V", pos))
  i <- which.min(abs(chr_map$POS - pos))
  chr_map$lab[i]
}

need_snap <- is.na(tg$gpos_label) | grepl("^V", tg$gpos_label)
if (any(need_snap)) {
  tg$gpos_label[need_snap] <- vapply(
    tg$POS[need_snap], snap_label, character(1)
  )
}

tg$gpos <- suppressWarnings(as.numeric(sub(".*:", "", tg$gpos_label)))
miss <- is.na(tg$gpos)
if (any(miss)) {
  tg$gpos[miss] <- as.numeric(tg$POS[miss])
  message("Note: ", sum(miss),
          " site(s) lack CDS genomic labels (using POS).")
}

message("Indel sites kept:")
print(as.data.frame(tg[, c("POS", "REF", "ALT", "AF", "gpos_label")]))


## =========================================================
## 5. Clusters / slots / direction / alleles
## =========================================================

tg <- tg %>% arrange(gpos)
tg$cluster <- cumsum(c(0, diff(tg$gpos) > 50)) + 1

tg <- tg %>%
  group_by(cluster) %>%
  arrange(gpos, .by_group = TRUE) %>%
  mutate(
    x   = row_number(),
    dir = rep(c("up", "down"), length.out = n())
  ) %>%
  ungroup() %>%
  mutate(
    short_allele = ifelse(nchar(REF) <= nchar(ALT), REF, ALT),
    long_allele  = ifelse(nchar(REF) >  nchar(ALT), REF, ALT),
    y_box = ifelse(dir == "up",  0.85, -0.85),
    y_end = ifelse(dir == "up",  1.00, -1.00),
    y_lab = ifelse(dir == "up",  1.15, -1.15)
  )


## =========================================================
## 6. Helper data frames
## =========================================================

col_ref <- "#5FC9B7"
col_alt <- "#DE5E6E"
leg_ref <- "Human reference (Sleep once)"
leg_alt <- "Indel (Sleep more than twice)"

cls    <- sort(unique(tg$cluster))
bar_df <- tg %>%
  group_by(cluster) %>%
  summarise(x1 = 0.6, x2 = max(x) + 0.4, .groups = "drop")
gene_df   <- bar_df %>% filter(cluster == cls[length(cls)])
legend_df <- data.frame(
  type = factor(c(leg_ref, leg_alt), levels = c(leg_ref, leg_alt)),
  cluster = cls[1],
  xmin = 1, xmax = 1.01,
  ymin = 1.55, ymax = 1.56
)


## =========================================================
## 7. Plot
## =========================================================

p1 <- ggplot() +
  geom_segment(data = bar_df,
               aes(x = x1, xend = x2, y = 0, yend = 0),
               linewidth = 7, colour = "grey90", lineend = "round") +
  geom_text(data = gene_df,
            aes(x = x2 + 0.35, y = 0, label = gene),
            hjust = 0, fontface = "italic", size = 4.2) +
  geom_segment(data = tg,
               aes(x = x, xend = x, y = 0, yend = y_end),
               linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_label(data = tg, aes(x = x, y = y_box, label = long_allele),
             fill = col_alt, colour = "white", fontface = "bold",
             size = 3, linewidth = 0, label.r = unit(0, "pt")) +
  geom_label(data = tg, aes(x = x, y = 0, label = short_allele),
             fill = col_ref, colour = "black", fontface = "bold",
             size = 3, linewidth = 0, label.r = unit(0, "pt")) +
  geom_text(data = tg, aes(x = x, y = y_lab, label = gpos_label),
            size = 3) +
  geom_rect(data = legend_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = type),
            alpha = 0, colour = NA) +
  scale_fill_manual(
    values = setNames(c(col_ref, col_alt), c(leg_ref, leg_alt)),
    breaks = c(leg_ref, leg_alt),
    name = NULL,
    guide = guide_legend(
      override.aes = list(alpha = 1, colour = NA)
    )
  ) +
  facet_wrap(~ cluster, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(add = 0.9)) +
  coord_cartesian(ylim = c(-1.8, 1.8), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text      = element_blank(),
    panel.grid      = element_blank(),
    panel.spacing.x = unit(2, "lines"),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "top",
    plot.margin     = margin(12, 60, 12, 12)
  )


## =========================================================
## 8. Save
## =========================================================

ggsave(file.path(OUT_DIR, "Fig4_H_NFIL3_sleepfrequency_indel.pdf"),
       plot = p1, width = 11, height = 3.4)
ggsave(file.path(OUT_DIR, "Fig4_H_NFIL3_sleepfrequency_indel.jpg"),
       plot = p1, width = 11, height = 3.4, dpi = 300, bg = "white")

message("Saved: Result/Fig4/Fig4_H_NFIL3_sleepfrequency_indel.jpg|.pdf")
