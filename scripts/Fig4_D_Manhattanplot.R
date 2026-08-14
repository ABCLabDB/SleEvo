############################################################
# Fig4_D_Manhattanplot.R
#
# Figure 4D. Sleep-timing Manhattan plot (publication style)
#   Visual hierarchy:
#     non-significant (pale grey)
#     -> significant (muted alternating blues)
#     -> target SNPs (dark accent + leader-line labels)
#
# Input:
#   - data/Result4/Sleeptiming_Manhattan_Dataset.tsv
#
# Output:
#   - Result/Fig4/Fig4_D_Manhattanplot.jpg
#   - Result/Fig4/Fig4_D_Manhattanplot.pdf
#
# Required working directory: repository root (SleEvo)
############################################################


library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

MANHATTAN_FILE <- file.path(
  PROJECT_DIR, "data", "Result4",
  "Sleeptiming_Manhattan_Dataset.tsv"
)

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Data
## =========================================================

TS_site1 <- fread(MANHATTAN_FILE)
## Expected columns: SNP, CHR, BP, P, Gene
setnames(TS_site1, c("SNP", "CHR", "BP", "P", "Gene"))
TS_site1[, CHR := as.numeric(CHR)]
TS_site1 <- TS_site1[!is.na(CHR) & !is.na(P)]

## cumulative x positions
don <- TS_site1 %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  left_join(TS_site1, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot,
         logP  = -log10(P))

axisdf <- don %>%
  group_by(CHR) %>%
  summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")


## =========================================================
## 3. Styling
## =========================================================

threshold   <- -log10(0.05)
target_snps <- c("Chr17:8144782", "Chr19:49932821", "Chr11:6320204")

## muted alternating palette (significant points only)
chr_levels <- sort(unique(don$CHR))
chr_colors <- setNames(
  rep(c("#2c5d7c", "#7fb2d1"), length.out = length(chr_levels)),
  chr_levels
)

grey_pt <- "#d9dcdf"   # below threshold
accent  <- "#c0392b"   # target SNPs
ymax    <- max(don$logP, na.rm = TRUE)

hl <- subset(don, SNP %in% target_snps)


## =========================================================
## 4. Plot
## =========================================================

P1 <- ggplot(don, aes(x = BPcum, y = logP)) +

  ## significance line first, so points sit on top of it
  geom_hline(yintercept = threshold, linetype = "dashed",
             colour = "grey55", linewidth = 0.35) +

  ## 1. everything, pale grey
  geom_point(colour = grey_pt, size = 0.9, shape = 16) +

  ## 2. significant points, muted alternating colour
  geom_point(data = subset(don, logP > threshold),
             aes(colour = factor(CHR)), size = 1.2, shape = 16) +
  scale_colour_manual(values = chr_colors, guide = "none") +

  ## 3. target SNPs — white halo then accent fill
  geom_point(data = hl, colour = "white", size = 3.4, shape = 16) +
  geom_point(data = hl, colour = accent,  size = 2.2, shape = 16) +

  ## labels with thin leader lines
  geom_text_repel(data = hl, aes(label = SNP),
                  size = 3.1, colour = "grey15",
                  nudge_y = 0.9, direction = "y",
                  segment.colour = "grey55", segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.45,
                  max.overlaps = Inf, seed = 1) +

  ## axes
  scale_x_continuous(breaks = axisdf$center, labels = axisdf$CHR,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, ymax + 1.6),
                     breaks = pretty_breaks(4),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +

  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    legend.position   = "none",
    panel.grid        = element_blank(),
    axis.line.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.line.y       = element_line(colour = "grey30", linewidth = 0.35),
    axis.ticks.y      = element_line(colour = "grey30", linewidth = 0.35),
    axis.text.x       = element_text(size = 7,  colour = "grey25",
                                     margin = margin(t = 2)),
    axis.text.y       = element_text(size = 8,  colour = "grey25"),
    axis.title.x      = element_text(size = 9,  colour = "grey15",
                                     margin = margin(t = 6)),
    axis.title.y      = element_text(size = 9,  colour = "grey15",
                                     margin = margin(r = 6)),
    plot.margin       = margin(6, 10, 4, 6)
  )


## =========================================================
## 5. Export
## =========================================================

out_pdf <- file.path(OUT_DIR, "Fig4_D_Manhattanplot.pdf")
out_jpg <- file.path(OUT_DIR, "Fig4_D_Manhattanplot.jpg")

## Prefer cairo_pdf when available; fall back to default pdf device.
pdf_ok <- FALSE
tryCatch({
  ggsave(out_pdf, P1, width = 180, height = 62, units = "mm",
         device = cairo_pdf)
  pdf_ok <- TRUE
}, error = function(e) {
  message("cairo_pdf unavailable; using default pdf device.")
})
if (!pdf_ok) {
  ggsave(out_pdf, P1, width = 180, height = 62, units = "mm")
}

ggsave(out_jpg, P1, width = 180, height = 62, units = "mm", dpi = 600)

message("Highlighted SNPs found: ", nrow(hl), " / ", length(target_snps))
message("Saved: Result/Fig4/Fig4_D_Manhattanplot.jpg|.pdf")

## ---------------------------------------------------------------------------
##  Tuning notes
##  - Point sizes are small on purpose: at 180 mm wide they render crisply and
##    dense regions stay readable. Bump to 1.1 / 1.5 if the plot looks sparse.
##  - accent  : colour of the target SNPs.  chr_colors : significant points.
##  - The white halo under each target point keeps it visible inside dense
##    clusters without adding a black outline.
## ---------------------------------------------------------------------------
