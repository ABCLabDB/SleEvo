############################################################
# Fig3D_and_Fig4C_evolution_pressure.R
#
# Manuscript panels for evolutionary pressure (dN/dS / site models).
# Replaces the older Fig3_D / Fig4_C selection-signature scripts.
#
# ── Fig3D ─────────────────────────────────────────────────
#   Sleep duration + NREM ratio
#   Manuscript Fig3 panel D style: gene-wide ω vs site-class
#   positive selection (-log10 FDR, M1a vs M2a), plus supporting
#   ω distribution (P1) and significant-gene summary (P7).
#
# ── Fig4C ─────────────────────────────────────────────────
#   Sleep timing + Sleep frequency
#   Manuscript Fig4 panel C: evolutionary pressure for the two
#   categorical sleep phenotypes (same three panel types).
#
# Panel types (each for Fig3D and Fig4C):
#   P1  Gene-wide dN/dS (M0) violin
#   P2  Constraint vs site-class positive selection scatter
#   P7  % genes with site-level positive selection (bar)
#
# Input (data/Evolution_Pressure/):
#   Sleep_duration.tsv, NREM_ratio.tsv,
#   Sleep_timing.tsv, Sleep_frequency.tsv
#
# Output (Result/Evolution_Pressure/):
#   Fig3D_P1_*.jpg|.pdf , Fig3D_P2_*.jpg|.pdf , Fig3D_P7_*.jpg|.pdf
#   Fig4C_P1_*.jpg|.pdf , Fig4C_P2_*.jpg|.pdf , Fig4C_P7_*.jpg|.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()
DATA_DIR <- file.path(PROJECT_DIR, "data", "Evolution_Pressure")
OUT_DIR  <- file.path(PROJECT_DIR, "Result", "Evolution_Pressure")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DPI <- 300
FMT <- c("jpg", "pdf")

## =========================================================
## 2. Phenotype settings
## =========================================================

PHS <- c("Sleep_duration", "NREM_ratio", "Sleep_timing", "Sleep_frequency")
LAB <- c(
  Sleep_duration  = "Sleep duration",
  NREM_ratio      = "NREM ratio",
  Sleep_timing    = "Sleep timing",
  Sleep_frequency = "Sleep frequency"
)

## Fig3D = continuous sleep traits; Fig4C = categorical sleep traits
GRP <- list(
  Fig3D = c("Sleep_duration", "NREM_ratio"),
  Fig4C = c("Sleep_timing", "Sleep_frequency")
)
FIG_TITLE <- c(
  Fig3D = "Fig3D  |  Sleep duration & NREM ratio",
  Fig4C = "Fig4C  |  Sleep timing & Sleep frequency"
)

PCOL <- c(
  `Sleep duration`  = "#2C6E8F",
  `NREM ratio`      = "#5FA8AF",
  `Sleep timing`    = "#9B3D5B",
  `Sleep frequency` = "#D08A6B"
)

## highlighted genes (diamonds) — three per phenotype
HL <- list(
  Sleep_duration  = c("ATF4", "ADRB1", "MAGEL2"),
  NREM_ratio      = c("ARNTL2", "CPT1A", "PPARA"),
  Sleep_timing    = c("PER1", "CRY2", "CAVIN3"),
  Sleep_frequency = c("ATF5", "HDAC1", "NFIL3")
)
HLCOL <- "#C2185B"

SZ_PT_HL <- 6.5
SZ_LAB   <- 4
SZ_LAB_HL <- 4.5
SZ_ANN   <- 4
SZ_MED   <- 4.5

## =========================================================
## 3. Load TSVs
## =========================================================

D <- rbindlist(lapply(PHS, function(p) {
  f <- file.path(DATA_DIR, paste0(p, ".tsv"))
  if (!file.exists(f)) stop("Missing: ", f)
  d <- fread(f)
  d$ph <- p
  d
}), fill = TRUE)

D$Phenotype <- factor(LAB[D$ph], levels = unname(LAB))
D$omega_M0  <- as.numeric(D$omega_M0)
D$FDR_M1a_vs_M2a <- as.numeric(D$FDR_M1a_vs_M2a)
D$N_BEB_sites_P99_M1a_vs_M2a <- as.numeric(D$N_BEB_sites_P99_M1a_vs_M2a)

if (!"Significant_M1a_vs_M2a_FDR05" %in% names(D)) {
  D$Significant_M1a_vs_M2a_FDR05 <- ifelse(
    !is.na(D$FDR_M1a_vs_M2a) & D$FDR_M1a_vs_M2a < 0.05, "YES", "NO"
  )
}
if (!"Status_M1a_vs_M2a" %in% names(D)) {
  D$Status_M1a_vs_M2a <- ifelse(!is.na(D$FDR_M1a_vs_M2a), "OK", NA_character_)
}

D$sigA <- D$Significant_M1a_vs_M2a_FDR05 == "YES" &
  !is.na(D$Significant_M1a_vs_M2a_FDR05)
D$hl <- mapply(function(g, p) g %in% HL[[p]], D$Gene_symbol, D$ph)

N_PH <- tapply(D$Gene_symbol, D$ph, function(x) length(unique(x)))
message("Loaded genes: ",
        paste(sprintf("%s=%d", names(N_PH), N_PH), collapse = ", "))

nlog10 <- function(p, cap = 250) {
  v <- -log10(pmax(as.numeric(p), 1e-300))
  pmin(v, cap)
}

th <- theme_minimal(base_size = 14) + theme(
  plot.title    = element_text(face = "bold", size = 16, margin = margin(b = 4)),
  plot.subtitle = element_text(size = 12, colour = "grey35", margin = margin(b = 8)),
  axis.title    = element_text(size = 16),
  axis.text     = element_text(size = 13, colour = "grey20"),
  panel.grid    = element_blank(),
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background  = element_rect(fill = "white", colour = NA),
  axis.line     = element_line(colour = "grey35", linewidth = 0.4),
  legend.position = "none",
  plot.margin   = margin(10, 12, 10, 10)
)
th_leg <- th + theme(
  legend.position = "top",
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.4, "cm"),
  legend.margin = margin(0, 0, 2, 0),
  legend.box.margin = margin(0, 0, 2, 0)
)

save_fig <- function(plot, name, width = 8, height = 6) {
  for (fm in FMT) {
    f <- file.path(OUT_DIR, paste0(name, ".", fm))
    if (fm == "pdf") {
      ggsave(f, plot, width = width, height = height, device = cairo_pdf,
             bg = "white")
    } else {
      ggsave(f, plot, width = width, height = height, dpi = DPI,
             device = "jpeg", quality = 96, bg = "white")
    }
  }
  message("Saved: Result/Evolution_Pressure/", name, ".jpg|.pdf")
}

pick <- function(z, by, n) {
  z <- z[order(-z[[by]]), ]
  z <- z[!duplicated(z$Gene_symbol), ]
  head(z, n)
}

## =========================================================
## 4. Build P1 / P2 / P7 for Fig3D or Fig4C
## =========================================================

mk <- function(fig_id) {
  ## fig_id: "Fig3D" or "Fig4C"
  phs  <- GRP[[fig_id]]
  labs <- unname(LAB[phs])
  d <- as.data.frame(D[D$ph %in% phs, ])
  d$Phenotype <- factor(LAB[d$ph], levels = labs)
  cols <- PCOL[labs]
  tagp <- fig_id
  fig_lab <- FIG_TITLE[[fig_id]]
  nlab <- paste(sprintf("%s n=%d", labs, as.integer(N_PH[phs])),
                collapse = "; ")

  message("\n========== ", fig_lab, " ==========")

  ## ── P1. Gene-wide omega(M0) ─────────────────────────────────────
  med <- aggregate(omega_M0 ~ Phenotype, d, median, na.rm = TRUE)
  x_max <- max(d$omega_M0, na.rm = TRUE) * 1.18
  p1 <- ggplot(d, aes(omega_M0, Phenotype, fill = Phenotype)) +
    geom_violin(colour = NA, alpha = .6, scale = "width", width = .8) +
    geom_boxplot(width = .15, outlier.size = 1.4, outlier.colour = "grey35",
                 fill = "white", colour = "grey25", linewidth = .45) +
    geom_text(data = med, aes(label = sprintf("%.3f", omega_M0)),
              vjust = -1.5, size = SZ_MED, colour = "grey20",
              fontface = "bold") +
    geom_point(data = d[d$hl, ], aes(omega_M0, Phenotype),
               inherit.aes = FALSE, shape = 23, fill = "white",
               colour = HLCOL, stroke = 1.1, size = SZ_PT_HL) +
    geom_text_repel(
      data = d[d$hl, ],
      aes(omega_M0, Phenotype, label = Gene_symbol),
      inherit.aes = FALSE, size = SZ_LAB_HL, fontface = "bold.italic",
      colour = HLCOL, segment.colour = HLCOL, segment.size = .4,
      min.segment.length = 0, direction = "y", nudge_y = -.30,
      box.padding = .4, max.overlaps = Inf, seed = 11
    ) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = pretty(c(0, x_max), n = 5),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(title = paste0(fig_id, "  |  Gene-wide dN/dS (M0)"),
         x = "ω  (dN/dS)", y = NULL,
         subtitle = paste0(fig_lab, "; diamonds = highlighted genes; ", nlab)) +
    th
  save_fig(p1, paste0(tagp, "_P1_omega_distribution"), width = 7, height = 5)

  ## ── P2. Constraint vs site-class selection (main Fig3D-style panel)
  d2 <- d[!is.na(d$FDR_M1a_vs_M2a), ]
  d2$y <- nlog10(d2$FDR_M1a_vs_M2a, cap = 220)
  thr <- -log10(0.05)
  labg <- pick(d2[d2$y > thr & !d2$hl, , drop = FALSE], "y", 6)

  y_top <- max(268, max(d2$y, na.rm = TRUE) * 1.05)

  p2 <- ggplot(d2, aes(omega_M0, y)) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#F6E7E2", alpha = .45) +
    annotate("rect", xmin = -Inf, xmax = 1, ymin = thr, ymax = Inf,
             fill = "#C8E0EA", alpha = .28) +
    geom_hline(yintercept = thr, linetype = "22", colour = "grey40",
               linewidth = .5) +
    geom_vline(xintercept = 1, colour = "#B4442E", linewidth = .7) +
    geom_point(
      data = d2[!d2$hl, ],
      aes(fill = Phenotype, size = N_BEB_sites_P99_M1a_vs_M2a),
      shape = 21, colour = "white", stroke = .35, alpha = .5
    ) +
    geom_text_repel(
      data = labg, aes(label = Gene_symbol), size = SZ_LAB,
      fontface = "italic", max.overlaps = 30, min.segment.length = 0,
      segment.size = .3, segment.colour = "grey55",
      box.padding = .4, seed = 1
    ) +
    geom_point(
      data = d2[d2$hl, ], aes(omega_M0, y),
      inherit.aes = FALSE, shape = 23, fill = "white",
      colour = HLCOL, stroke = 1.1, size = SZ_PT_HL
    ) +
    geom_text_repel(
      data = d2[d2$hl, ],
      aes(omega_M0, y, label = Gene_symbol),
      inherit.aes = FALSE, size = SZ_LAB_HL, fontface = "bold.italic",
      colour = HLCOL, segment.colour = HLCOL, segment.size = .4,
      min.segment.length = 0, box.padding = .55, point.padding = .35,
      force = 18, force_pull = .4, max.overlaps = Inf, seed = 11
    ) +
    scale_fill_manual(values = cols) +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, .25)) +
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, 200, 50),
      expand = expansion(mult = c(0.03, 0.06))
    ) +
    annotate("text", x = 0.02, y = y_top * 0.96, hjust = 0, vjust = 1,
             size = SZ_ANN, colour = "#215C7A", lineheight = .95,
             fontface = "bold",
             label = "site-class\npositive selection") +
    annotate("text", x = 1.02, y = 8, hjust = 0, vjust = 0,
             size = SZ_ANN, colour = "#B4442E", fontface = "bold",
             label = "ω = 1") +
    labs(
      title = paste0(
        fig_id,
        "  |  Gene-wide constraint versus site-class positive selection"
      ),
      subtitle = paste0(
        fig_lab, ". Above the dashed line: M1a vs M2a FDR < 0.05. ",
        "Diamonds = highlighted genes. ", nlab
      ),
      x = "gene-wide ω  (M0)",
      y = expression(-log[10]*"(FDR)   M1a vs M2a")
    ) +
    th_leg +
    guides(fill = guide_legend(override.aes = list(size = SZ_PT_HL)))
  save_fig(p2, paste0(tagp, "_P2_quadrant_scatter"), width = 8, height = 6)

  ## ── P7. Significant gene counts ──────────────────────────────────
  has_m7 <- "Significant_M7_vs_M8_FDR05" %in% names(d) &&
    any(!is.na(d$Significant_M7_vs_M8_FDR05) &
          d$Significant_M7_vs_M8_FDR05 != "")

  b <- do.call(rbind, lapply(phs, function(p) {
    z <- d[d$ph == p, , drop = FALSE]
    rows <- list(
      data.frame(
        Phenotype = unname(LAB[p]),
        test = "M1a vs M2a",
        sig = sum(z$Significant_M1a_vs_M2a_FDR05 == "YES", na.rm = TRUE),
        run = sum(!is.na(z$FDR_M1a_vs_M2a)),
        stringsAsFactors = FALSE
      )
    )
    if (has_m7) {
      rows[[2]] <- data.frame(
        Phenotype = unname(LAB[p]),
        test = "M7 vs M8",
        sig = sum(z$Significant_M7_vs_M8_FDR05 == "YES", na.rm = TRUE),
        run = sum(z$Status_M7_vs_M8 == "OK", na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }))
  b$Phenotype <- factor(b$Phenotype, levels = labs)
  b$pct <- ifelse(b$run > 0, b$sig / b$run * 100, 0)

  fill_vals <- c(`M1a vs M2a` = "#2C6E8F", `M7 vs M8` = "#C9A227")
  sub7 <- if (has_m7) {
    paste0(fig_lab, "; FDR < 0.05; M7/M8 incomplete and prone to false positives")
  } else {
    paste0(fig_lab, "; FDR < 0.05 (M1a vs M2a); M7/M8 not present in input tables")
  }

  p7 <- ggplot(b, aes(Phenotype, pct, fill = test)) +
    geom_col(position = position_dodge(.72), width = .62,
             colour = "white", linewidth = .25) +
    geom_text(aes(label = sprintf("%d/%d", sig, run)),
              position = position_dodge(.72), vjust = -0.5,
              size = SZ_MED) +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25)) +
    labs(title = paste0(fig_id, "  |  Genes with site-level positive selection"),
         subtitle = sub7, x = NULL, y = "% of analysed genes") +
    th_leg
  save_fig(p7, paste0(tagp, "_P7_significant_counts"), width = 7, height = 5)
}

## =========================================================
## 5. Run Fig3D then Fig4C
## =========================================================

## Fig3D — Sleep duration & NREM ratio (manuscript Fig3 panel D)
mk("Fig3D")

## Fig4C — Sleep timing & Sleep frequency (manuscript Fig4 panel C)
mk("Fig4C")

message("\nDone. Outputs under Result/Evolution_Pressure/")
message("  Fig3D_* = Sleep duration + NREM ratio")
message("  Fig4C_* = Sleep timing + Sleep frequency")
