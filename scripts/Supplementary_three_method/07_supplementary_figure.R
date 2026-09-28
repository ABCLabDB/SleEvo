############################################################
## Supplementary Figure 10 — 3기준 재현성 (A/B/C/D)
## 규격: Nature 더블컬럼 180 mm × 168 mm, 600 dpi
##  - 패널 태그는 패널 영역 바깥(좌상단 여백)에 배치
##  - 최소 글자 5 pt 이상, sans-serif
############################################################
.libPaths(c("/disk2/bijsy/Test/Evolution_pressure/0.Code/Rlib", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(grid); library(ggrepel)})

MP  <- "/disk2/bijsy/Test/Sleep_association_Test/Result/Maximum_Parsimony_for_Cladistics"
OUT <- file.path(MP, "3method_comparison")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")

A  <- fread(file.path(TAB, "SupplementaryTable_three_method_comparison.tsv"))
S  <- fread(file.path(TAB, "SupplementaryTable_reproducibility_by_phenotype.tsv"))
PHL  <- c("NREM ratio","Sleep duration","Sleep timing","Sleep frequency")
PHL2 <- c("NREM\nratio","Sleep\nduration","Sleep\ntiming","Sleep\nfrequency")
A[, Phenotype := factor(Phenotype, levels = PHL)]
S[, Phenotype := factor(Phenotype, levels = PHL)]

## ── 팔레트 ────────────────────────────────────────────────
INK   <- "#33393F"; INK2 <- "#6E767E"; GRID <- "#E8EBED"
MCOL  <- c(NJ = "#7FA9C9", MP = "#DFA475", ML = "#93B892")
SUPC  <- c(`3/3` = "#2F5F7D", `2/3` = "#7FA9C9", `1/3` = "#C7DCE8", `0/3` = "#F1F4F6")
SEQ   <- c("#F7FAFB","#D7E6EE","#A8C7DA","#6E9DBC","#2F5F7D")
POS   <- "#3C6F94"; NEG <- "#F1F4F6"

## ── 공통 테마 ─────────────────────────────────────────────
base <- theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    text            = element_text(colour = INK),
    plot.title      = element_text(size = 7.8, face = "bold", colour = INK,
                                   hjust = 0, margin = margin(b = 1, l = 14)),
    plot.subtitle   = element_text(size = 5.8, colour = INK2, hjust = 0,
                                   margin = margin(b = 5, l = 14), lineheight = 1.15),
    plot.title.position = "plot",
    plot.tag        = element_text(size = 10, face = "bold", colour = INK),
    plot.tag.position = c(0, 1),
    axis.title      = element_text(size = 6.4, colour = INK2),
    axis.text       = element_text(size = 6.2, colour = INK),
    axis.text.x     = element_text(size = 6.2, colour = INK, lineheight = .95),
    axis.ticks      = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = GRID, linewidth = .3),
    legend.title    = element_text(size = 6, colour = INK2),
    legend.text     = element_text(size = 5.9, colour = INK),
    legend.key.size = unit(6, "pt"),
    legend.margin   = margin(0,0,0,0),
    legend.box.spacing = unit(2, "pt"),
    plot.margin     = margin(13, 4, 2, 2))

flat <- base + theme(panel.grid.major.y = element_blank())

## ══ A ═════════════════════════════════════════════════════
Ad <- melt(S[, .(Phenotype, NJ = NJ_significant, MP = MP_significant, ML = ML_significant)],
           id.vars = "Phenotype", variable.name = "Method", value.name = "n")
pA <- ggplot(Ad, aes(Phenotype, n, fill = Method)) +
  geom_col(position = position_dodge(.76), width = .66,
           colour = NA, linewidth = 0, linejoin = "round") +
  geom_text(aes(label = n), position = position_dodge(.76), vjust = -0.55,
            size = 1.95, colour = INK) +
  scale_fill_manual(values = MCOL, name = NULL) +
  scale_x_discrete(labels = PHL2) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40),
                     expand = expansion(mult = c(0, .02))) +
  labs(title = "Significant genes under each criterion",
       subtitle = "143 genes tested per phenotype", x = NULL, y = "Genes") +
  base +
  theme(legend.position = "inside",
        legend.position.inside = c(.99, 1.02), legend.justification = c(1, 1),
        legend.direction = "vertical",
        legend.background = element_blank(), legend.key = element_blank())

## ══ B ═════════════════════════════════════════════════════
LEV <- c("0/3","1/3","2/3","3/3")
Bd  <- A[, .N, by = .(Phenotype, Support)]
Bd  <- Bd[CJ(Phenotype = factor(PHL, levels = PHL), Support = LEV),
          on = .(Phenotype, Support)][is.na(N), N := 0L]
Bd[, Support := factor(Support, levels = LEV)]
setorder(Bd, Phenotype, -Support)                       # 3/3 이 아래
Bd[, `:=`(top = cumsum(N), mid = cumsum(N) - N/2), by = Phenotype]
Bd[, inside := N >= 10]

pB <- ggplot(Bd, aes(Phenotype, N, fill = Support)) +
  geom_col(width = .62, colour = "white", linewidth = .35) +
  geom_text(data = Bd[inside == TRUE],
            aes(y = mid, label = N, colour = Support %in% c("3/3","2/3")),
            size = 1.95, fontface = "bold", show.legend = FALSE) +
  geom_text_repel(data = Bd[inside == FALSE & N > 0],
                  aes(y = mid, label = N), colour = INK2,
                  size = 1.95, fontface = "bold",
                  nudge_x = .40, direction = "y", hjust = 0,
                  min.segment.length = 0, segment.size = .2,
                  segment.colour = "grey70", box.padding = .08,
                  point.padding = 0, seed = 319, show.legend = FALSE) +
  scale_fill_manual(values = SUPC, name = "Criteria\nsignificant",
                    guide = guide_legend(reverse = TRUE, byrow = TRUE,
                                         override.aes = list(colour = NA, linewidth = 0))) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK2)) +
  scale_x_discrete(labels = PHL2, expand = expansion(add = c(.62, .72))) +
  scale_y_continuous(breaks = seq(0, 150, 50), expand = expansion(mult = c(0, .04))) +
  labs(title = "Cross-criterion support", subtitle = "all 143 genes per phenotype",
       x = NULL, y = "Genes") +
  base

## ══ C ═════════════════════════════════════════════════════
Cd <- melt(S[, .(Phenotype, `NJ vs MP` = Spearman_NJ_MP, `NJ vs ML` = Spearman_NJ_ML,
                 `MP vs ML` = Spearman_MP_ML)],
           id.vars = "Phenotype", variable.name = "Pair", value.name = "rho")
pC <- ggplot(Cd, aes(Pair, Phenotype, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.6) +
  geom_text(aes(label = sprintf("%.2f", rho), colour = rho > .55),
            size = 2.3, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradientn(colours = SEQ, limits = c(0.15, 0.80),
                       breaks = c(0.2, 0.4, 0.6, 0.8),
                       name = "Spearman ρ",
                       guide = guide_colourbar(barwidth = unit(4, "pt"),
                                               barheight = unit(26, "pt"),
                                               ticks = FALSE, frame.colour = NA)) +
  scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "white")) +
  scale_y_discrete(limits = rev(PHL), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(title = "Rank concordance of gene-level associations",
       subtitle = "Spearman \u03c1 across all 143 genes",
       x = NULL, y = NULL) +
  flat + theme(plot.margin = margin(13, 4, 2, 2))

## ══ D ═════════════════════════════════════════════════════
## 유의 판정 기준과 표시값을 일치시킨다:
##   연속형(NREM ratio, Sleep duration) -> BH 보정 P (Q 컬럼)
##   범주형(Sleep timing, Sleep frequency) -> raw permutation P (P 컬럼)
CONT <- c("NREM ratio","Sleep duration")
CAND <- data.table(
  Gene = c("ARNTL2","CPT1A","ADRB1","ATF4","CAVIN3","CRY2","PER1","ATF5"),
  Phenotype = c("NREM ratio","NREM ratio","Sleep duration","Sleep duration",
                "Sleep timing","Sleep timing","Sleep timing","Sleep frequency"))
D  <- merge(CAND, A, by = c("Gene","Phenotype"))
Dd <- melt(D[, .(Gene, Phenotype, NJ = Sig_NJ, MP = Sig_MP, ML = Sig_ML)],
           id.vars = c("Gene","Phenotype"), variable.name = "Method", value.name = "Sig")
Dq <- melt(D[, .(Gene, Phenotype, NJ = Q_NJ, MP = Q_MP, ML = Q_ML)],
           id.vars = c("Gene","Phenotype"), variable.name = "Method", value.name = "Q")
Dp <- melt(D[, .(Gene, Phenotype, NJ = P_NJ, MP = P_MP, ML = P_ML)],
           id.vars = c("Gene","Phenotype"), variable.name = "Method", value.name = "P")
Dd <- merge(merge(Dd, Dq, by = c("Gene","Phenotype","Method")),
            Dp, by = c("Gene","Phenotype","Method"))
Dd[, val := ifelse(as.character(Phenotype) %in% CONT, Q, P)]
stopifnot(all((Dd$val < 0.05) == Dd$Sig))               # 표시값 = 판정 기준

Dd[, lab := ifelse(val < 1e-3, formatC(val, format = "e", digits = 1),
                               formatC(val, format = "f", digits = 3))]
Dd[, lab := sub("e-0*([0-9]+)$", "\u00b710\u207b\\1", lab)]
sup <- c("1"="\u00b9","2"="\u00b2","3"="\u00b3","4"="\u2074","5"="\u2075",
         "6"="\u2076","7"="\u2077","8"="\u2078","9"="\u2079")
for (k in names(sup)) Dd[, lab := gsub(paste0("\u207b", k, "$"), paste0("\u207b", sup[[k]]), lab)]

Dd[, Phenotype := factor(as.character(Phenotype), levels = PHL)]
Dd[, Gene   := factor(Gene, levels = rev(CAND$Gene))]
Dd[, Method := factor(Method, levels = c("NJ","MP","ML"))]   # 전 패널 동일 순서

pD <- ggplot(Dd, aes(Method, Gene, fill = Sig)) +
  geom_tile(colour = "white", linewidth = 1.6) +
  geom_text(aes(label = lab, colour = Sig), size = 1.85, show.legend = FALSE) +
  facet_grid(Phenotype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(`TRUE` = POS, `FALSE` = NEG), name = NULL,
                    breaks = c(TRUE, FALSE), labels = c("significant", "not significant"),
                    guide = guide_legend(override.aes = list(colour = "grey80", linewidth = .3),
                                         keywidth = unit(8, "pt"), keyheight = unit(6, "pt"))) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK2)) +
  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Candidate-gene association significance across criteria",
       subtitle = paste0("Values are BH-adjusted P for continuous traits and raw ",
                         "permutation P for categorical traits"),
       x = NULL, y = NULL) +
  flat +
  theme(axis.text.y = element_text(face = "italic", size = 6.2),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 5.8, colour = INK2, hjust = 1),
        panel.spacing.y = unit(3, "pt"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key = element_rect(colour = NA),
        plot.margin = margin(13, 2, 0, 2))

## ══ 조립 ══════════════════════════════════════════════════
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.02), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = INK),
        plot.tag.position = c(0.004, 0.992))

MM <- 1/25.4
for (fm in c("jpg","pdf")) {
  f <- file.path(FIG, paste0("SupplementaryFigure10_reproducibility.", fm))
  if (fm == "pdf")
    ggsave(f, fig, width = 180*MM, height = 168*MM, device = cairo_pdf, bg = "white")
  else
    ggsave(f, fig, width = 180*MM, height = 168*MM, dpi = 600,
           device = "jpeg", quality = 98, bg = "white")
}
## 투고용 TIFF (LZW 압축, 600 dpi)
ft <- file.path(FIG, "SupplementaryFigure10_reproducibility.tif")
tiff(ft, width = 180*MM, height = 168*MM, units = "in", res = 600,
     compression = "lzw", type = "cairo", bg = "white")
print(fig); invisible(dev.off())

cat("저장 완료: jpg / pdf / tif  (180 x 168 mm, 600 dpi)\n[DONE]\n")
