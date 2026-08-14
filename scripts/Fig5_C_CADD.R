############################################################
# Fig5_C_CADD.R
#
# Figure 5C. Non-synonymous variants split by CADD score
#   Horizontal two-segment bar: CADD >= 20 (teal) vs < 20 (grey)
#
# Input:
#   - data/Cross_validation/TableS23.txt  (variant annotation)
#
# Output:
#   - Result/Fig5/Fig5_C_CADD.jpg
#   - Result/Fig5/Fig5_C_CADD.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

TABLE_S23 <- file.path(PROJECT_DIR, "data", "Cross_validation",
                       "TableS23.txt")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Data and categories
## =========================================================

supple_23 <- fread(TABLE_S23) |> as.data.frame()

supple_23$Category <- ifelse(supple_23$MAF == 0, "Monomorphic", "Non-monomorphic")
supple_23$Subcategory <- ifelse(
  supple_23$Category == "Monomorphic",
  "Monomorphic",
  supple_23$ExonicFunc.refGene
)
supple_23$Subcategory[is.na(supple_23$Subcategory)] <- "Not exonic"


## =========================================================
## 3. Non-synonymous only + CADD split
## =========================================================

df_ns <- supple_23 %>%
  filter(Category == "Non-monomorphic",
         Subcategory == "nonsynonymous SNV") %>%
  mutate(CADD_group = ifelse(CADD_phred >= 20, "CADD ≥ 20", "CADD < 20"))

total_ns <- nrow(df_ns)
n_hi     <- sum(df_ns$CADD_group == "CADD ≥ 20", na.rm = TRUE)


## =========================================================
## 4. Panel C-style horizontal bar (geom_rect, fixed layout)
## =========================================================

col_hi <- "#3AB0A2"   # teal (CADD >= 20)
col_lo <- "#D9D9D9"   # rest (grey)

bar_y <- 1.1
bar_h <- 0.8

rect_df <- data.frame(
  xmin = c(0, n_hi),
  xmax = c(n_hi, total_ns),
  ymin = bar_y,
  ymax = bar_y + bar_h,
  grp  = c("hi", "lo")
)

p2 <- ggplot() +
  ## two-segment bar
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = grp),
            colour = NA) +
  ## label inside teal segment
  annotate("text", x = n_hi / 2, y = bar_y + bar_h / 2,
           label = paste0(n_hi, " variants (CADD ≥ 20)"),
           colour = "white", fontface = "bold", size = 4.2) +
  ## top bracket + total count
  annotate("segment", x = 0, xend = total_ns, y = 2.25, yend = 2.25,
           colour = "black", linewidth = 0.5) +
  annotate("text", x = total_ns / 2, y = 2.42,
           label = paste0(total_ns, " non-synonymous variants"),
           fontface = "bold", size = 4.8) +
  ## "Functionally significant" callout (down, then arrow right)
  annotate("segment", x = n_hi * 0.15, xend = n_hi * 0.15,
           y = bar_y, yend = 0.55, colour = "black", linewidth = 0.5) +
  annotate("segment", x = n_hi * 0.15, xend = n_hi * 0.6,
           y = 0.55, yend = 0.55, colour = "black", linewidth = 0.5,
           arrow = arrow(length = unit(0.10, "inches"), type = "closed")) +
  annotate("text", x = n_hi * 0.65, y = 0.55, hjust = 0,
           label = "Functionally significant",
           fontface = "bold", size = 4.2) +
  scale_fill_manual(values = c(hi = col_hi, lo = col_lo), guide = "none") +
  coord_cartesian(xlim = c(0, total_ns), ylim = c(0, 3), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(6, 30, 6, 6))


## =========================================================
## 5. Save (wide and short, panel-C aspect)
## =========================================================

ggsave(file.path(OUT_DIR, "Fig5_C_CADD.pdf"),
       plot = p2, width = 6.5, height = 2.4)
ggsave(file.path(OUT_DIR, "Fig5_C_CADD.jpg"),
       plot = p2, width = 6.5, height = 2.4, dpi = 300, bg = "white")

message("Non-synonymous variants: ", total_ns,
        " | CADD >= 20: ", n_hi)
message("Saved: Result/Fig5/Fig5_C_CADD.jpg|.pdf")
