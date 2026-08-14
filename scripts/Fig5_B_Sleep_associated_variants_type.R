############################################################
# Fig5_B_Sleep_associated_variants_type.R
#
# Figure 5B. Fraction of sleep-associated variant types
#   Lollipop bar plot (monomorphic vs exonic subcategories)
#
# Input:
#   - data/Cross_validation/TableS22.txt  (sleep-associated SNP list)
#   - data/Cross_validation/TableS23.txt  (variant annotation)
#
# Output:
#   - Result/Fig5/Fig5_B_Sleep_associated_variants_type.jpg
#   - Result/Fig5/Fig5_B_Sleep_associated_variants_type.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)


## =========================================================
## 1. Paths
## =========================================================

PROJECT_DIR <- getwd()

TABLE_S22 <- file.path(PROJECT_DIR, "data", "Cross_validation",
                       "TableS22.txt")
TABLE_S23 <- file.path(PROJECT_DIR, "data", "Cross_validation",
                       "TableS23.txt")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


## =========================================================
## 2. Data
## =========================================================

supple_22 <- fread(TABLE_S22) |> as.data.frame()
supple_23 <- fread(TABLE_S23) |> as.data.frame()

useGene <- unique(supple_22$Gene)


## =========================================================
## 3. Variant categories
## =========================================================

supple_23$Category <- ifelse(supple_23$MAF == 0, "Monomorphic", "Non-monomorphic")
supple_23$Subcategory <- ifelse(
  supple_23$Category == "Monomorphic",
  "Monomorphic",
  supple_23$ExonicFunc.refGene
)
supple_23$Subcategory[is.na(supple_23$Subcategory)] <- "Not exonic"

df <- supple_23 %>%
  group_by(Category, Subcategory) %>%
  summarise(N = n(), .groups = "drop")

total <- sum(df$N)
df$frac <- df$N / total

df <- df %>%
  mutate(label = paste0(
    N, "/", sum(N), "\n(",
    round(frac * 100, 1), "%)"
  ))

df$Subcategory <- factor(df$Subcategory, levels = c(
  "Monomorphic",
  "nonsynonymous SNV",
  "synonymous SNV",
  "stopgain",
  "startloss",
  "Not exonic"
))

cols <- c(
  "Monomorphic"       = "#E64B35",
  "nonsynonymous SNV" = "#4DBBD5",
  "synonymous SNV"    = "#00A087",
  "stopgain"          = "#3C5488",
  "startloss"         = "#F39B7F",
  "Not exonic"        = "#999999"
)


## =========================================================
## 4. Plot
## =========================================================

p1 <- ggplot(df, aes(x = Subcategory, y = frac,
                     fill = Subcategory, colour = Subcategory)) +
  geom_segment(aes(x = Subcategory, xend = Subcategory,
                   y = 0, yend = frac),
               linewidth = 10) +
  geom_point(shape = 21, size = 9.2, stroke = 0.5) +
  geom_text(aes(label = label, y = frac + 0.05),
            vjust = 0, size = 4.5, colour = "black", lineheight = 0.92) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  scale_y_continuous(labels = percent,
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1.22),
                     expand = c(0, 0)) +
  theme_classic(base_size = 18) +
  labs(x = NULL, y = "Fraction of variants") +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    axis.text.y  = element_text(margin = margin(r = 4)),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin  = margin(t = 20, r = 14, b = 10, l = 18)
  )


## =========================================================
## 5. Save
## =========================================================

ggsave(file.path(OUT_DIR, "Fig5_B_Sleep_associated_variants_type.jpg"),
       p1, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Fig5_B_Sleep_associated_variants_type.pdf"),
       p1, width = 8, height = 7, bg = "white")

message("Sleep-associated genes (TableS22): ", length(useGene))
message("Variants plotted (TableS23): ", total)
message("Saved: Result/Fig5/Fig5_B_Sleep_associated_variants_type.jpg|.pdf")
