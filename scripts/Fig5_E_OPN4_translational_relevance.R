############################################################
# Fig5_E_OPN4_translational_relevance.R
#
# Figure 5E. Translational relevance — OPN4 Chr10:86658604
#   Panel-E style: curated ATLAS species (night C / daytime T),
#   5-bp window, colored focal blocks, Human GWAS allele row.
#
# Input:
#   - data/Fasta/OPN4_muscle.fasta
#   - data/Result1/species_sleep_metadata.txt
#   - data/Result4/SleepTiming_SNP_for_Manhattan.tsv
#   - data/Cross_validation/TableS25.txt
#
# Output:
#   - Result/Fig5/Fig5_E_OPN4_Chr10_86658604.jpg
#   - Result/Fig5/Fig5_E_OPN4_Chr10_86658604.pdf
#
# Required working directory: repository root (SleEvo)
############################################################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


## =========================================================
## 1. Paths / constants
## =========================================================

PROJECT_DIR <- getwd()

FASTA_FILE <- file.path(PROJECT_DIR, "data", "Fasta", "OPN4_muscle.fasta")
META_FILE  <- file.path(PROJECT_DIR, "data", "Result1",
                        "species_sleep_metadata.txt")
SNP_MAP    <- file.path(PROJECT_DIR, "data", "Result4",
                        "SleepTiming_SNP_for_Manhattan.tsv")
TABLE_S25  <- file.path(PROJECT_DIR, "data", "Cross_validation",
                        "TableS25.txt")

OUT_DIR <- file.path(PROJECT_DIR, "Result", "Fig5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

gene         <- "OPN4"
target_label <- "Chr10:86658604"
aln_idx      <- 667L                 # Real_Position in muscle alignment
## Panel E shows five coordinates around the focal SNP
win_idx      <- (aln_idx - 2L):(aln_idx + 2L)   # 665..669 -> 86658602..86658606


## =========================================================
## 2. Curated Panel-E species (display order: top -> bottom)
## =========================================================

## Night group (allele C) then daytime group (allele T)
panel_species <- data.frame(
  display = c(
    "Human", "Chimpanzee", "Gibbon", "Vervet-agm",
    "Bolivian squirrel monkey", "Squirrel", "Zebra finch", "Fruit fly",
    "Siberian musk deer", "Pig", "Horse", "Ferret",
    "Microbat", "Armadillo", "Opossum", "Japanese quail"
  ),
  symbol = c(
    "homo_sapiens", "pan_troglodytes", "nomascus_leucogenys",
    "chlorocebus_sabaeus", "saimiri_boliviensis_boliviensis",
    "ictidomys_tridecemlineatus", "taeniopygia_guttata",
    "drosophila_melanogaster",
    "moschus_moschiferus", "sus_scrofa", "equus_caballus",
    "mustela_putorius_furo", "myotis_lucifugus",
    "dasypus_novemcinctus", "monodelphis_domestica",
    "coturnix_japonica"
  ),
  group = c(rep("Sleep at night", 8), rep("Sleep at daytime", 8)),
  stringsAsFactors = FALSE
)
panel_species$y <- rev(seq_len(nrow(panel_species)))   # Human at top


## =========================================================
## 3. Colors (Panel E)
## =========================================================

col_night <- "#B7D7A8"     # green block (C / sleep at night)
col_day   <- "#D5C6E0"     # purple block (T / sleep at daytime)
col_human_row <- "#D6EAF8" # Human ATLAS row
col_gwas_bg   <- "#FCE4D6" # Human GWAS strip
col_gwas_dot  <- "#7B5EA7" # GWAS allele circle
col_grey_txt  <- "grey65"
col_focus_txt <- "grey10"


## =========================================================
## 4. Load alignment bases from FASTA (base R, no Biostrings)
## =========================================================

read_muscle_fasta <- function(path) {
  lines <- readLines(path)
  out <- list(); cur <- NULL; seq <- character()
  for (ln in lines) {
    if (startsWith(ln, ">")) {
      if (!is.null(cur)) out[[cur]] <- paste(seq, collapse = "")
      cur <- sub("^>", "", ln)
      seq <- character()
    } else {
      seq <- c(seq, gsub("\\s+", "", ln))
    }
  }
  if (!is.null(cur)) out[[cur]] <- paste(seq, collapse = "")
  out
}

fa <- read_muscle_fasta(FASTA_FILE)

## map ensembl symbol -> sequence
sym2seq <- list()
for (nm in names(fa)) {
  parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
  sym <- if (length(parts) >= 2) tolower(parts[2]) else NA_character_
  if (!is.na(sym)) sym2seq[[sym]] <- toupper(fa[[nm]])
}

missing <- setdiff(panel_species$symbol, names(sym2seq))
if (length(missing) > 0) {
  stop("Missing FASTA sequences for: ", paste(missing, collapse = ", "))
}


## =========================================================
## 5. Genomic labels for window columns
## =========================================================

snp_map <- fread(SNP_MAP) |> as.data.frame()
snp_map <- snp_map[snp_map$Gene == gene, ]

lab_for <- function(rp) {
  hit <- snp_map$Column_Index[snp_map$Real_Position == rp][1]
  if (!is.na(hit) && grepl("^Chr", hit)) return(hit)
  if (!is.na(hit) && grepl("^CDS_", hit)) {
    return(paste0("Chr10:", sub("^CDS_", "", hit)))
  }
  ## interpolate from neighbors if needed
  prev <- snp_map$Column_Index[snp_map$Real_Position == rp - 1][1]
  if (!is.na(prev) && grepl("866586", prev)) {
    g <- as.integer(sub(".*:", "", prev)) + 1L
    return(paste0("Chr10:", g))
  }
  paste0("V", rp)
}

pos_labs <- vapply(win_idx, lab_for, character(1))
names(pos_labs) <- as.character(win_idx)
## ensure focal label
pos_labs[as.character(aln_idx)] <- target_label


## =========================================================
## 6. Build long data for ATLAS panel
## =========================================================

atlas_rows <- list()
for (i in seq_len(nrow(panel_species))) {
  sym <- panel_species$symbol[i]
  seq <- sym2seq[[sym]]
  for (j in seq_along(win_idx)) {
    rp <- win_idx[j]
    atlas_rows[[length(atlas_rows) + 1]] <- data.frame(
      Species  = panel_species$display[i],
      group    = panel_species$group[i],
      y        = panel_species$y[i],
      Position = pos_labs[as.character(rp)],
      x        = j,
      Base     = substr(seq, rp, rp),
      is_focus = (rp == aln_idx),
      stringsAsFactors = FALSE
    )
  }
}
atlas_df <- bind_rows(atlas_rows)
atlas_df$Species <- factor(atlas_df$Species,
                           levels = rev(panel_species$display))

## focal fill by phenotype group
atlas_df$Fill <- "white"
atlas_df$Fill[atlas_df$is_focus & atlas_df$group == "Sleep at night"] <-
  col_night
atlas_df$Fill[atlas_df$is_focus & atlas_df$group == "Sleep at daytime"] <-
  col_day

atlas_df$TxtCol <- ifelse(atlas_df$is_focus, col_focus_txt, col_grey_txt)
atlas_df$TxtFace <- ifelse(atlas_df$is_focus, "bold", "plain")


## =========================================================
## 7. Human GWAS row (TableS25)
## =========================================================

s25 <- fread(TABLE_S25) |> as.data.frame()
s25_hit <- s25[s25$Gene == gene &
                 as.integer(s25$Position) == 86658604L, ][1, ]
gwas_alt <- as.character(s25_hit$ALT)
message("GWAS ALT=", gwas_alt, " | p=", signif(s25_hit$`p-value`, 4),
        " | ", s25_hit$Phenotype_human_GWAS)

gwas_df <- data.frame(
  Position = pos_labs,
  x = seq_along(pos_labs),
  Base = ifelse(pos_labs == target_label, gwas_alt, ""),
  stringsAsFactors = FALSE
)
gwas_y <- 0.15


## =========================================================
## 8. Geometry for blocks / arrows
## =========================================================

x_focus <- which(pos_labs == target_label)
n_pos   <- length(pos_labs)

night_ys <- panel_species$y[panel_species$group == "Sleep at night"]
day_ys   <- panel_species$y[panel_species$group == "Sleep at daytime"]
human_y  <- panel_species$y[panel_species$display == "Human"]

block_df <- data.frame(
  xmin = x_focus - 0.48,
  xmax = x_focus + 0.48,
  ymin = c(min(night_ys) - 0.45, min(day_ys) - 0.45),
  ymax = c(max(night_ys) + 0.45, max(day_ys) + 0.45),
  fill = c(col_night, col_day),
  row.names = NULL
)

arrow_df <- data.frame(
  x = x_focus + 0.55,
  xend = n_pos + 0.55,
  y = c(mean(night_ys), mean(day_ys), gwas_y),
  lab = c("Sleep at night", "Sleep at daytime", "Sleep overall"),
  col = c("#2c6e49", "#7B5EA7", "#7B5EA7"),
  row.names = NULL
)


## =========================================================
## 9. Plot
## =========================================================

p <- ggplot() +
  ## Human ATLAS row background
  annotate("rect",
           xmin = 0.5, xmax = n_pos + 0.5,
           ymin = human_y - 0.48, ymax = human_y + 0.48,
           fill = col_human_row, colour = NA) +
  ## focal phenotype blocks
  geom_rect(data = block_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill),
            colour = NA, inherit.aes = FALSE) +
  ## nucleotides
  geom_text(data = atlas_df,
            aes(x = x, y = y, label = Base,
                colour = TxtCol, fontface = TxtFace),
            size = 4.2, inherit.aes = FALSE) +
  ## Human GWAS strip background
  annotate("rect",
           xmin = 0.5, xmax = n_pos + 0.5,
           ymin = gwas_y - 0.42, ymax = gwas_y + 0.42,
           fill = col_gwas_bg, colour = NA) +
  ## GWAS allele circle + letter
  annotate("point",
           x = x_focus, y = gwas_y,
           size = 9, shape = 21,
           fill = col_gwas_dot, colour = col_gwas_dot) +
  annotate("text",
           x = x_focus, y = gwas_y,
           label = gwas_alt, colour = "white",
           fontface = "bold", size = 4.2) +
  ## left labels: species + Human GWAS
  scale_y_continuous(
    breaks = c(panel_species$y, gwas_y),
    labels = c(panel_species$display, "Human GWAS"),
    expand = expansion(add = 0.35)
  ) +
  scale_x_continuous(
    breaks = seq_len(n_pos),
    labels = pos_labs,
    position = "top",
    limits = c(0.4, n_pos + 2.6),
    expand = c(0, 0)
  ) +
  ## bold focal x label via overlay annotation
  annotate("text",
           x = x_focus, y = max(panel_species$y) + 0.95,
           label = target_label, fontface = "bold", size = 3.3,
           colour = "grey10") +
  ## arrows + phenotype labels
  geom_segment(data = arrow_df,
               aes(x = x, xend = xend, y = y, yend = y, colour = col),
               linewidth = 1.1, lineend = "round",
               arrow = arrow(length = unit(0.14, "inches"),
                             type = "closed"),
               inherit.aes = FALSE) +
  geom_text(data = arrow_df,
            aes(x = xend + 0.12, y = y, label = lab, colour = col),
            hjust = 0, fontface = "bold", size = 4.0,
            inherit.aes = FALSE) +
  ## Evolution based ATLAS side caption
  annotate("text",
           x = n_pos + 2.35,
           y = mean(range(panel_species$y)),
           label = "Evolution\nbased\nATLAS",
           colour = "grey55", size = 3.3, fontface = "italic",
           lineheight = 0.95) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(clip = "off",
                  ylim = c(gwas_y - 0.7, max(panel_species$y) + 1.15)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(
      angle = 40, hjust = 0, vjust = 0, size = 9, colour = "grey55"
    ),
    axis.text.y = element_text(size = 10, hjust = 0, colour = "grey20"),
    plot.margin = margin(t = 28, r = 70, b = 14, l = 10),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


## =========================================================
## 10. Save
## =========================================================

stem <- "Fig5_E_OPN4_Chr10_86658604"
ggsave(file.path(OUT_DIR, paste0(stem, ".jpg")),
       p, width = 7.2, height = 7.0, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")),
       p, width = 7.2, height = 7.0, bg = "white")

message("Panel E species: ", nrow(panel_species),
        " | window: ", paste(pos_labs, collapse = ", "))
message("Saved: Result/Fig5/", stem, ".jpg|.pdf")
