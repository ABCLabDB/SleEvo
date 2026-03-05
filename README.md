# Sleep Evolution

This repository contains analysis code and data structures used to study the evolutionary association between circadian genes and sleep phenotypes across species.

The analyses link phylogenetic clustering of gene sequences to species-level sleep traits: **total sleep time**, **NREM sleep ratio**, **sleep timing**, and **sleep frequency**. Downstream scripts produce phylogenetic trees, boxplots, heatmaps, and mutation maps.

---

## Repository structure

All paths are relative to the **project root** (the folder that contains `data/`, `scripts/`, and `figures/`).

```
Sleep_Evolution/
├── data/                          # Input data and intermediate results
│   ├── Circadian_gene_Nucleotide_Matrix/   # Gene-wise nucleotide matrices (Result3–5)
│   ├── Fasta/                              # MUSCLE-aligned CDS FASTA (per gene)
│   ├── Heatmap/                             # Phenotype summary matrices
│   ├── Result1/                             # Species-level sleep metadata
│   ├── Result2/                             # Total sleep time (incl. NucleotideMatrix/)
│   ├── Result3/                             # NREM ratio association
│   ├── Result4/                             # Sleep timing association
│   ├── Result5/                             # Sleep frequency association
│   ├── Result6/                             # RNA-seq expression & DEG (see data/Result6/README.md)
│   └── README.md
│
├── figures/                       # Generated plots (created by scripts)
│   ├── Result2/
│   ├── Result3/
│   ├── Result4/
│   ├── Result5/
│   └── Result6/
│
├── scripts/
│   ├── 00_setup.R
│   ├── Result1_SuppFig1to6_sleep_architecture_metadata.R
│   ├── Result2_visualization.R
│   ├── Result2_AA_mutation_analysis.R
│   ├── Result2_AA_PML_visualization.R
│   ├── Result2and3_anova_enrichment.R
│   ├── Result3_visualization.R
│   ├── Result4_visualization.R
│   ├── Result5_SNP_Cochran_sleep_frequency.R
│   ├── Result5_Sleepfrequency_associated_AA_effect_from_SNP.R
│   ├── Result5_Figures_visualization.R
│   ├── Result6_visualization_DEG_boxplot.R
│   └── README.md
│
├── .gitignore
└── README.md
```

- **Data and file formats:** [data/README.md](data/README.md)
- **Scripts and run order:** [scripts/README.md](scripts/README.md)

---

## Overview of analyses

| Result   | Phenotype   | Description |
|----------|-------------|-------------|
| **Result1** | Metadata     | Species-level sleep architecture (total sleep time, NREM ratio, sleep frequency, timing). Supplementary Figures 1–6. |
| **Result2** | Total sleep time | Phylogenetic clustering and association; SNP/AA mutation analysis; figures 2A–2J. |
| **Result3** | NREM ratio   | Clustering and association; phylogenetic signal; Manhattan/SNP profiles; ARNTL2 mutation map. |
| **Result4** | Sleep timing | Cochran–Armitage enrichment; dendrograms; Manhattan; PER1 amino acid landscape; SuppFig12. |
| **Result5** | Sleep frequency | SNP-level Cochran–Armitage test; AA effect from SNPs; figure panels (lollipop, trees, mosaic). |
| **Result6** | RNA-seq DEG  | Phenotype-specific expression matrices and linear-model DEG boxplots. |

---

## Requirements

- **R** (scripts are written in R).
- R packages (installed by `00_setup.R`): `data.table`, `dplyr`, `ggplot2`, `ggpubr`, `DescTools`. Additional packages are loaded by individual scripts (e.g. `ape`, `Biostrings`, `pheatmap`, `Cairo`); install as needed if not already present.

---

## Quick start

1. **Clone or copy** the repository so that the folder containing `data/`, `scripts/`, and `figures/` is your working directory (project root).
2. **Place input data** as described in [data/README.md](data/README.md) (metadata, FASTA, nucleotide matrices, and any precomputed Result files you need).
3. **Run from project root in R:**
   ```r
   setwd("path/to/Sleep_Evolution")   # project root
   source("scripts/00_setup.R")
   ```
4. **Run analyses** in an order that matches data dependencies (see [scripts/README.md](scripts/README.md)). For example:
   - `Result2and3_anova_enrichment.R` produces Result2/Result3 association tables.
   - Then run visualization scripts (Result2_visualization.R, Result3_visualization.R, etc.) as needed.

Figures are written under `figures/`; script-specific outputs are documented in [scripts/README.md](scripts/README.md).

---

## Reproducibility

- All scripts are intended to be run from the **repository (project) root**.
- Run `scripts/00_setup.R` first to install and load core packages.
- External data and file formats are documented in [data/README.md](data/README.md).
- Figure and table outputs are written to `figures/` and under `data/` as described in each script’s header and in [scripts/README.md](scripts/README.md).

---

## Scripts

Key R scripts in `scripts/` and their roles in the **final** pipeline:

- **00_setup.R**: Install and load core R packages (`data.table`, `dplyr`, `ggplot2`, `ggpubr`, `DescTools`, etc.). Always run once from the project root before other scripts.
- **Result1_SuppFig1to6_sleep_architecture_metadata.R**: Uses `data/Result1/species_sleep_metadata.txt` to generate Supplementary Figures 1–6 (sleep architecture overview).
- **Result2and3_anova_enrichment.R**: Builds phylogenetic clusters from MUSCLE-aligned CDS FASTA (`data/Fasta/*_muscle.fasta`) and tests association with total sleep time and NREM ratio (outputs `data/Result2/Total_sleep_time_Anova_Result.tsv` and `data/Result3/NREM_ratio_Anova_Result.tsv`).
- **Result2_visualization.R**: Visualizes Result2 (trees, box/strip plots, phylogenetic signal, SNP views, heatmaps) using `data/Result2/NucleotideMatrix/`, SNP tables, and `data/Heatmap/Total_sleep_time.tsv` (writes to `figures/Result2/`).
- **Result2_AA_mutation_analysis.R** / **Result2_AA_PML_visualization.R**: Analyze and visualize amino-acid mutations for total sleep time–associated SNPs (uses permutation/SNP tables in `data/Result2/` and writes mutation tables/plots).
- **Result3_visualization.R**: Visualizes NREM ratio results (trees, phylogenetic signal, Manhattan, ARNTL2 mutation map) using `data/Circadian_gene_Nucleotide_Matrix/` and `data/Heatmap/NREM_ratio.tsv` (writes to `figures/Result3/`).
- **Result4_visualization.R**: Sleep timing enrichment and PER1 amino-acid landscape (uses `data/Result4/*Sleeptiming*.tsv`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Heatmap/Sleep_timing.tsv`; writes to `figures/Result4/`).
- **Result5_SNP_Cochran_sleep_frequency.R**, **Result5_Sleepfrequency_associated_AA_effect_from_SNP.R**, **Result5_Figures_visualization.R**: SNP-level Cochran–Armitage tests for sleep frequency, amino-acid effect annotation, and main figure panels (inputs under `data/Result5/`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Fasta/`; outputs to `figures/Result5/`).
- **Result6_visualization_DEG_boxplot.R**: Reads phenotype-specific RNA-seq matrices (`data/Result6/*_Count.tsv`, `*_Metadata.tsv`) and produces DEG boxplots under `figures/Result6/` using a linear model `Expression ~ Group + Species:Sequencer`.

For full details (e.g. exact input file names and outputs), see [scripts/README.md](scripts/README.md).

---
