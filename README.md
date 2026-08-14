# Sleep Evolution (SleEvo)

Analysis and figure code for evolutionary associations between circadian genes and sleep phenotypes across species.

**`Result*` numbers** organize analyses and data by **sleep phenotype** (metadata, total sleep time, NREM ratio, sleep timing, sleep frequency, RNA-seq DEG).

**`Fig*` scripts** generate the **main manuscript figures**, grouped by biological theme (not by Result number):

| Figure | Theme | Primary phenotypes |
|--------|--------|--------------------|
| **Fig3** | Sleep duration & NREM ratio | Total sleep time, NREM ratio |
| **Fig4** | Sleep timing & sleep frequency | Sleep timing, sleep frequency |
| **Fig5** | Cross-validation | Human constraint / GWAS / transcriptomics |

---

## Repository structure

All paths are relative to the **project root** (the folder that contains `data/`, `scripts/`, and `Result/`).

```
SleEvo/
├── data/                          # Inputs & intermediate analysis tables
│   ├── Circadian_gene_Nucleotide_Matrix/
│   ├── Cross_validation/          # Fig5 tables (TableS22–S28, LOEUF, etc.)
│   ├── Fasta/                     # MUSCLE-aligned CDS FASTA
│   ├── Heatmap/
│   ├── Result1/ … Result6/        # Phenotype-organized analysis inputs
│   └── README.md
│
├── Result/                        # Main figure outputs (from Fig* scripts)
│   ├── Fig3/
│   ├── Fig4/
│   └── Fig5/
│
├── scripts/
│   ├── 00_setup.R
│   ├── Fig3_*.R … Fig5_*.R        # Main figure generation (preferred)
│   ├── Result*_*.R                # Phenotype analysis pipelines
│   ├── Supplementary_Fig1to6_*.R
│   └── README.md
│
└── README.md
```

- **Data formats:** [data/README.md](data/README.md)
- **Script details:** [scripts/README.md](scripts/README.md)

---

## How Result* relates to Fig*

- **Result1–6** = analysis modules / data folders by phenotype (and supporting tables).
- **Fig3–5** = publication figure panels. They **reuse** Result* outputs (and related tables) but are organized by the figure story:

| Fig | Uses (typical) |
|-----|----------------|
| **Fig3** | Result2 (total sleep time) + Result3 (NREM ratio) |
| **Fig4** | Result4 (sleep timing) + Result5 (sleep frequency) |
| **Fig5** | Cross-validation tables (`data/Cross_validation/`) + Result6 DEG expression |

Older monolithic Result* *visualization* scripts live under `scripts/99.previous/` for reference. Prefer the modular **`Fig*`** scripts for regenerating main figures.

---

## Main figures (`Fig*` scripts)

Run from the **project root**. Outputs go to `Result/Fig3`, `Result/Fig4`, or `Result/Fig5` (jpg + pdf unless noted).

### Fig3 — Sleep duration & NREM ratio

| Script | Panel / content |
|--------|------------------|
| `Fig3_AB_Sleepduration_tree.R` | **3A/B** Phylogenetic trees for total-sleep–associated genes (e.g. ADRB1, ATF4) with sleep-duration bars |
| `Fig3_C_Phylogenetic_signal.R` | **3C** Blomberg’s K / Moran’s I across genes (alphabetical top 20) |
| `Fig3_D_Selection_signature_heatmap.R` | **3D** dN/dS vs Tajima’s D selection-signature scatter |
| `Fig3_E_NREMratio.R` | **3E** NREM-ratio phylogenetic tree (e.g. ARNTL2) |
| `Fig3_F_NREMratio_Manhattanplot.R` | **3F** Manhattan plot of NREM-associated SNPs |
| `Fig3_G_ARNTL2_NREMratio_variants.R` | **3G** ARNTL2 focal SNP alignment + logo/boxplot |

### Fig4 — Sleep timing & sleep frequency

| Script | Panel / content |
|--------|------------------|
| `Fig4_A_Sleeptiming_enrichment.R` | **4A** Cochran enrichment lollipops (timing + frequency; also combined) |
| `Fig4_B_PER1_sleeptiming_tree.R` | **4B** PER1 tree + sleep-timing heatmap column |
| `Fig4_C_Selection_signature_heatmap.R` | **4C** Selection signature for sleep-timing candidates |
| `Fig4_D_Manhattanplot.R` | **4D** Sleep-timing Manhattan plot |
| `Fig4_E_PER1_AminoAcid_mutation.R` | **4E** PER1 nonsynonymous lollipop (UniProt domains) |
| `Fig4_F_Sleepfrequency_tree.R` | **4F** ATF5 / NFIL3 sleep-frequency trees (+ combined) |
| `Fig4_G_Sleepfrequency_mutation.R` | **4G** ATF5 & CARTPT SNP alignment windows |
| `Fig4_H_NFIL3_sleepfrequency_indel.R` | **4H** NFIL3 indel schematic |

### Fig5 — Cross-validation

| Script | Panel / content |
|--------|------------------|
| `Fig5_B_Sleep_associated_variants_type.R` | **5B** Fraction of variant types (monomorphic vs exonic classes) |
| `Fig5_C_CADD.R` | **5C** Non-synonymous variants split by CADD ≥ 20 |
| `Fig5_D_LOEUF.R` | **5D** Gene-level LOEUF intolerance bar + permutation null |
| `Fig5_E_OPN4_translational_relevance.R` | **5E** OPN4 Chr10:86658604 ATLAS vs Human GWAS |
| `Fig5_H_DEG_boxplot.R` | **5H** DEG expression boxplots → `Result/Fig5/Fig5_H_DEG_boxplot/` |

---

## Phenotype analysis scripts (`Result*`)

These produce **tables / statistics** (or supporting panels) organized by phenotype. They feed Fig* scripts or supplementary analyses.

| Script | Role |
|--------|------|
| `Supplementary_Fig1to6_sleep_architecture_metadata.R` | Supp. Figs 1–6 from species sleep metadata (Result1) |
| `Result2and3_anova_enrichment.R` | Phylogenetic clustering + ANOVA/Kruskal for total sleep & NREM |
| `Result2_sleep_duration_variants.R` | Total-sleep–associated SNPs |
| `Result2_AA_mutation_analysis.R` | Amino-acid mutations for total-sleep SNPs |
| `Result3_nrem_ratio_variants.R` | NREM-associated SNPs |
| `Result3_AA_mutation_analysis.R` | Amino-acid mutations for NREM SNPs |
| `Result3_ARNTL2_variants_LD.R` | LD-like co-evolution heatmap among NREM nonsynonymous SNPs |
| `Result4_Cochran_enrichment.R` | Sleep-timing gene-level Cochran enrichment |
| `Result4_Sleeptiming_optimal_cluster_group.R` | Optimal cluster grouping for sleep timing |
| `Result5_SNP_Cochran_sleep_frequency.R` | Sleep-frequency SNP Cochran–Armitage tests |
| `Result5_Sleepfrequency_associated_AA_effect_from_SNP.R` | AA effects of frequency-associated SNPs |
| `Result5_Figures_visualization.R` | Legacy Result5 figure bundle (prefer modular `Fig4_*` / `Fig5_*` where applicable) |

---

## Requirements

- **R** (≥ 4.x recommended)
- Core packages via `scripts/00_setup.R` (`data.table`, `dplyr`, `ggplot2`, …)
- Additional packages used by Fig* scripts as needed, e.g. `ape`, `Biostrings`, `ggtree`, `ggrepel`, `patchwork`, `DescTools`

---

## Quick start

```r
setwd("path/to/SleEvo")   # project root
source("scripts/00_setup.R")

# Example: regenerate Fig3C
source("scripts/Fig3_C_Phylogenetic_signal.R")
```

1. Place inputs as described in [data/README.md](data/README.md).
2. Run phenotype analyses (`Result*`) when you need to refresh tables.
3. Run **`Fig*`** scripts to write publication panels under `Result/Fig3|Fig4|Fig5/`.

---

## Reproducibility

- Run all scripts from the **repository root**.
- Fig* scripts write to `Result/FigN/` with relative paths only.
- Phenotype data remain under `data/Result1`–`Result6` and `data/Cross_validation/`.
