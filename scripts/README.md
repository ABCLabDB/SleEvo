# Scripts overview

R scripts for sleep-evolution analyses and **main figure** generation.

All paths are relative to the **project root** (folder containing `data/`, `scripts/`, `Result/`).

- **`Result*` numbers** = phenotype-organized **analysis** (and some supporting plots).
- **`Fig3` / `Fig4` / `Fig5`** = **manuscript figures**:
  - **Fig3** — Sleep duration & NREM ratio  
  - **Fig4** — Sleep timing & sleep frequency  
  - **Fig5** — Cross-validation  

Run `00_setup.R` first. Prefer **`Fig*`** scripts to regenerate main panels (`Result/Fig3|Fig4|Fig5/`).

---

## Recommended run order

1. `00_setup.R`
2. Phenotype analyses as needed (`Result2and3_anova_enrichment.R`, SNP/AA scripts, Cochran scripts, …)
3. Main figures: `Fig3_*.R` → `Fig4_*.R` → `Fig5_*.R`

---

## Setup

### 00_setup.R

Install/load core packages (`data.table`, `dplyr`, `ggplot2`, `ggpubr`, `DescTools`, …). Run from project root.

---

## Fig3 — Sleep duration & NREM ratio

Outputs: `Result/Fig3/`

| Script | Description |
|--------|-------------|
| `Fig3_AB_Sleepduration_tree.R` | Phylogenetic trees + sleep-duration bars (e.g. ADRB1, ATF4) |
| `Fig3_C_Phylogenetic_signal.R` | Blomberg’s K / Moran’s I (alphabetical first 20 genes) |
| `Fig3_D_Selection_signature_heatmap.R` | dN/dS vs Tajima’s D selection signature |
| `Fig3_E_NREMratio.R` | NREM-ratio tree (e.g. ARNTL2) |
| `Fig3_F_NREMratio_Manhattanplot.R` | NREM SNP Manhattan plot |
| `Fig3_G_ARNTL2_NREMratio_variants.R` | ARNTL2 focal SNP alignment + logo/boxplot |

Typical inputs: `data/Result1/`, `data/Result2/`, `data/Result3/`, `data/Fasta/`, `data/Heatmap/`.

---

## Fig4 — Sleep timing & sleep frequency

Outputs: `Result/Fig4/`

| Script | Description |
|--------|-------------|
| `Fig4_A_Sleeptiming_enrichment.R` | Cochran enrichment lollipops (timing + frequency; combined panel) |
| `Fig4_B_PER1_sleeptiming_tree.R` | PER1 tree + categorical sleep-timing heatmap |
| `Fig4_C_Selection_signature_heatmap.R` | Selection signature for timing candidates |
| `Fig4_D_Manhattanplot.R` | Sleep-timing Manhattan plot |
| `Fig4_E_PER1_AminoAcid_mutation.R` | PER1 nonsynonymous protein lollipop (UniProt domains) |
| `Fig4_F_Sleepfrequency_tree.R` | ATF5 / NFIL3 frequency trees (+ side-by-side combined) |
| `Fig4_G_Sleepfrequency_mutation.R` | ATF5 & CARTPT SNP alignment windows |
| `Fig4_H_NFIL3_sleepfrequency_indel.R` | NFIL3 indel schematic |

Typical inputs: `data/Result4/`, `data/Result5/`, `data/Fasta/`, `data/Circadian_gene_Nucleotide_Matrix/`.

---

## Fig5 — Cross-validation

Outputs: `Result/Fig5/` (DEG plots under `Result/Fig5/Fig5_H_DEG_boxplot/`)

| Script | Description |
|--------|-------------|
| `Fig5_B_Sleep_associated_variants_type.R` | Variant-type fractions (TableS23) |
| `Fig5_C_CADD.R` | Non-synonymous CADD ≥ 20 horizontal bar |
| `Fig5_D_LOEUF.R` | LOEUF &lt; 0.6 bar + permutation null (TableS22/S24 + gnomAD) |
| `Fig5_E_OPN4_translational_relevance.R` | OPN4 Chr10:86658604 ATLAS vs Human GWAS (TableS25) |
| `Fig5_H_DEG_boxplot.R` | DEG boxplots (tst / nrem / timing / frequency) from `data/Result6/` |

Typical inputs: `data/Cross_validation/`, `data/Result6/`, `data/Fasta/`.

---

## Phenotype analysis (`Result*`)

### Clustering / association

- **`Result2and3_anova_enrichment.R`** — CDS distance → NJ → clustering; ANOVA/Kruskal for total sleep & NREM.  
  Output: `data/Result2/Total_sleep_time_Anova_Result.tsv`, `data/Result3/NREM_ratio_Anova_Result.tsv`.

### Total sleep time (Result2)

- **`Result2_sleep_duration_variants.R`** — total-sleep SNPs → `data/Result2/Total_sleep_time_SNP.tsv` (or similar).
- **`Result2_AA_mutation_analysis.R`** — AA mutations for those SNPs.

### NREM ratio (Result3)

- **`Result3_nrem_ratio_variants.R`** — NREM SNPs → `data/Result3/NREM_Ratio_SNP.tsv`.
- **`Result3_AA_mutation_analysis.R`** — AA mutations for NREM SNPs.
- **`Result3_ARNTL2_variants_LD.R`** — LD-like r² heatmap among NREM nonsynonymous SNPs.

### Sleep timing (Result4)

- **`Result4_Sleeptiming_optimal_cluster_group.R`** — optimal K / group assignment.
- **`Result4_Cochran_enrichment.R`** — gene-level Cochran enrichment → `data/Result4/Sleeptiming_Cochran_Result.tsv`.

**SNP testing note:** biallelic sites use Cochran–Armitage; multiallelic sites use χ².

### Sleep frequency (Result5)

- **`Result5_SNP_Cochran_sleep_frequency.R`** — site-level Cochran–Armitage → significant SNP table.
- **`Result5_Sleepfrequency_associated_AA_effect_from_SNP.R`** — AA effects from those SNPs.
- **`Result5_Figures_visualization.R`** — legacy Result5 figure bundle (prefer `Fig4_*` / `Fig5_*` for main panels).

### Metadata supplementary figures

- **`Supplementary_Fig1to6_sleep_architecture_metadata.R`** — Supp. Figs 1–6 from `data/Result1/species_sleep_metadata.txt`.

---

## Notes

- Run from project root; Fig* scripts create `Result/FigN/` as needed.
- See [data/README.md](../data/README.md) for input layouts.
- Indel calling for frequency analyses used external tools (bcftools / jvarkit); see project documentation / Result5 notes.
