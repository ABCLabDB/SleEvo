# Scripts overview 

R scripts for phylogenetic and statistical analyses of sleep-related traits. All paths are relative to the **project root** (the folder containing `data/`, `scripts/`, and `figures/`).

**Run `00_setup.R` first** to install and load core packages. Then run scripts in an order consistent with data dependencies (see recommended order below).

---

## Recommended run order

1. **00_setup.R** — install/load packages.
2. **Result2and3_anova_enrichment.R** — produces Result2 and Result3 association tables (requires `data/Fasta/`, `data/Result1/species_sleep_metadata.txt`).
3. Other Result2/Result3 analyses and visualizations as needed (some require precomputed permutation/SNP files in `data/Result2/` and `data/Result3/`).
4. **Result4_visualization.R**, **Result5_***, **Result6_visualization_DEG_boxplot.R** — require corresponding input data under `data/Result4/`, `data/Result5/`, and `data/Result6/` as described in [data/README.md](../data/README.md).

Scripts are modular and can be run independently provided the required inputs exist.

---

## Setup

### 00_setup.R

- **Purpose:** Install (if missing) and load packages used across analyses.
- **Packages:** `data.table`, `dplyr`, `ggplot2`, `ggpubr`, `DescTools`.
- **Run from:** Project root.

---

## Result1 — Sleep architecture metadata

### Result1_SuppFig1to6_sleep_architecture_metadata.R

- **Purpose:** Generate Supplementary Figures 1–6 (total sleep time, NREM ratio, sleep frequency, sleep timing, primate timing–frequency, global timing–frequency mosaic).
- **Input:** `data/Result1/species_sleep_metadata.txt`.
- **Output:** Figures displayed in R session (save manually or adapt script).
- **See:** [data/README.md](../data/README.md).

---

## Result2 & Result3 — Clustering and association

### Result2and3_anova_enrichment.R

- **Purpose:** Phylogenetic clustering of MUSCLE-aligned CDS (genetic distance → NJ tree → cophenetic distance → hierarchical clustering) and association with sleep phenotypes (ANOVA / Kruskal–Wallis).
- **Input:** `data/Fasta/*_muscle.fasta`, `data/Result1/species_sleep_metadata.txt`.
- **Output:** `data/Result2/Total_sleep_time_Anova_Result.tsv`, `data/Result3/NREM_ratio_Anova_Result.tsv`.

### Result2_visualization.R

- **Purpose:** Result2 figures: phylogenetic trees (2A), box/strip plots (2B), phylogenetic signal (2E), SNP-centered nucleotide view (2F), selection heatmap (2H).
- **Input:** `data/Result1/species_sleep_metadata.txt`, `data/Result2/TST_key_12_optimal_Kruskal.tsv`, `data/Result2/Phylogenetic_Signal_Data_permutation.tsv`, `data/Result2/Total_Sleep_Time_SNP.tsv`, `data/Result2/NucleotideMatrix/<Gene>.tsv`, `data/Fasta/`, `data/Heatmap/Total_sleep_time.tsv`.
- **Output:** `figures/Result2/PhyloTree/`, `figures/Result2/Boxplot/`, `figures/Result2/Phylogenetic_Signal/`, `figures/Result2/SNP/`, `figures/Result2/Heatmap/`.

### Result2_AA_mutation_analysis.R

- **Purpose:** Amino-acid mutation analysis for sleep-associated SNPs.
- **Input:** `data/Result1/species_sleep_metadata.txt`, `data/Result2/clustering.permutation200.Sleep_real.tsv`, `data/Result2/ALL_SNP.tsv`, `data/Result2/NucleotideMatrix/<Gene>.tsv`.
- **Output:** `data/Result2/Total_sleep_time_AA_Mutation_Data.tsv`.

### Result2_AA_PML_visualization.R

- **Purpose:** PML amino-acid mutation visualization (e.g. Figure 2I/2J).
- **Input:** `data/Result2/Total_sleep_time_AA_Mutation_Data.tsv`.
- **Output:** `figures/Result2/AminoAcid/PML_Mutation.pdf`.

### Result3_visualization.R

- **Purpose:** Result3 figures: NREM phylogenetic tree (3A), phylogenetic signal (3C), heatmap (3D), Manhattan (3E), ARNTL2 mutation map.
- **Input:** `data/Result1/species_sleep_metadata.txt`, `data/Result3/NREM_key_12_optimal_Kruskal.tsv`, `data/Result3/NREM_ratio_Anova_Result.tsv`, `data/Result3/NREM_ratio_phylogeneticsignal.tsv`, `data/Result3/NREM_ratio_SNP.tsv`, `data/Result3/NREM_ratio_AminoAcid_mutation.tsv`, `data/Fasta/`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Heatmap/NREM_ratio.tsv`.
- **Output:** `figures/Result3/`.

---

## Result4 — Sleep timing

### Result4_visualization.R

- **Purpose:** Result4 figures: lollipop enrichment (4A), sleep timing dendrogram (4B), evolution quadrant (4C), Manhattan (4D), PER1 amino acid landscape (4E), SuppFig12 SNP heatmap.
- **Input:** `data/Result4/Sleeptiming_Cochran_Result.tsv`, `Sleeptiming_Manhattan_Dataset.tsv`, `Sleeptiming_AA_Mutation.tsv`, `data/Result1/species_sleep_metadata.txt`, `data/Fasta/`, `data/Heatmap/Sleep_timing.tsv`, `data/Circadian_gene_Nucleotide_Matrix/`.
- **Output:** `figures/Result4/`.

---

## Result5 — Sleep frequency

### Result5_SNP_Cochran_sleep_frequency.R

- **Purpose:** Site-level Cochran–Armitage trend test between nucleotide variants and sleep frequency.
- **Input:** `data/Result1/species_sleep_metadata.txt`, `data/Result5/Sleep_frequency_Cochran.tsv`, `data/Circadian_gene_Nucleotide_Matrix/<Gene>.tsv`.
- **Output:** `figures/Result5/SNP_Significant.tsv`.

### Result5_Sleepfrequency_associated_AA_effect_from_SNP.R

- **Purpose:** Amino-acid-level effects of sleep frequency–associated SNPs (translate, classify synonymous/non-synonymous, annotate property change).
- **Input:** `data/Result1/species_sleep_metadata.txt`, `figures/Result5/SNP_Significant.tsv`, `data/Circadian_gene_Nucleotide_Matrix/<Gene>.tsv`.
- **Output:** `figures/Result5/AminoAcid_Mutation_DF.tsv`.

### Result5_Figures_visualization.R

- **Purpose:** Result5 figure panels: Cochran lollipop (5A), phylogenetic trees for ATF5/NFIL3 (5B/E), mosaic plots (5C/D).
- **Input:** `data/Result1/species_sleep_metadata.txt`, `data/Result5/Sleep_frequency_Cochran.tsv`, `data/Fasta/<Gene>_muscle.fasta`.
- **Output:** `figures/Result5/*.pdf`.

### Indel identification workflow (external tools)

To complement SNP analyses in Result5, indel variants were identified
using external command-line tools.

Tool versions:
- bcftools v1.10.2
- Jvarkit (Git commit dea54f1a2)

Workflow:

1. Multiple sequence alignments (FASTA) were converted to VCF format
   using:
   `jvarkit msa2vcf --ignore-n-bases`

2. Indel variants were extracted using:
   `bcftools view -v indels`

3. Variant statistics and site-level annotations were generated using:
   `bcftools stats` and `bcftools query`

4. Indels were classified as insertions or deletions based on
   the relative sequence lengths of the REF (human reference)
   and ALT alleles.

This workflow documents the exact software versions and parameters
used to ensure reproducibility of the indel identification step.

---

## Result6 — DEG visualization

### Result6_visualization_DEG_boxplot.R

- **Purpose:** Boxplot visualization for DEGs (e.g. Figure 6B). Linear model: `Expression ~ Group + Species:Sequencer`.
- **Input:** `data/Result6/*_Metadata.tsv`, `data/Result6/*_Count.tsv` (e.g. `Total_sleep_time_Metadata.tsv`, `Total_sleep_time_Count.tsv`, and similarly for NREM_ratio, Sleep_timing, Sleep_frequency).
- **Output:** `figures/Result6/<gene>_<suffix>.pdf`.

---

## Permutation-based cluster association (sleep timing & sleep frequency)

For categorical sleep traits, an empirical permutation test was used:

1. Species clustered by hierarchical clustering (average linkage) on T92 genetic distance.
2. Cochran–Armitage test statistic computed between cluster labels and phenotype categories.
3. Cluster labels permuted (e.g. N = 1,000,000); empirical p-value: `P = (|T_perm| >= |T_obs| + 1) / (N + 1)`.

Implemented in R (e.g. DescTools). Indel workflow used external tools (bcftools, Jvarkit) as documented in the project.

---

## Notes

- All scripts assume the **project root** contains `data/` and (for most) `figures/`; directories are created as needed.
- Scripts are modular and can be run independently per result section, provided required inputs exist as in [data/README.md](../data/README.md).
