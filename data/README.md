# Data directory

Input files and intermediate tables for Sleep Evolution analyses. All paths are relative to the project root.

**Convention:** `data/Result1`–`Result6` are organized by **phenotype**. Main manuscript panels are written by `Fig*` scripts to `Result/Fig3|Fig4|Fig5/` (see root [README.md](../README.md)).

---

## Metadata (Result1)

Species-level sleep phenotype metadata is required for all analyses.

| Item     | Location |
|----------|----------|
| File     | `species_sleep_metadata.txt` |
| Directory| `data/Result1/` |

The file contains species identifiers and sleep-related phenotypes (total sleep time, NREM ratio, sleep frequency, sleep timing) used by the scripts. Optional: `circadian_Gene_list.tsv` for gene lists where referenced.

---

## FASTA (aligned CDS)

MUSCLE-aligned CDS FASTA files are required per gene for phylogenetic and visualization scripts.

| Item     | Location / format |
|----------|--------------------|
| Directory| `data/Fasta/` |
| Naming   | `[GENE]_muscle.fasta` |

Each file contains CDS alignments across species (e.g. from an ortholog selection and MUSCLE alignment pipeline).

**Alignment provenance**

The aligned CDS FASTA files provided in `data/Fasta/` correspond to the final filtered alignment set used in the manuscript. Multiple sequence alignments were generated using MUSCLE v5 and subsequently filtered using GUIDANCE2. 

---

## Nucleotide matrices

Gene-wise nucleotide (and CDS) matrices are used by Result2, Result3, Result4, and Result5. In the **final version** there are two distinct locations:

| Use      | Directory | Scripts |
|----------|------------|--------|
| Result2  | `data/Result2/NucleotideMatrix/` | `Result2_AA_mutation_analysis.R`, `Fig3_*` |
| Result3–5| `data/Circadian_gene_Nucleotide_Matrix/` | Result3–5 analysis scripts, `Fig3_*` / `Fig4_*` |

- One matrix file per gene, e.g. `[GENE].tsv`.
- Ensure the correct directory is populated for the scripts you run (Result2 uses `Result2/NucleotideMatrix`, Result3–5 use `Circadian_gene_Nucleotide_Matrix`).

---

## Result directories (inputs and outputs)

### Result1 — `data/Result1/`

| Role   | File / note |
|--------|-------------|
| Input  | `species_sleep_metadata.txt` (required) |
| Optional| `circadian_Gene_list.tsv` |

---

### Result2 — `data/Result2/`

| Role   | File / path |
|--------|-------------|
| Input  | `species_sleep_metadata.txt` (via Result1) |
| Input  | FASTA in `data/Fasta/` |
| Input  | `data/Result2/NucleotideMatrix/<Gene>.tsv` for visualization and AA mutation |
| Input  | `clustering.permutation200.Sleep_real.tsv` (for `Result2_AA_mutation_analysis.R`) |
| Input  | `ALL_SNP.tsv`, `Total_Sleep_Time_SNP.tsv`, `TST_key_12_optimal_Kruskal.tsv`, `Phylogenetic_Signal_Data_permutation.tsv` as produced by upstream steps |
| Output | `Total_sleep_time_Anova_Result.tsv` (from `Result2and3_anova_enrichment.R`) |
| Output | `Total_sleep_time_AA_Mutation_Data.tsv` (from `Result2_AA_mutation_analysis.R`) |

---

### Result3 — `data/Result3/`

| Role   | File / path |
|--------|-------------|
| Input  | `data/Fasta/`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Heatmap/NREM_ratio.tsv` |
| Output | `NREM_ratio_Anova_Result.tsv`, `NREM_key_12_optimal_Kruskal.tsv` (from `Result2and3_anova_enrichment.R`) |
| Output | `NREM_ratio_phylogeneticsignal.tsv`, `NREM_ratio_SNP.tsv`, `NREM_ratio_AminoAcid_mutation.tsv` (from other analyses) |

---

### Result4 — `data/Result4/`

| Role   | File |
|--------|------|
| Input / output | `Sleep_Timing_optimal_K.tsv` (optimal cluster K) |
| Output | `Sleeptiming_Fisher_Result.tsv` — gene-level Fisher permutation (`Gene`, `Cluster`, `P_Value`) from `Result4_Sleeptiming_tree_association_test.R` |
| Input  | `Sleeptiming_Manhattan_Dataset.tsv`, `Sleeptiming_AA_Mutation.tsv` (from upstream analyses) |
| Input  | `data/Fasta/`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Heatmap/Sleep_timing.tsv` |

Fig4 timing panels (A/B/C/E) read `Sleeptiming_Fisher_Result.tsv` (not the older Cochran enrichment table).

---

### Result5 — `data/Result5/`

| Role   | File / path |
|--------|-------------|
| Input / output | `Number_of_sleep_optimal_K.tsv` (optimal cluster K) |
| Output | `Sleep_frequency_Fisher.tsv` — gene-level Fisher permutation (`Gene`, `Cluster`, `P_Value`) from `Result5_Sleepfrequency_tree_association_test.R` |
| Input  | `SNP.tsv`, `All_Indel_Split.tsv` |
| Input  | `data/Circadian_gene_Nucleotide_Matrix/<Gene>.tsv` (or nucleotide matrices as configured) |
| Output | Significant SNPs / AA tables from `Result5_SNP_Cochran_sleep_frequency.R` and `Result5_Sleepfrequency_associated_AA_effect_from_SNP.R` (paths as set in those scripts) |

Fig4 frequency panels (A/F) and Result5 SNP gene filtering read `Sleep_frequency_Fisher.tsv`.

---

### Cross-validation — `data/Cross_validation/`

Cross-validation / human constraint inputs used by **Fig5** scripts.

| File | Role |
|------|------|
| `TableS22.txt` | Sleep-associated SNP / gene list (Fig5B, Fig5D) |
| `TableS23.txt` | Variant annotation (Func / ExonicFunc, CADD, etc.; Fig5B, Fig5C) |
| `TableS24.txt` | Gene-level LOEUF / pLI for target genes (Fig5D) |
| `TableS25.txt` | Human GWAS overlap for translational relevance (Fig5E) |
| `TableS27.txt` | RNA-seq run / BioProject metadata |
| `TableS28.txt` | DEG summary across four sleep phenotypes |
| `gnomAD_LOEUF_result_all.tsv` | Background gnomAD LOEUF table (Fig5D) |
| `supplementaryS2.tsv` | Circadian gene annotation reference |
| `gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` | Raw gnomAD LoF metrics (source archive) |

---

### Result6 — `data/Result6/`

Phenotype-specific expression matrices, metadata, and DEG results. See **[data/Result6/README.md](Result6/README.md)** for file structure, matrix construction, and DEG method.

| Role   | File pattern |
|--------|--------------|
| Input  | `*_Count.tsv`, `*_Metadata.tsv` (and optionally `*_DEGs.tsv`) |
| Output | DEG boxplots via `Fig5_H_DEG_boxplot.R` → `Result/Fig5/Fig5_H_DEG_boxplot/` |

**RNA-seq data provenance**

Raw RNA-seq data were processed as described in the Methods section (adapter trimming, alignment, sorting, and gene-level quantification). The count matrices provided in this directory correspond to the finalized processed datasets used for differential expression and regression analyses in the manuscript. These files enable full reproduction of transcriptomic-level statistical analyses without requiring reprocessing of raw FASTQ files.

---

## Heatmap data — `data/Heatmap/`

Summary matrices used by visualization scripts:

| File                     | Used by            |
|--------------------------|--------------------|
| `Total_sleep_time.tsv`   | Fig3 / Result2     |
| `NREM_ratio.tsv`         | Fig3 / Result3     |
| `Sleep_timing.tsv`       | Fig4 / Result4     |
| `Number_of_sleep_times.tsv` | Fig4 / Result5 as needed |

**Evolutionary statistics provenance**

Summary statistics reflecting evolutionary pressure, including nucleotide diversity (π), Tajima’s D, and dN/dS ratios, were calculated using MEGA-CC v10.2.6 as described in the Methods section. The values provided in this directory correspond to the finalized outputs used for figure generation in the manuscript.

---
Exact column names and formats are defined by the scripts that read these files (see [scripts/README.md](../scripts/README.md)).

## Reproducibility note

This repository is designed for manuscript-level reproducibility. All aligned sequences, nucleotide matrices, association result tables, and processed RNA-seq count matrices correspond to the finalized datasets used in the published analyses. Software versions and statistical procedures are described in the Methods section. Permutation-based analyses were performed with fixed random seeds to ensure reproducibility.
