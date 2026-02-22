# Data directory

This document describes the **final version** of the data layout: input files and output locations for the Sleep Evolution analyses. All paths are relative to the project root.

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

---

## Nucleotide matrices

Gene-wise nucleotide (and CDS) matrices are used by Result2, Result3, Result4, and Result5. In the **final version** there are two distinct locations:

| Use      | Directory | Scripts |
|----------|------------|--------|
| Result2  | `data/Result2/NucleotideMatrix/` | `Result2_visualization.R`, `Result2_AA_mutation_analysis.R` |
| Result3–5| `data/Circadian_gene_Nucleotide_Matrix/` | `Result3_visualization.R`, `Result4_visualization.R`, `Result5_SNP_Cochran_sleep_frequency.R`, `Result5_Sleepfrequency_associated_AA_effect_from_SNP.R` |

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
| Input  | `Sleeptiming_Cochran_Result.tsv`, `Sleeptiming_Manhattan_Dataset.tsv`, `Sleeptiming_AA_Mutation.tsv` (from upstream analyses) |
| Input  | `data/Fasta/`, `data/Circadian_gene_Nucleotide_Matrix/`, `data/Heatmap/Sleep_timing.tsv` |

---

### Result5 — `data/Result5/`

| Role   | File / path |
|--------|-------------|
| Input  | `Sleep_frequency_Cochran.tsv` |
| Input  | `data/Circadian_gene_Nucleotide_Matrix/<Gene>.tsv` |
| Output | Significant SNPs are written to `figures/Result5/SNP_Significant.tsv` by `Result5_SNP_Cochran_sleep_frequency.R` |
| Output | `figures/Result5/AminoAcid_Mutation_DF.tsv` from `Result5_Sleepfrequency_associated_AA_effect_from_SNP.R` |

---

### Result6 — `data/Result6/`

Phenotype-specific expression matrices, metadata, and DEG results. See **[data/Result6/README.md](Result6/README.md)** for file structure, matrix construction, and DEG method.

| Role   | File pattern |
|--------|--------------|
| Input  | `*_Count.tsv`, `*_Metadata.tsv` (and optionally `*_DEGs.tsv`) |
| Output | Figures from `Result6_visualization_DEG_boxplot.R` are saved under `figures/Result6/` |

---

## Heatmap data — `data/Heatmap/`

Summary matrices used by visualization scripts:

| File                     | Used by            |
|--------------------------|--------------------|
| `Total_sleep_time.tsv`   | Result2 scripts    |
| `NREM_ratio.tsv`         | Result3 scripts    |
| `Sleep_timing.tsv`       | Result4 scripts    |
| `Number_of_sleep_times.tsv` | Optional / referenced as needed |

Exact column names and formats are defined by the scripts that read these files (see [scripts/README.md](../scripts/README.md)).

