## Scripts overview

This directory contains R scripts used for phylogenetic and statistical
analyses of sleep-related traits.

---

### 00_setup.R

Initial setup script for loading required R packages and defining
common settings used across analyses.

---

### Result2and3_anova_enrichment.R

Performs phylogenetic clustering of MUSCLE-aligned CDS sequences and
tests associations between cluster assignments and sleep phenotypes.

Analyses included:
- Total sleep time (Result2)
- NREM sleep ratio (Result3)

Outputs gene-level association statistics to:
- `data/Result2/`
- `data/Result3/`

---

### Result2_phylogenetic_tree_plot.R

Generates phylogenetic dendrograms for genes significantly associated
with total sleep time.

Trees are annotated by sleep phenotype categories
(Short Sleep, Others, Long Sleep) based on species-level metadata.

This script is intended for exploratory visualization and figure
generation and does not automatically save plots by default.

---

## Result5 — Sleep frequency–associated analyses

### Result5_SNP_Cochran_sleep_frequency.R

Performs site-level association testing between nucleotide variants
and sleep frequency (number of sleep times per day) using the
Cochran–Armitage trend test.

Outputs significant SNPs to:
- `figures/Result5/SNP_Significant.tsv`

Note:
In Result5, very few significant SNPs were detected at the nucleotide level.
Therefore, additional analyses were conducted to investigate whether
structural variants (insertions and deletions; indels) contribute to
sleep frequency differences.

---

### Result5_FigBCDE_sleeptiming_associated_gene.R

Generates figure panels for sleep frequency–associated genes
(ATF5 and NFIL3):

- Fig5B/E: Phylogenetic trees colored by sleep frequency phenotype
- Fig5C/D: Mosaic plots comparing cluster majority vs observed phenotype

Input:
- `data/Result1/species_sleep_metadata.txt`
- `data/Fasta/<Gene>_muscle.fasta`
- `data/Result5/Sleep_frequency_Cochran.tsv`

Output:
- `figures/Result5/Figure5B_Tree_<Gene>.pdf`
- `figures/Result5/Figure5C_Mosaic_<Gene>.pdf`

---

### Indel identification workflow (external tools)

To complement SNP analyses in Result5, indel variants were identified
using external command-line tools:

1. Multiple sequence alignments (FASTA) were converted to VCF format
   using `jvarkit msa2vcf --ignore-n-bases`.

2. Indel variants were extracted using:
   `bcftools view -v indels`.

3. Variant statistics and site-level annotations were generated using:
   `bcftools stats` and `bcftools query`.

4. Indels were classified as insertions or deletions based on
   relative sequence lengths of REF (human reference) and ALT alleles.

This workflow is documented for reproducibility but relies on
standard bioinformatics tools rather than custom R implementations.

---

## Notes

- All analyses assume the project root structure:
  `Sleep_Evolution/`
- Input data are located under `data/`
- Output figures and tables are written to `figures/`
- Scripts are modular and can be executed independently
  for each result section.
