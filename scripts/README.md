## Scripts overview

This directory contains R scripts used for phylogenetic and statistical
analyses of sleep-related traits.


### 00_setup.R

Initial setup script for loading required R packages and defining
common settings used across analyses.


### Result2and3_anova_enrichment.R

Performs phylogenetic clustering of MUSCLE-aligned CDS sequences and
tests associations between cluster assignments and sleep phenotypes.

Analyses included:
- Total sleep time (Result2)
- NREM sleep ratio (Result3)

Outputs gene-level association statistics to:
- `data/Result2/`
- `data/Result3/`


### Result2_phylogenetic_tree_plot.R

Generates phylogenetic dendrograms for genes significantly associated
with total sleep time.

Trees are annotated by sleep phenotype categories
(Short Sleep, Others, Long Sleep) based on species-level metadata.

This script is intended for exploratory visualization and figure
generation and does not automatically save plots by default.
