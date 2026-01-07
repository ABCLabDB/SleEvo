# Sleep_Evolution

This repository contains analysis code and data structures used to study
the evolutionary association between circadian genes and sleep phenotypes
across species.

The analyses focus on linking phylogenetic clustering of gene sequences
to species-level sleep traits, including total sleep time and NREM sleep ratio.


## Repository structure
Sleep_Evolution/
├── data/
│ ├── Fasta/ # MUSCLE-aligned CDS FASTA files (gene-wise)
│ ├── Result1/ # Species-level sleep metadata
│ ├── Result2/ # Total sleep time association results
│ ├── Result3/ # NREM sleep ratio association results
│ └── README.md
│
├── scripts/
│ ├── 00_setup.R
│ ├── Result2_phylogenetic_tree_plot.R
│ ├── Result2and3_anova_enrichment.R
│ └── README.md
│
├── figures/
│ └── Result2/ # Phylogenetic tree figures (optional)
│
└── README.md


## Overview of analyses

- **Result2**  
  Phylogenetic clustering of circadian genes and association with
  total sleep time across species.

- **Result3**  
  Phylogenetic clustering of circadian genes and association with
  NREM sleep ratio across species.

Downstream visualization scripts generate phylogenetic tree plots
annotated with sleep phenotype categories.


## Reproducibility

All scripts are written in R and are intended to be executed from
the repository root directory.

External data requirements and file formats are documented in
`data/README.md`.
