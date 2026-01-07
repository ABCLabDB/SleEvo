## Metadata

Species-level sleep phenotype metadata is required to reproduce the analyses.

Place the following file under the directory below:

- species_sleep_metadata.txt

Directory:
- data/Result1/

This metadata file contains species identifiers and sleep-related phenotypes
(e.g. total sleep time, NREM sleep ratio) used in downstream analyses.


## External FASTA data

This analysis requires MUSCLE-aligned CDS FASTA files for each gene.

Expected directory:
- data/Fasta/

File naming format:
- [GENE]_muscle.fasta

Each FASTA file contains CDS alignments across species.
These files were generated using a custom ortholog selection and
MUSCLE alignment pipeline and are provided in this repository.


## Result2: Total sleep time association

This directory contains gene-level association results between
phylogenetic clustering and total sleep time.

Directory:
- data/Result2/

Key output file:
- TST_key_12_optimal_Kruskal.tsv

This file is used as direct input for downstream visualization
(e.g. phylogenetic tree annotation and summary plots).


## Result3: NREM sleep ratio association

This directory contains gene-level association results between
phylogenetic clustering and NREM sleep ratio.

Directory:
- data/Result3/

Key output file:
- NREM_key_12_optimal_Kruskal.tsv

This file is used as direct input for downstream visualization
(e.g. phylogenetic tree annotation and summary plots).
