# Supplementary three-method comparison

Supplementary Tables 12–13, the four supporting supplementary tables, and the reproducibility Supplementary Figure, together with every input and script needed to regenerate them.

Criteria are ordered **NJ → MP → ML** in every table, figure panel, and column name. NJ is neighbour joining on T92 distances (the published analysis), MP is maximum parsimony, and ML is maximum likelihood.

In this repository the inputs and outputs live under `data/Supplementary_three_method/`, and the scripts live under `scripts/Supplementary_three_method/`.

---

## `03_supplementary_table` — six supplementary tables

| File | Rows × columns | Content |
|---|---|---|
| `SupplementaryTable12_three_method_comparison.tsv` | 572 × 14 | **Table 12.** One row per gene × phenotype. `P_*`, `Q_*`, and `Sig_*` for each criterion, and how many criteria are significant |
| `SupplementaryTable13_tree_quality_and_reproducibility.tsv` | 4 × 30 | **Table 13.** One row per phenotype. MP tree quality (PS, CI, RI, RC, bootstrap, equally best trees, consensus) + ML tree quality (parsimony-informative sites, logL, SH-aLRT, UFBoot) + three-criterion reproducibility |
| `SupplementaryTable_MP_tree_statistics.tsv` | 572 × 27 | Gene-level source table for the MP part of Table 13 |
| `SupplementaryTable_ML_tree_statistics.tsv` | 572 × 13 | Gene-level source table for the ML part of Table 13 (selected model, logL, AIC, BIC, site counts) |
| `SupplementaryTable_reproducibility_by_phenotype.tsv` | 4 × 20 | Reproducibility, three-criterion overlap, minimum possible overlap, and Spearman ρ. Input to Table 13 and figure panels A and C |
| `SupplementaryTable_three_method_summary.tsv` | 4 × 18 | Summary by criterion combination (NJ∩MP, NJ∩ML, MP∩ML, union, and so on) |

## figure

`SupplementaryFigure10_reproducibility.{jpg,pdf,tif}` is 180 × 168 mm at 600 dpi. The jpg is 600 dpi, the pdf is vector, and the tif is LZW-compressed for submission. `SupplementaryFigure10_legend.txt` is the figure legend.

- **A** — number of significant genes per phenotype
- **B** — reproducibility level (3/3, 2/3, 1/3, 0/3)
- **C** — pairwise Spearman ρ
- **D** — cross-criterion consistency for representative candidate genes

---

## Layout

```
data/Supplementary_three_method/
├── 00_input/                 starting inputs
│   ├── alignments/                   143 × <GENE>_muscle.fasta (MUSCLE; not realigned)
│   ├── METADATA_fixedVersion.txt     68 species × sleep phenotypes and taxonomy
│   ├── supplementaryS2.tsv           gene list
│   ├── Sleep_Timing_optimal_K.tsv    fixed K for sleep timing (carried over from the original analysis)
│   └── Number_of_sleep_optimal_K.tsv fixed K for sleep frequency
├── 01_trees/                 rebuilt gene trees
│   ├── MP/              572   <phenotype>__<gene>_MP.nwk, node labels = bootstrap %
│   ├── MP_equalbest/    513   bundles of equally most-parsimonious trees (median 12 per file; 15,235 trees in total)
│   └── ML/<phenotype>/  572   <gene>_ML.nwk, node labels = SH-aLRT/UFBoot
├── 02_association/           association results by criterion (12 files)
│   ├── NJ/   4   original analysis
│   ├── MP/   4   optimalK ×2 (continuous), Fisher_perm20k ×2 (categorical)
│   └── ML/   4   same layout as MP
├── 03_supplementary_table/   six supplementary tables (table above)
└── figure/                   jpg, pdf, tif, and legend

scripts/Supplementary_three_method/   scripts 01–07, in run order
```

**572 = 143 genes × 4 phenotypes.** The species set differs by phenotype (sleep duration, 68 species; sleep timing, 55; NREM ratio, 45; sleep frequency, 42), so trees were built separately for each phenotype. A single tree cannot be subset later and still yield K groups of the species actually analysed.

**Why `MP_equalbest` has 513 files rather than 572:** for 59 gene–phenotype combinations the most-parsimonious tree was unique, so no equal-best file was written. PER1 (sleep timing) is one of them. This folder is used for the categorical sensitivity analysis (Fisher *P* recomputed on up to 20 equally best trees; `EqualBest_*` columns).

---

## Run order

| # | Script | Reads | Writes |
|---|---|---|---|
| 01 | `01_mp_trees.R` | `00_input/` | `01_trees/MP/`, `01_trees/MP_equalbest/`, MP tree statistics |
| 02 | `02_mp_association.R` | `01_trees/MP*`, `00_input/` | `02_association/MP/` |
| 03 | `03_ml_trees_iqtree.py` | `00_input/` | `01_trees/ML/`, ML tree statistics |
| 04 | `04_ml_association.R` | `01_trees/ML/`, `00_input/` | `02_association/ML/` |
| 05 | `05_three_method_comparison.R` | `02_association/{NJ,MP,ML}/` | **Table 12**, `three_method_summary`, `reproducibility_by_phenotype` |
| 06 | `06_tree_quality_table.py` | MP and ML tree statistics, `01_trees/ML/`, output of 05 | **Table 13** |
| 07 | `07_supplementary_figure.R` | Table 12, `reproducibility_by_phenotype` | **Supplementary Figure** |

Scripts 01–04 rebuild trees and take about 10 hours (on 128 cores: MP about 8 h, ML about 20 min, and each association step about 1 h). Because `02_association/` is already included, **05 → 06 → 07 alone regenerates the supplementary tables and the figure** (a few seconds).

```bash
OMP_NUM_THREADS=1 Rscript scripts/Supplementary_three_method/05_three_method_comparison.R
python3            scripts/Supplementary_three_method/06_tree_quality_table.py
OMP_NUM_THREADS=1 Rscript scripts/Supplementary_three_method/07_supplementary_figure.R
```

Checked on 2026-09-28: re-running 05 → 06 → 07 reproduced the four supplementary tables above byte for byte.

R library path used for the original run: `/disk2/bijsy/Test/Evolution_pressure/0.Code/Rlib` (`data.table`, `ape`, `phangorn`, `Biostrings`, `ggplot2`, `patchwork`, `ggrepel`).

---

## Paths inside the scripts

The scripts were copied **as they were actually run**, as a reproducibility record. Their internal paths still point at the original locations. Read them with the mapping below. Paths in this table are relative to `data/Supplementary_three_method/`.

| Path inside the script | Location in this bundle |
|---|---|
| `/disk2/bijsy/Test/Muscle/` | `00_input/alignments/` |
| `…/Sleep_association_Test/METADATA_fixedVersion.txt` | `00_input/METADATA_fixedVersion.txt` |
| `…/Sleep_association_Test/supplementaryS2.tsv` | `00_input/supplementaryS2.tsv` |
| `…/Result/Sleeptiming/Sleep_Timing_optimal_K.tsv` | `00_input/Sleep_Timing_optimal_K.tsv` |
| `…/Result/Sleepfrequency/Number_of_sleep_optimal_K.tsv` | `00_input/Number_of_sleep_optimal_K.tsv` |
| `…/v2/trees/`, `…/v2/trees_equalbest/` | `01_trees/MP/`, `01_trees/MP_equalbest/` |
| `…/ML_based/trees/<phenotype>/` | `01_trees/ML/<phenotype>/` |
| `…/Result/NREM_optimalK_ANOVA_selectedP_1M_FAST.tsv` and the other three NJ files | `02_association/NJ/` |
| `…/v2/tables/MP_v2_*`, `…/ML_based/tables/ML_*` | `02_association/{MP,ML}/` |
| `…/v2/tables/SupplementaryTable_MP_tree_statistics_v2.tsv` | `03_supplementary_table/SupplementaryTable_MP_tree_statistics.tsv` |
| `…/ML_based/tables/SupplementaryTable_ML_tree_statistics.tsv` | `03_supplementary_table/` (same file name) |
| `…/3method_comparison/tables/` | `03_supplementary_table/` |
| `…/3method_comparison/figures/` | `figure/` |

---

## Statistical definitions

All three criteria use the same downstream procedure: tree → `cophenetic.phylo` → average-linkage hierarchical clustering → `cutree(K)`. Only the tree changes.

- **Continuous traits** (sleep duration, NREM ratio): one-way ANOVA on cluster membership, *K* = 2–12, and the *K* with the largest *F* is kept. Significance is **Benjamini–Hochberg adjusted *P* < 0.05**, so `Q_*` determines `Sig_*`.
- **Categorical traits** (sleep timing, sleep frequency): Fisher’s exact test on the 2 × *K* table at the fixed *K* carried over from the original analysis, with a 20,000-permutation correction. Significance is the **uncorrected permutation *P* < 0.05**, so `P_*` determines `Sig_*`. The adjusted value is also stored in `Q_*`.

This distinction is stated in the subtitle of Supplementary Figure panel D and in the figure legend.

---

## Deliberately excluded

Anything not used by the three deliverables above was left out.

| Excluded | Original location | Reason |
|---|---|---|
| Strict consensus trees | `v2/consensus/` | Resolution is already summarised in the MP tree-statistics table |
| Raw IQ-TREE output | `ML_based/iqtree_raw/` | The needed values are in the ML tree statistics and in the Newick node labels |
| Checkpoints | `v2/checkpoint/`, `ML_based/checkpoint/` | Only for resuming a re-run |
| Node-level bootstrap (13 MB) | `v2/tables/node_level_bootstrap_MP_v2.tsv` | Only the summary enters Table 13 |
| Fig1–Fig4 and `Figure_three_method_summary` | `3method_comparison/figures/` | Replaced by the final Supplementary Figure |
| v1 MP results | `99.Pre/`, `tables/` | Replaced by v2 |
| `Code/04_topology_comparison.R` and its output | `Maximum_Parsimony_for_Cladistics/` | Bug: `as.splits()` includes terminal splits, so bootstrap matching is broken. Do not cite |
