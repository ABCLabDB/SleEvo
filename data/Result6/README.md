# Result6: Phenotype-specific RNA-seq matrices and DEG analysis

This directory holds the **final version** of phenotype-specific expression matrices, sample metadata, and linear-model-based DEG results used by the Sleep Evolution project.

---

## File structure

For each sleep phenotype, the following files are expected under `data/Result6/`:

| Suffix            | Description |
|-------------------|-------------|
| `*_Count.tsv`     | Log2(FPKM+1) expression matrix (genes × samples). |
| `*_Metadata.tsv`  | Sample-level metadata for modeling (e.g. sample ID, Group, Species, Instrument). |
| `*_DEGs.tsv`      | Differentially expressed genes from the linear model (optional for visualization). |

### Phenotypes (file name prefixes)

- `Total_sleep_time_Count.tsv`, `Total_sleep_time_Metadata.tsv`, `Total_sleep_time_DEGs.tsv`
- `NREM_ratio_Count.tsv`, `NREM_ratio_Metadata.tsv`, `NREM_ratio_DEGs.tsv`
- `Sleep_timing_Count.tsv`, `Sleep_timing_Metadata.tsv`, `Sleep_timing_DEGs.tsv`
- `Sleep_frequency_Count.tsv`, `Sleep_frequency_Metadata.tsv`, `Sleep_frequency_DEGs.tsv`

The script `Result6_visualization_DEG_boxplot.R` reads the `*_Count.tsv` and `*_Metadata.tsv` files and writes boxplot PDFs to **`figures/Result6/`** (not under `data/Result6/`).

---

## Expression matrix construction

1. Gene-level counts were generated with **featureCounts**.
2. DESeq2 size-factor normalization was performed:

```r
dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData   = sample_metadata,
    design    = ~ 1
)

dds <- DESeq(dds)

mcols(dds)$basepairs <- gene_length_vector
fpkm_matrix <- fpkm(dds)
```
3. FPKM values were merged across species using ortholog-aligned gene IDs.
4. Values were transformed as `log2(FPKM + 1)` and saved as `*_Count.tsv`.
Expression values were transformed as:

```r
log2(FPKM + 1)
```

---

## DEG identification

Differential expression was evaluated with a linear model per gene:

```r
lm_fit <- lm(Expression ~ Group + Species:Sequencer, data = data_to_plot)
```

- **Group**: sleep phenotype category (e.g. short vs long sleep).
- **Species:Sequencer**: interaction term to control species-specific sequencing platform effects.

For each gene, the coefficient of `Group` and its p-value were extracted; FDR was applied (`p.adjust(..., method = "BH")`). Genes with FDR < 0.05 are reported in `*_DEGs.tsv`.

---

## Statistical rationale

- **DESeq2** was used only for normalization.
- **Linear modeling** was used for phenotype association.
- Sequencing platform effects are controlled via the species–sequencer interaction.
- Multiple testing correction: Benjamini–Hochberg FDR.

This setup supports cross-species transcriptomic association analysis for sleep phenotypes.
