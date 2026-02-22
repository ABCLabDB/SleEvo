## Result6: Phenotype-specific RNA-seq matrices and DEG analysis

This directory contains phenotype-specific expression matrices,
metadata tables, and linear model-based DEG results.

---

## File structure

For each sleep phenotype:

- `*_Count.tsv`  
  Log2-transformed FPKM expression matrix (Gene × Sample)

- `*_Metadata.tsv`  
  Sample-level metadata used for modeling

- `*_DEGs.tsv`  
  Differentially expressed genes identified by linear modeling

Phenotypes included:

- Total_sleep_time
- NREM_ratio
- Sleep_frequency
- Sleep_timing

---

## Expression matrix construction

Gene-level counts were generated using featureCounts.
For each species:

1. Raw counts were loaded.
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

3. FPKM values were merged across species
   using ortholog-aligned gene IDs.
4. Expression values were transformed as:

```r
log2(FPKM + 1)
```

The resulting matrix was saved as `*_Count.tsv`.

---

## DEG identification

Differential expression was evaluated using a linear model
applied gene-by-gene:

```r
lm_fit <- lm(Expression ~ Group + Species:Sequencer,
             data = data_to_plot)
```

Model terms:

- **Group**: sleep phenotype category
  (e.g., short vs long sleep)
- **Species:Sequencer**: interaction term to control for
  species-specific sequencing platform effects

For each gene:

- The coefficient of `Group` was extracted
- P-values were obtained from the linear model
- Multiple testing correction was performed using FDR

```r
FDR <- p.adjust(p_value, method = "BH")
```

Genes with FDR < 0.05 were reported in `*_DEGs.tsv`.

---

## Statistical rationale

- DESeq2 was used only for normalization.
- Linear modeling was used for phenotype association testing.
- Sequencing platform effects were controlled via
  species-specific interaction terms.
- Multiple testing correction was applied using
  Benjamini–Hochberg FDR.

This framework enables cross-species transcriptomic
association analysis for sleep phenotypes.
