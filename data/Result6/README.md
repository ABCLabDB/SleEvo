## RNA-seq normalization and phenotype-specific expression matrices

Gene-level read counts were generated using featureCounts and processed
separately for each sleep phenotype:

- Total sleep time
- NREM ratio
- Sleep timing per day
- Sleep frequency

Only species with available phenotype annotations were included
for each phenotype-specific matrix.

---

### 1. Species-specific preprocessing

For each species:

- Raw featureCounts output was loaded.
- Gene IDs were filtered using ortholog mapping to ensure
  cross-species alignment.
- Only samples present in phenotype-specific metadata were retained.
- Gene length information (basepairs) was extracted.

---

### 2. DESeq2-based normalization

A DESeq2 dataset was constructed using a null design (~1),
as no differential testing was performed at this stage:

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

This step performs:

- Library size normalization (size factors)
- Length correction using gene basepair information
- Conversion to FPKM values

---

### 3. Cross-species matrix integration

Species-specific FPKM matrices were inserted into a unified
gene × sample matrix aligned by ortholog ID.

Only samples with valid expression values were retained.

---

### 4. Log transformation

Final expression values were transformed as:

```r
log2(FPKM + 1)
```

This stabilizes variance and reduces skewness for
downstream linear modeling.

---

### Output

For each sleep phenotype:

- Log2-transformed FPKM expression matrix (Gene × Sample)
- Cleaned phenotype-specific metadata (Sample × Covariates)

These matrices are used for downstream regression
and comparative cross-species analyses.
