# FFPE Artifact Correction - Visualizations

## Overview

This code provides a reproducible **R/Sweave analysis pipeline** for evaluating and comparing methods designed to correct **formalin‑fixed paraffin‑embedded (FFPE) sequencing artifacts** in cancer genomics data.

The workflow integrates variant datasets from multiple cohorts and systematically benchmarks several artifact‑correction tools by assessing their effects on variant calls, mutational signatures, and downstream analytical performance.

The pipeline was developed as part of a large‑scale comparative study of FFPE artifact filtering strategies and is suitable for whole‑genome sequencing (WGS) and whole‑exome sequencing (WES) analyses.


## Output

The pipeline generates:

- Integrated variant datasets
- Performance evaluation tables
- Article figures (Fig.2-4)
- Supplementary tables


## Requirements

### Software

- **R** (≥ 4.1 recommended)

### Required R Packages

Core packages:

- tidyverse
- data.table
- ggplot2
- ComplexHeatmap
- pROC
- gridExtra
- reshape2

Bioinformatics packages:

- VariantAnnotation
- MutationalPatterns
- BSgenome.Hsapiens.UCSC.hg38
- sigminer
- deconstructSigs


## Input Data

The pipeline expects:

- VCF files (raw and artifact‑corrected)
- Variand-specific TSV tables obtained using the FFPE_filtration scripts, e.g.: 

```
| chr  | pos    | ref | alt | start  | end    | sample_name | Tumor.VAF | Tumor.DP | Tumor.ALT_COUNT | Tumor.REF_COUNT | Classification_label | FFPE.DP | FF.DP |
|------|--------|-----|-----|--------|--------|-------------|-----------|----------|-----------------|-----------------|----------------------|---------|-------|
| chr1 | 182120 | C   | A   | 182119 | 182120 | Lu06_P1     | 1         | 2        | 2               | 0               | FN                   | 2       | 8     |
| chr1 | 268716 | T   | C   | 268715 | 268716 | Lu06_P1     | 0.666667  | 4        | 2               | 1               | TP                   | 4       | 2     |
```

---
