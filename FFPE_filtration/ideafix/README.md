## ideafix VCF Filtering Pipeline
This script filters FFPE-induced deamination artifacts from VCF files (Mutect2 PASS variants) using the ideafix package and an XGBoost classification model.

It processes each VCF file by extracting variant features, classifying potential artifacts, and removing predicted deaminations to produce a cleaned VCF for downstream analysis.

## Requirements
- R (≥4.1 recommended)
- bcftools
- samtools
  
## Required R packages:
 - ideafix
 - VariantAnnotation
 - dplyr
 - tidyr
 - purrr
 - tibble

## Input
- VCF files
- FASTA reference genome

## Output
- Filtered VCF files


