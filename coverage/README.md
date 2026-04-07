
## Coverage Extraction & Aggregation
This folder provides a two-step pipeline for extracting genomic coverage from BAM files and integrating it with variant-level data derived from Mutect2 VCF files.

## Overview
This pipeline consists of two main components:
 - coverage.sh - 
This step processes merged TSV files (see folder *evaluation*) and computes per-sample coverage using bedtools.
 - coverage_merge.R - 
This step merges coverage data back into variant tables and optionally generates QC plots.

## Input
- TSV files generated from Mutect2 VCFs via som.py
- BAM files

## Output
- BED files with variant coordinates
- Coverage files (*_coverage.txt) per sample
- Merged TSV files with added coverage columns
- (Optional) boxplots of DP and coverage distributions

## Dependencies

Bash pipeline:
- bedtools
- awk
- sort

R pipeline:
- tidyr
- ggplot2
