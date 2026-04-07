## Evaluation pipeline
This folder contains scripts to run som.py on paired samples and aggregate the results.

## Instruction
1. Make a list of paired vcf files (original + filtered)
2. Run Run_som.sh
3. Run merge.sh

## Input
- PAIRS_FILE – CSV/TXT file listing paired VCF paths (format: vcf1,vcf2).
- Reference genome – FASTA file (e.g., hg38.fa).
- VCF files – one per sample.
  
## Output
SOM.py Pipeline
- Organized output directories containing results of som.py per pair.
- Log files per pair, e.g., ffpe_frozen_<sample>.log.
- Files generated here are used in the TP/FN/FP extraction pipeline.

TP/FN/FP Extraction and Merging Pipeline
- Merged TSV files per sample (labeled with FN, FP, TP1, TP2) in OUTPUT_FOLDER.
Each TSV includes columns:
```
sample_name, chromosome, position, ref, alt, VAF, DP, ALT_COUNT, REF_COUNT, label
```
These TSVs serve as input for coverage calculation and further analyses in the *coverage* folder.

## Dependencies
- bash
- awk, sed, head, tail, zcat
- bcftools
- python2

