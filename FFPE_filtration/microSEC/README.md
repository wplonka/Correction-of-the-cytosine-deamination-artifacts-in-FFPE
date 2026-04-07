## microSEC filtration pipeline
This folder contains a collection of scripts for preparing input data, running microSEC, and filtering sequencing artifacts from somatic variant calls.

## Overview
1. Run microSEC_mutation_info_file_preparation.R
2. Run microSEC_sample_info_file_preparation.R (If you are working with chromosomes separetly edit this code)
3. Make a list of info sample files and save eg. info_files_list.txt
4. Run microSEC_run.sh
   
In case of errors related to incorrect data types in the chromStart and chromEnd columns (e.g., coercion errors, non-numeric values, or issues converting to integer), the following lines should be added to the original microSEC code to explicitly enforce numeric-to-integer conversion:

download_region$chromStart <- as.integer(as.numeric(download_region$chromStart))
download_region$chromEnd   <- as.integer(as.numeric(download_region$chromEnd))

before:
   write_tsv(x = download_region,
              file = paste0(bam_file,".bed"),
              progress = F,
              col_names = F)


5. (Optional) merge results to one tsv if you worked with chromosomes separetly.
6. Make a list of vcf and output tsv from MicroSEC files eg. MicroSEC_paired_files.txt
7. Run microSEC_filtration.sh

## Requirements
- R (≥4.1 recommended)
- samtools

## Required R packages
- VariantAnnotation
- openxlsx
- BSgenome.Hsapiens.UCSC.hg38 # Or your reference genome
- GenomicRanges
- rtracklayer
- stringr
- dplyr
- Biostrings
- GenomicAlignments
- Rsamtools
- remotes
- GenomeInfoDb

## Input
- VCF files
- BAM files
- reference genome FASTA
- simple repead BED
- read lengths
- adapter sequences

## Output
- Filtered VCF files
