This repository enables full reproduction of the analyses performed to evaluate methods for filtering FFPE-induced sequencing artifacts and benchmark their performance against matched frozen samples used as the gold standard.


📁 **Repository structure**
 - _data_ - Contains information about data availability and metadata required to reproduce the analysis.
 - _preprocessing_ - raw input data --> VCF files generated with GAT Mutect2
 - _FFPE_filtration_ - Execution of FFPE artifact filtering methods and generation of filtered VCF files.
  Each method is implemented in a separate subdirectory:
    - microSEC
    - DeepOmics
    - DeepOmics_PLUS
    - ideafix
    - FFPolish
    - SOBdetector
    - FFPErase
  
- _evaluation_ - The comparison of filtered FFPE VCF files with the matched frozen-sample VCFs was performed using som.py, and the resulting TP, FP, and FN output VCF files were subsequently merged into a single TSV file containing the variant classification labels together with the coverage extracted from the VCFs.
- _coverage_ - Extraction of coverage values from BAM files for FFPE samples and matched frozen samples. These data are then merged with the TSV files generated during the evaluation step.
- visualization - Scripts used to generate all figures and plots presented in the manuscript.



💻 **Programming languages**

The codebase consists of scripts written in:
- R
- Python
- Bash


⚙️ **Dependencies**

R packages
- VariantAnnotation
- openxlsx
- ideafix
- dplyr
- tidyr
- purrr
- tibble
- BSgenome.Hsapiens.UCSC.hg38 # Or your reference genome
- GenomicRanges
- rtracklayer
- stringr
- Biostrings
- GenomicAlignments
- Rsamtools
- remotes
- GenomeInfoDb
- ggplot2
- tidyverse
- data.table
- ComplexHeatmap
- pROC
- gridExtra
- reshape2
- MutationalPatterns
- sigminer
- deconstructSigs
- microSEC.R file

Command-line Tools
- bcftools
- samtools
- GATK — Genome Analysis Toolkit
- bwa
- Nextflow (≥21.x)

 Python & Scripts
- hap.py
  
Bash Utilities
- Standard Bash commands: awk, sed, head, tail, zcat, mkdir, sort, etc.

📄 **Data availability**

Detailed information on data access is provided in the data directory.



📖 **Citation**

If you use this repository in your research, please cite:

_Manuscript details to be added upon publication._
