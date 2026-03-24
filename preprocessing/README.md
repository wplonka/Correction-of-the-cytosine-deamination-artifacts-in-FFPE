# mutation_study.py

## Usage

mutation_study.py samples_table pair_nr [--mc=mt2] [--ctrl=name] [--refgenome=GRCh37|GRCm38|GRCz11|hg38] [--pr=N] [--novaseq] [--smallbam] [-s]


---

## Description

Identifies somatic mutations from WGS/WES data starting from FASTQ files.

---

## Input

**samples_table**  
Tab-delimited table containing samples and associated files.

Accepted format:
patient_ID sample_ID Lane R1.fa [R2.fa] [Sex:M/F]


**pair_nr**  
Sample number to be analyzed (intended for use with SLURM job arrays or GNU parallel).

---

## Output

MuTect output files:

- VCF files containing detected variants  
- Coverage files

---

## Parameters

**--ctrl=name**  
Name of the control sample for somatic calling  
Default: none

**--mc=mt2 | none**  
Mutation caller type  
Default: mt2 (MuTect2)

**--mcpon=file_name | none**  
Panel of Normals (PoN) file used with MuTect2  
Default: none

**--refgenome=GRCh37 | GRCm38 | GRCz11 | hg38**  
Reference genome version  
Default: hg38

**--pr=N**  
Number of processors used per task  
Default: maximum available

**--smallbam**  
Removes the original BAM file and keeps only a reduced BAM containing reads overlapping variants detected by MuTect2

**--novaseq**  
Enables additional quality trimming in TrimGalore adapted to 2-color chemistry used in NovaSeq/NextSeq platforms (ignores quality scores of G bases)

**--adapter1=SEQUENCE**  
Custom 3' adapter sequence  
Default: auto-detection by TrimGalore

**--adapter2=SEQUENCE**  
Custom 5' adapter sequence  
Default: auto-detection by TrimGalore

**-s**  
Simulation mode — prints commands without executing them  
Default: OFF

---

## Requirements

### HTS Processing Tools

- bedtools  
- samtools  
- bwa  
- GATK — Genome Analysis Toolkit
