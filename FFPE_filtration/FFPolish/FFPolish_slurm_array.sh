#!/bin/bash
#SBATCH --job-name=FFPolish_FFPE
#SBATCH --cpus-per-task=8
#SBATCH --array=1-100

# ==========================================================
# FFPolish FFPE Artifact Filtering — SLURM Array Version
# ==========================================================
#
# DESCRIPTION
# Each SLURM array task processes one sample listed in a TSV file.
# The workflow runs FFPolish to detect FFPE-induced sequencing artifacts,
# annotates retained and removed variants, and produces a final sorted VCF.
#
# IMPORTANT ASSUMPTION
# Input VCF files MUST already be bgzipped (.vcf.gz).
# The script will NOT compress VCF files.
#
# INPUT FILE FORMAT (TSV)
# sample_id    input_vcf.gz    input_bam
#
# OUTPUT
# One file per sample:
# sample_id.FFPolish.filtered.vcf
#
# TEMPORARY FILES
# Written to:
# tmp/
#
# SOFTWARE REQUIREMENTS
# • FFPolish
# • bcftools
# • bedtools
# • Conda environment containing FFPolish
#
# USAGE
# 1. Count number of samples:
#       wc -l samples.tsv
# 2. Adjust array range:
#       #SBATCH --array=1-N
# 3. Submit:
#       sbatch ffpolish_array.slurm
#
# ==========================================================

set -euo pipefail
mkdir -p logs tmp

SAMPLE_LIST="samples.tsv"
FFPOLISH="ffpolish"
BCFTOOLS="bcftools"
BEDTOOLS="bedtools"
GENOME="/path/to/genome.fa"

# Activate conda environment
conda activate ffpolish

# Read sample line corresponding to array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

sample_id=$(echo "$LINE" | awk '{print $1}')
vcf_gz=$(echo "$LINE" | awk '{print $2}')
bam=$(echo "$LINE" | awk '{print $3}')

PREFIX="tmp/${sample_id}"

echo "Processing ${sample_id}"
echo "Input VCF: ${vcf_gz}"
echo "Input BAM: ${bam}"

# 1. Run FFPolish filtering
${FFPOLISH} filter -o . -p ${PREFIX} ${GENOME} ${vcf_gz} ${bam}

# 2. Annotate retained variants
HEADER="##INFO=<ID=${sample_id}.ffpe_filter,Number=1,Type=String,Description=\"FFPolish FFPE filtering\">"

awk -F'\t' 'BEGIN{OFS=FS} !/^#/ {$8=$8";'"${sample_id}"'.ffpe_filter=snv"}1' \
  ${PREFIX}_filtered.vcf | \
  ${BCFTOOLS} annotate --header-line "${HEADER}" -o ${PREFIX}.annot.vcf

# 3. Extract removed variants
${BEDTOOLS} intersect -v -a ${vcf_gz} -b ${PREFIX}.annot.vcf > ${PREFIX}.missing.vcf

# 4. Annotate removed variants
awk -F'\t' 'BEGIN{OFS=FS} !/^#/ {$8=$8";'"${sample_id}"'.ffpe_filter=artifact"}1' \
  ${PREFIX}.missing.vcf > ${PREFIX}.missing.annot.vcf

# 5. Merge retained and removed variants
cat ${PREFIX}.annot.vcf ${PREFIX}.missing.annot.vcf > ${PREFIX}.merged.vcf

# 6. Sort final VCF
${BCFTOOLS} sort -o ${sample_id}.FFPolish.filtered.vcf ${PREFIX}.merged.vcf

echo "Finished ${sample_id}"