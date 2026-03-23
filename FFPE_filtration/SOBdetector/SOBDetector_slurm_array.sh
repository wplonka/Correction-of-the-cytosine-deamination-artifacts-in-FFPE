#!/bin/bash
#SBATCH --job-name=SOBDetector_FFPE
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100

# ==========================================================
# SOBDetector FFPE Artifact Filtering — SLURM Array Version
# ==========================================================
#
# Each SLURM task processes ONE sample from samples.tsv
#
# INPUT (TSV):
# sample_id    input_vcf    input_bam
#
# OUTPUT:
# sample_id.SOBDetector.filtered.vcf
#
# USAGE:
# 1. Adjust --array range to number of samples
# 2. sbatch sobdetector_array.slurm
#
# ==========================================================

set -euo pipefail
mkdir -p logs

SAMPLE_LIST="samples.tsv"
SOBDETECTOR="/path/to/sobdetector"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

sample_id=$(echo "$LINE" | awk '{print $1}')
vcf_path=$(echo "$LINE" | awk '{print $2}')
bam_path=$(echo "$LINE" | awk '{print $3}')

temp_vcf="${sample_id}.SOBDetector.tmp"
final_vcf="${sample_id}.SOBDetector.filtered.vcf"

echo "Processing ${sample_id}"

rm -f ${temp_vcf}

${SOBDETECTOR} \
  --input-type VCF \
  --input-variants ${vcf_path} \
  --input-bam ${bam_path} \
  --sample-name ${sample_id} \
  --output-variants ${temp_vcf} \
  --only-passed false

sed -i "s/artiStatus/${sample_id}.ffpe_filter/g" ${temp_vcf}
mv ${temp_vcf} ${final_vcf}

echo "Finished ${sample_id}"