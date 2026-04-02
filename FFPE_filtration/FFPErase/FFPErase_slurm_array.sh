#!/bin/bash
#SBATCH --job-name=FFPErase_FFPE
#SBATCH --array=1-100
#SBATCH --cpus-per-task=8


# ==========================================================
# FFPErase SLURM Array Pipeline — FFPE Artifact Filtering
# ==========================================================
#
# DESCRIPTION
# This pipeline runs FFPErase to identify and correct FFPE-induced
# sequencing artifacts. Each SLURM array task processes one sample
# from a TSV table.
#
# INPUT FILE FORMAT (TSV)
# Columns:
# DATASET SAMPLE VCF BAM INSERT COVERAGE
#
# OUTPUT
# Output directory per sample:
# OUTBASE/DATASET/SAMPLE
# The pipeline produces filtered VCFs and reports in this directory.
#
# TEMPORARY FILES
# Work directories are created under WORKBASE/DATASET/SAMPLE
# and removed after successful completion.
#
# SOFTWARE REQUIREMENTS
# • Nextflow
# • Conda environment with FFPErase
#
# USAGE
# 1. Adjust --array to the number of samples in the table
# 2. Submit the job:
#       sbatch FFPErase_slurm_array.sh
# ==========================================================

set -euo pipefail

# ===============================
# CONFIGURATION
# ===============================

TABLE="sample_details_v1-CGCI-Tumor.tsv"

REFERENCE="/library/GENOMES/hg38/Homo_sapiens_assembly38.fasta"
BED="/library/GENOMES/hg38/wgs_calling_regions.hg38.bed"
MODEL="/library/FFPErase/model.snvs.pkl"

OUTBASE="FFPErase"
WORKBASE="FFPErase/workdirs"

PIPELINE="~/.nextflow/assets/papaemmelab/nf-ffperase"

# ===============================
# LOAD ENVIRONMENT
# ===============================

conda activate nextflow

# ===============================
# READ SAMPLE LINE FOR ARRAY TASK
# ===============================

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ${TABLE})

IFS=$'\t' read -r DATASET SAMPLE VCF BAM INSERT COVERAGE <<< "${LINE}"

# ===============================
# SANITY CHECK
# ===============================

echo "===================================="
echo "Running sample: ${SAMPLE}"
echo "Dataset: ${DATASET}"
echo "VCF: ${VCF}"
echo "BAM: ${BAM}"
echo "Insert size: ${INSERT}"
echo "Coverage: ${COVERAGE}"
echo "===================================="

# ===============================
# CREATE OUTPUT / WORK DIRECTORIES
# ===============================

OUTDIR="${OUTBASE}/${DATASET}/${SAMPLE}"
WORKDIR="${WORKBASE}/${DATASET}/${SAMPLE}"

mkdir -p "${OUTDIR}"
mkdir -p "${WORKDIR}"

cd "${OUTDIR}" || exit 1

# Nextflow configuration for local python links
cat > nextflow.config << 'EOF'
process {
    beforeScript = '''
    ln -sf /usr/bin/python3 ./python
    ln -sf /usr/bin/python3 ./python33
    export PATH=$(pwd):$PATH
    '''
}
EOF

# ===============================
# NEXTFLOW HOME (per-sample)
# ===============================

NXF_BASE="${OUTBASE}/nxf_home"
NXF_HOME="${NXF_BASE}/${DATASET}/${SAMPLE}"
mkdir -p "${NXF_HOME}"
export NXF_HOME

# ===============================
# RUN FFPErase PIPELINE
# ===============================

nextflow run "${PIPELINE}" \
    --step full \
    --vcf "${VCF}" \
    --bam "${BAM}" \
    --reference "${REFERENCE}" \
    --outdir "${OUTDIR}" \
    --coverage "${COVERAGE}" \
    --medianInsert "${INSERT}" \
    --mutationType snvs \
    --bed "${BED}" \
    --model "${MODEL}" \
    --modelName snvs \
    -work-dir "${WORKDIR}" \
    -resume

EXIT_CODE=$?

# ===============================
# CLEANUP WORKDIR
# ===============================

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Pipeline finished successfully"
    echo "Removing work directory: ${WORKDIR}"
    rm -rf "${WORKDIR}"
else
    echo "ERROR: Pipeline failed (exit code ${EXIT_CODE})"
    echo "Work directory preserved for debugging: ${WORKDIR}"
fi

exit ${EXIT_CODE}