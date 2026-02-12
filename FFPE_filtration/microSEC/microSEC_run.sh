#!/bin/bash
#SBATCH --partition=k40 # Adjust
#SBATCH --output=/path/to/log_files/slurm-microsec-%A_%a.out # Adjust
#SBATCH --array=1-24%12  # Adjust array range and concurrency as needed

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------
output_dir="/path/to/output_directory"
r_script="/path/to/MicroSEC.R"
info_file_list="/path/to/info_files_list.txt"

# -----------------------------
# START LOG
# -----------------------------
echo "Job started"
date "+%Y-%m-%d %H:%M:%S"

# -----------------------------
# GET FILE FOR THIS ARRAY TASK
# -----------------------------
# Extract the file corresponding to this SLURM_ARRAY_TASK_ID
arrayfile=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $1}' "$info_file_list")
echo "Processing file: $arrayfile"

# -----------------------------
# RUN R SCRIPT
# -----------------------------
Rscript "$r_script" "$output_dir" "$arrayfile" Y

# -----------------------------
# END LOG
# -----------------------------
echo "Job finished"
date "+%Y-%m-%d %H:%M:%S"