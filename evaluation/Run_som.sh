#!/bin/bash

# ============================================================
# SOM.py Pipeline
# ------------------------------------------------------------
# This script:
# 1. Reads pairs of VCF files from a CSV/TXT file
# 2. Runs som.py on each pair
# 3. Saves results and logs in organized directories
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# File containing paired VCF paths (format: vcf1,vcf2)
PAIRS_FILE="path/to/your/paired_files.txt"

# Path to som.py script
SOM_PY="path/to/your/som.py"

# Reference genome fasta
REFERENCE_FA="path/to/your/hg38.fa"

# Output directory for results
OUTPUT_DIR="path/to/your/output_directory"
mkdir -p "$OUTPUT_DIR"

# Log directory
LOG_DIR="path/to/your/log_directory"
mkdir -p "$LOG_DIR"

LOG_PREFIX="eg. ffpe_frozen"


# -----------------------------
# MAIN LOOP
# -----------------------------
while IFS=', ' read -r vcf1 vcf2; do

  # Extract base name from second VCF (Tumor sample) without extension
  base_name=$(basename "$vcf2" ".vcf")

  # Log file for this pair
  log_file="${LOG_DIR}/${LOG_PREFIX}_${base_name}.log"

  # Output prefix for som.py
  output_file="${OUTPUT_DIR}/results_${base_name}_"

  echo "Running som.py on pair: $vcf1 and $vcf2"

  # Run som.py with desired arguments
  python2 "$SOM_PY" "$vcf1" "$vcf2" \
    -o "$output_file" \
    -r "$REFERENCE_FA" \
    --keep-scratch \
    --scratch-prefix "$output_file" \
    -P \
    --include-nonpass \
    --logfile "$log_file" > "$log_file" 2>&1

  echo "Completed analysis for: $vcf1 and $vcf2"

done < "$PAIRS_FILE"

echo "All analyses completed."