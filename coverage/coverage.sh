#!/bin/bash
#SBATCH --partition=your_partition
#SBATCH --array=0-9%5   # Adjust according to number of tasks and concurrency

# ============================================================
# BED Coverage Extraction per Sample
# ------------------------------------------------------------
# This script:
# 1. Converts merged TSV files to BED format
# 2. Sorts and filters BED
# 3. Calculates per-sample coverage using bedtools coverage
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------
TSV_FOLDER="path/to/merged_tsvs"
OUTPUT_FOLDER="path/to/bed_output"
COVERAGE_FOLDER="path/to/coverage_output"
BAM_FOLDER="path/to/bam_files"

mkdir -p "$OUTPUT_FOLDER"
mkdir -p "$COVERAGE_FOLDER"

# -----------------------------
# SELECT TSV FILE FOR THIS SLURM ARRAY TASK
# -----------------------------
TSV_FILES=($TSV_FOLDER/*.tsv)
tsv_file="${TSV_FILES[$SLURM_ARRAY_TASK_ID]}"
tsv_name=$(basename "$tsv_file" .tsv)
echo "Processing TSV: $tsv_file"

# -----------------------------
# CONVERT TSV TO BED
# -----------------------------
bed_file="$OUTPUT_FOLDER/all_${tsv_name}.bed"

if [[ ! -f "$bed_file" ]]; then
    awk 'BEGIN{OFS="\t"} NR>1 {
        ref_length = length($4);
        alt_length = length($5);
        start = $3 - 1;
        end = start + ref_length;
        if (start == end) end += 1;
        print $2, start, end;
    }' "$tsv_file" | sort -u -k1,1 -k2,2n | sort -k1,1V -k2,2n > "${bed_file}.sorted"
    mv "${bed_file}.sorted" "$bed_file"
    echo "  Saved sorted BED: $bed_file"
else
    echo "  BED file already exists, skipping: $bed_file"
fi

# -----------------------------
# EXTRACT SAMPLE NAMES
# -----------------------------
sample_names=$(awk 'NR>1 {print $1}' "$tsv_file" | sort | uniq)
echo "  Processing samples: $sample_names"

# -----------------------------
# CALCULATE COVERAGE FOR EACH BAM
# -----------------------------
for bam_file in "$BAM_FOLDER"/*.bam; do
    bam_name=$(basename "$bam_file" | sed 's/\.bam$//')  # adjust if needed

    if echo "$sample_names" | grep -q "^$bam_name$"; then
        echo " Matched BAM: $bam_file to sample $bam_name"
        temp_coverage_file="$COVERAGE_FOLDER/${bam_name}_${tsv_name}_coverage.txt"

        if [[ -f "$temp_coverage_file" ]]; then
            echo " Coverage file exists, skipping: $temp_coverage_file"
            continue
        fi

        echo " Running bedtools coverage: $bed_file -> $temp_coverage_file"
        bedtools coverage -sorted -a "$bed_file" -b "$bam_file" > "$temp_coverage_file"
        echo "  Coverage saved: $temp_coverage_file"
    fi
done

echo "All tasks completed for TSV: $tsv_file"