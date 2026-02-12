#!/bin/bash

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Input file containing paired TSV and VCF paths
input_file="MicroSEC_paired_files.txt"

# Output directory for filtered VCF files
output_dir="path/to/output_directory"

# Output summary file
summary_file="microSEC_filtered_vcf_summary.txt"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Clear the summary file if it exists
> "$summary_file"

# -----------------------------
# PROCESS EACH TSV/VCF PAIR
# -----------------------------
while IFS=$'\t' read -r tsv_file vcf_file; do
    echo "Processing: $tsv_file and $vcf_file"

    # -----------------------------
    # Step 1: Identify rows to remove from TSV
    # -----------------------------
    temp_file=$(mktemp)

    # Assumes that the second column contains the "Artifact suspicious" annotation
    awk 'NR > 1 && $0 ~ /Artifact suspicious/ {print NR-1}' "$tsv_file" > "$temp_file"

    echo "Rows to remove:"
    cat "$temp_file"

    # -----------------------------
    # Step 2: Count variants before filtering
    # -----------------------------
    original_count=$(grep -v "^#" "$vcf_file" | wc -l)
    echo "Original variant count: $original_count"

    # -----------------------------
    # Step 3: Filter VCF based on rows to remove
    # -----------------------------
    intermediate_vcf="${output_dir}/$(basename "${vcf_file%.vcf}_microsEC_filtered.vcf")"

    if [ -s "$temp_file" ]; then
        awk '
        BEGIN {
            while (getline line < "'"$temp_file"'") {
                to_remove[line] = 1;
            }
        }
        $0 ~ /^#/ || !(FNR in to_remove)
        ' "$vcf_file" > "$intermediate_vcf"
    else
        # No rows to remove; copy the original VCF
        cp "$vcf_file" "$intermediate_vcf"
    fi

    # -----------------------------
    # Step 4: Count variants after filtering
    # -----------------------------
    filtered_count=$(grep -v "^#" "$intermediate_vcf" | wc -l)
    echo "Filtered variant count: $filtered_count"

    # -----------------------------
    # Step 5: Save counts to summary file
    # -----------------------------
    echo -e "${vcf_file}\t${original_count}\t${filtered_count}" >> "$summary_file"

    # Clean up temporary file
    rm -f "$temp_file"

done < "$input_file"

echo "Filtering complete. Summary saved to: $summary_file."