#!/bin/bash

# ============================================================
# TP/FN/FP Extraction and Merging Pipeline
# ------------------------------------------------------------
# This script:
# 1. Iterates over all subfolders in a specified input folder
# 2. Looks for a 'tpfn' subfolder containing 4 VCF files
# 3. Extracts sample and variant information from each VCF
# 4. Adds labels (FN, FP, TP1, TP2)
# 5. Merges all TSVs into a single output file per sample
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Input folder containing subfolders with 'tpfn' directories
INPUT_FOLDER="path/to/your/FFPErase2_input"

# Output folder for merged TSVs
OUTPUT_FOLDER="path/to/your/merged_output"
mkdir -p "$OUTPUT_FOLDER"

# Labels for the 4 VCF files (assumes order 0000, 0001, 0002, 0003)
LABELS=("FN" "FP" "TP1" "TP2")

# Suffix to append to output filenames
OUTPUT_SUFFIX="eg. _FFpolish_labeled.tsv"

# -----------------------------
# MAIN LOOP
# -----------------------------
for dir in "$INPUT_FOLDER"/*/; do
    echo "Processing folder: $dir"

    # Check if 'tpfn' subfolder exists
    TPFN_FOLDER="$dir/tpfn"
    if [ -d "$TPFN_FOLDER" ]; then
        echo "Found tpfn folder: $TPFN_FOLDER"

        # Get base name of parent folder and prepare output file name
        BASE_NAME=$(basename "$dir")
        OUTPUT_NAME=$(echo "$BASE_NAME" | sed 's/\.Mutect2.*//')${OUTPUT_SUFFIX}

        # Collect VCF files
        FILES=("$TPFN_FOLDER"/*.vcf.gz)
        if [ ${#FILES[@]} -eq 4 ]; then

            TEMP_FILES=()

            for i in {0..3}; do
                INPUT_FILE="$TPFN_FOLDER/000$i.vcf.gz"
                OUTPUT_TEMP_FILE="$TPFN_FOLDER/000${i}_temp.tsv"
                TEMP_FILES+=("$OUTPUT_TEMP_FILE")

                # Extract sample name using bcftools
                SAMPLE_NAME=$(bcftools query -l "$INPUT_FILE" | head -n 1)

                # Process VCF to TSV with label
                zcat "$INPUT_FILE" | awk -v label="${LABELS[$i]}" -v sample_name="$SAMPLE_NAME" 'BEGIN {
                    OFS="\t"
                    print "sample_name","chromosome","position","ref","alt","VAF","DP","ALT_COUNT","REF_COUNT","label"
                }
                /^#CHROM/ {
                    split($0, header, "\t")
                    for (j = 1; j <= length(header); j++) {
                        if (header[j] == "FORMAT") format_col = j
                        if (j > format_col) sample_col = j  # assume single sample
                    }
                }
                !/^#/ {
                    split($format_col, format_fields, ":")
                    split($sample_col, sample_values, ":")

                    dp="NA"; af="NA"; alt_count="NA"; ref_count="NA"; vaf="NA"
                    for (k in format_fields) {
                        if (format_fields[k] == "DP") dp = sample_values[k]
                        if (format_fields[k] == "AF") af = sample_values[k]
                        if (format_fields[k] == "AD") {
                            split(sample_values[k], ad_values, ",")
                            ref_count = ad_values[1]
                            alt_count = ad_values[2]
                        }
                    }

                    if (dp != "NA" && alt_count != "NA") {
                        vaf = alt_count / dp
                    }

                    print sample_name, $1, $2, $4, $5, vaf, dp, alt_count, ref_count, label
                }' > "$OUTPUT_TEMP_FILE"
            done

            # Merge all labeled TSV files into one
            head -n 1 "${TEMP_FILES[0]}" > "$OUTPUT_FOLDER/$OUTPUT_NAME"
            for file in "${TEMP_FILES[@]}"; do
                tail -n +2 "$file" >> "$OUTPUT_FOLDER/$OUTPUT_NAME"
            done

            echo "Merged TSV saved as $OUTPUT_NAME in $OUTPUT_FOLDER"

            # Remove temporary TSV files
            rm "${TEMP_FILES[@]}"

        else
            echo "Exactly 4 .vcf.gz files not found in $TPFN_FOLDER"
        fi
    else
        echo "tpfn folder does not exist in $dir"
    fi
done

echo "All folders processed."