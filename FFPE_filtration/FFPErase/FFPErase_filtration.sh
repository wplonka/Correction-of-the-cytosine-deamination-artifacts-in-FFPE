#!/bin/bash
# ============================================================
# FFPErase VCF Filtering Pipeline
# ------------------------------------------------------------
# This script:
# 1. Loops over FFPErase TSV files
# 2. Matches corresponding VCF files
# 3. Filters variants according to TSV "True" labels
# 4. Saves filtered VCFs and a summary table
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Folder containing FFPErase TSV files
FFPErase_FOLDER="/path/to/FFPErase_tsv"

# Folder containing original VCF files
VCF_FOLDER="/path/to/original_vcf"

# Output folder for filtered VCFs and summary
OUT_FOLDER="/path/to/output_folder"

# -----------------------------
# SETUP
# -----------------------------

mkdir -p "$OUT_FOLDER"

SUMMARY="$OUT_FOLDER/FFPErase_filtration_summary.txt"
LOG="$OUT_FOLDER/filtering_log.txt"

echo -e "Sample\tVariants_before\tVariants_after\tTSV_True\tTSV_False" > "$SUMMARY"

# -----------------------------
# LOOP OVER TSV FILES
# -----------------------------

for tsv in "$FFPErase_FOLDER"/*.classified_df_snvs.tsv
do
    base=$(basename "$tsv" .classified_df_snvs.tsv)
    echo "Processing $base ..." | tee -a "$LOG"

    # Find matching VCF
    vcf=$(ls "$VCF_FOLDER"/"$base".Mutect2*.vcf 2>/dev/null)
    if [ ! -f "$vcf" ]; then
        echo "  ERROR: No matching VCF found for $base" | tee -a "$LOG"
        continue
    fi

    echo "  TSV: $tsv" | tee -a "$LOG"
    echo "  VCF: $vcf" | tee -a "$LOG"

    # Count True / False in TSV
    tsv_true=$(awk -F'\t' '$NF=="True"{c++} END{print c+0}' "$tsv")
    tsv_false=$(awk -F'\t' '$NF=="False"{c++} END{print c+0}' "$tsv")

    echo "  TSV counts: True=$tsv_true, False=$tsv_false" | tee -a "$LOG"

    # Prepare keep list
    keep_file="$OUT_FOLDER/$base.keep"
    awk -F'\t' -v strip_chr="$STRIP_CHR_PREFIX" '
    NR>1 && $NF=="True" {
        chr=$1
        if(strip_chr=="true") sub(/^chr/,"",chr)
        print chr"\t"$2"\t"$4"\t"$5
    }' "$tsv" > "$keep_file"

    echo "  Keep list saved to $keep_file (first 5 lines):" | tee -a "$LOG"
    head -n5 "$keep_file" | tee -a "$LOG"

    # Count variants before
    before=$(grep -vc "^#" "$vcf")
    echo "  Variants before filtering: $before" | tee -a "$LOG"

    # Filter VCF
    filtered_vcf="$OUT_FOLDER/${base}_FFPErase_filtered.vcf"
    awk -v keep="$keep_file" -v '
    BEGIN{
        FS=OFS="\t"
        while((getline<keep)>0){
            key=$1 FS $2 FS $3 FS $4
            K[key]=1
        }
    }
    /^#/ { print; next }
    {
        chr=$1
        if(strip_chr=="true") sub(/^chr/,"",chr)
        pos=$2
        ref=$4
        alt=$5
        key=chr FS pos FS ref FS alt
        if(key in K) print
    }
    ' "$vcf" > "$filtered_vcf"

    # Count variants after
    after=$(grep -vc "^#" "$filtered_vcf")
    echo "  Variants after filtering: $after" | tee -a "$LOG"

    # Append summary
    echo -e "$base\t$before\t$after\t$tsv_true\t$tsv_false" >> "$SUMMARY"

    # Cleanup
    rm "$keep_file"
    echo "  Done $base" | tee -a "$LOG"

done

echo "FFPErase filtering finished at $(date)" | tee -a "$LOG"
echo "Summary saved to $SUMMARY" | tee -a "$LOG"











FFPErase_folder="/scratch/scratch-ssd/rjaksik/FFPErase/SUT_WES"
VCF_folder="/scratch/scratch-hdd/4all/wplonka/FFPE_results/LUAD/SEP/original"
OUT_folder="/scratch/scratch-hdd/4all/wplonka/FFPE_results/LUAD/SEP/FFPErase2"


#!/bin/bash

# ----------------------------
# Folders
# ----------------------------

mkdir -p "$OUT_folder"

SUMMARY="$OUT_folder/LUAD_FFPErase_filtration_summary_v2.txt"
LOG="$OUT_folder/filtering_log.txt"

echo "FFPErase filtering started at $(date)" > "$LOG"
echo -e "name\tbefore\tafter\ttsv_true\ttsv_false" > "$SUMMARY"

# ----------------------------
# Loop over TSV files
# ----------------------------
for tsv in "$FFPErase_folder"/*.classified_df_snvs.tsv
do
    base=$(basename "$tsv" .classified_df_snvs.tsv)
    echo "Processing $base ..." | tee -a "$LOG"

    # Find matching VCF
    vcf=$(ls "$VCF_folder"/"$base".Mutect2*.vcf 2>/dev/null)
    if [ ! -f "$vcf" ]; then
        echo "  ERROR: No matching VCF found for $base" | tee -a "$LOG"
        continue
    fi

    echo "  TSV: $tsv" | tee -a "$LOG"
    echo "  VCF: $vcf" | tee -a "$LOG"

    # Count True / False in TSV
    tsv_true=$(awk -F'\t' '$NF=="True"{c++} END{print c+0}' "$tsv")
    tsv_false=$(awk -F'\t' '$NF=="False"{c++} END{print c+0}' "$tsv")

    echo "  TSV counts: True=$tsv_true, False=$tsv_false" | tee -a "$LOG"

    # Prepare keep list
    keep_file="$OUT_folder/$base.keep"
    awk -F'\t' '
    NR>1 && $NF=="True" {
        chr=$1
        sub(/^chr/,"",chr)
        print chr"\t"$2"\t"$4"\t"$5
    }' "$tsv" > "$keep_file"

    echo "  Keep list saved to $keep_file (first 5 lines):" | tee -a "$LOG"
    head -n5 "$keep_file" | tee -a "$LOG"

    # Count variants before
    before=$(grep -vc "^#" "$vcf")
    echo "  Variants before filtering: $before" | tee -a "$LOG"

    # Filter VCF
    filtered_vcf="$OUT_folder/${base}_FFPErase_filtered.vcf"
    awk -v keep="$keep_file" '
    BEGIN{
        FS=OFS="\t"
        while((getline<keep)>0){
            key=$1 FS $2 FS $3 FS $4
            K[key]=1
        }
    }
    /^#/ { print; next }
    {
        chr=$1
        sub(/^chr/,"",chr)
        pos=$2
        ref=$4
        alt=$5
        key=chr FS pos FS ref FS alt
        if(key in K) print
    }
    ' "$vcf" > "$filtered_vcf"

    # Count variants after
    after=$(grep -vc "^#" "$filtered_vcf")
    echo "  Variants after filtering: $after" | tee -a "$LOG"

    # Append summary
    echo -e "$base\t$before\t$after\t$tsv_true\t$tsv_false" >> "$SUMMARY"

    # Cleanup
    rm "$keep_file"
    echo "  Done $base" | tee -a "$LOG"
done

echo "Summary saved to $SUMMARY" | tee -a "$LOG"
