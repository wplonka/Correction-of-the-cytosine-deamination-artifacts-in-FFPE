# ============================================================
# SOBDetector VCF Artifact Filtering Pipeline
# ------------------------------------------------------------
# Reads VCF files, removes "artifact" variants, writes filtered VCFs
# Generates Excel summary of removed/remaining variants
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Folder containing VCF files after SObdetector
INPUT_FOLDER <- "path/to/SOBdetector_vcf_folder"

# Output folder for filtered VCFs
OUTPUT_DIR <- "path/to/filtered_vcf_folder"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Output Excel summary
SUMMARY_FILE <- file.path(OUTPUT_DIR, "SOBdetector_artifact_summary.xlsx")

# Genome reference
GENOME <- "hg38"

# -----------------------------
# LOAD REQUIRED LIBRARIES
# -----------------------------

library(VariantAnnotation)
library(openxlsx) 

# -----------------------------
# PROCESS VCF FILES
# -----------------------------

vcf_files <- list.files(INPUT_FOLDER, pattern = ".vcf", full.names = TRUE)

artifact_summary <- data.frame(
  Sample = character(), 
  Removed_Artifacts = integer(), 
  Remaining_Variants = integer(), 
  stringsAsFactors = FALSE
)

for (vcf_file in vcf_files) {
  
  message("Processing: ", basename(vcf_file))
  
  # Load VCF
  vcf <- readVcf(vcf_file, "hg38")
  
  # Get sample name(s)
  sample_names <- colnames(geno(vcf)$GT)
  
  for (sample_name in sample_names) {
    
    sample_col_name <- paste0(sample_name, ".ffpe_filter")
    
    # Extract INFO fields
    info <- info(vcf)@listData
    info_names <- names(info)
    
    # Logical mask: keep variants that are not "artifact"
    mask <- data.frame(info[grepl(sample_col_name, info_names)]) != "artifact"
    
    new_vcf <- vcf[mask, ]
    
    # Artifact counts
    removed_artifacts <- sum(!mask)
    remaining_variants <- sum(mask)
    
    artifact_summary <- rbind(
      artifact_summary,
      data.frame(
        Sample = sample_name,
        Removed_Artifacts = removed_artifacts,
        Remaining_Variants = remaining_variants
      )
    )
    
    # Write filtered VCF
    output_file <- file.path(OUTPUT_FOLDER, paste0(sample_name, "_SOBdetector_filtered.vcf"))
    writeVcf(new_vcf, output_file)
  }
  
  message("Finished: ", basename(vcf_file))
}

# -----------------------------
# WRITE SUMMARY TO EXCEL
# -----------------------------

write.xlsx(artifact_summary, SUMMARY_FILE, rowNames = FALSE)

message("All files processed. Summary written to ", SUMMARY_FILE)