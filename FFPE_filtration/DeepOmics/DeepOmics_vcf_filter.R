# ============================================================
# DeepOmics FFPE VCF Filtering Pipeline
# ------------------------------------------------------------
# This script:
# 1. Reads VCF files
# 2. Uses INFO field (IS_VARIANT) to detect FFPE artifacts
# 3. Removes variants predicted as artifacts
# 4. Writes filtered VCF files
# 5. Creates a summary table (removed vs remaining variants)
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

input_dir     <- "path/to/your/DeepOmics_vcf_directory"
output_dir    <- "path/to/your/output_directory"
summary_file  <- "path/to/your/output_directory/Deepomics_artifact_summary.xlsx"
genome_build  <- "hg38"


# -----------------------------
# LOAD LIBRARIES
# -----------------------------

library(VariantAnnotation)
library(openxlsx)

# -----------------------------
# LIST VCF FILES
# -----------------------------

vcf_files <- list.files(input_dir, pattern = "\\.vcf$", full.names = TRUE)

# -----------------------------
# MAIN LOOP
# -----------------------------

#summary data frame
artifact_summary <- data.frame(
  Sample = character(),
  Removed_Artifacts = integer(),
  Remaining_Variants = integer(),
  stringsAsFactors = FALSE
)

for (vcf_file in vcf_files) {
  
  tryCatch({
    
    # Read vcf
    vcf <- readVcf(vcf_file, genome_build)
    
    sample_name <- colnames(geno(vcf)$GT)
    
    # Extract INFO field
    info_df <- data.frame(IS_VARIANT = info(vcf)$IS_VARIANT)
    

    # Filter artifacts
    filtered_vcf <- vcf[info_df$IS_VARIANT != "artifact", ]
    
    # Count statistics
    removed_artifacts  <- sum(info_df$IS_VARIANT == "artifact", na.rm = TRUE)
    remaining_variants <- sum(info_df$IS_VARIANT != "artifact", na.rm = TRUE)

    artifact_summary <- rbind(
      artifact_summary,
      data.frame(
        Sample = sample_name,
        Removed_Artifacts = removed_artifacts,
        Remaining_Variants = remaining_variants
      )
    )
    
    # Write filtered vcf
    output_file <- file.path(
      output_dir,
      paste0(sample_name, "_DeepOmics_filtered.vcf")
    )
    
    writeVcf(filtered_vcf, output_file)
    
    message("Processed: ", basename(vcf_file))
    
  }, error = function(e) {
    
    message("Error processing: ", basename(vcf_file))
    message("Reason: ", e$message)
    
  })
}


# Write summary table
write.xlsx(artifact_summary, summary_file, rowNames = FALSE)
