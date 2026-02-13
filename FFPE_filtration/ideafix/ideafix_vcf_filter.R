# ============================================================
# ideafix VCF Filtering Pipeline
# ------------------------------------------------------------
# This script:
# 1. Reads VCF files (Mutect2 pass)
# 2. Extracts variant descriptors using ideafix
# 3. Classifies variants with XGBoost model
# 4. Removes predicted deaminations
# 5. Writes filtered VCF files
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Input directory containing VCF files
INPUT_VCF_DIR <- "path/to/your/vcf_directory"

# Pattern to match VCF files (used in list.files())
VCF_PATTERN <- "pass\\.vcf$"

# Reference genome (FASTA)
REFERENCE_FASTA <- "path/to/reference_genome.fa"

# Output directory
OUTPUT_DIR <- "path/to/output_directory"


# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Optional: output folder for TSV annotations
ANNOTATION_OUTPUT_DIR <- OUTPUT_DIR

# Overwrite existing output files?
OVERWRITE <- FALSE

# -----------------------------
# LOAD REQUIRED LIBRARIES
# -----------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ideafix)
library(VariantAnnotation)

# -----------------------------
# FIX one_hot_encoding BUG IN ideafix
# -----------------------------

one_hot_encoding2 <- function(variant_descriptors) {
  categorical_vars <- variant_descriptors %>%
    select_if(function(x) is.factor(x) | is.character(x)) %>%
    colnames()
  categorical_vars <- setdiff(categorical_vars, c("id", "complete_id"))
  if (length(categorical_vars)) {
    dedup_X_by_var <- lapply(categorical_vars, function(var) {
      variant_descriptors %>%
        mutate(yesno = 1) %>%
        spread(var, yesno, fill = 0, sep = "_")
    })
    dedup_X_by_var %>%
      purrr::reduce(.f = full_join) -> dedup_X
    # This final step discards id/complete_id column
    dedup_X <- dedup_X %>%
      select_if(is.numeric)     
    
    #ideafix-FIX: fill in the missing columns
    XGBoost_model_filename <- system.file("extdata", "XGBoost_final_model.RDS", package = "ideafix")
    XGBoost_model <- readRDS(XGBoost_model_filename)
    expected_cols  = colnames(XGBoost_model$trainingData)
    missing_cols = expected_cols[!expected_cols %in% colnames(dedup_X)]
    dedup_X_mod <- dedup_X %>%
      bind_cols(tibble(!!!setNames(rep(list(0), length(missing_cols)), missing_cols))) %>%
      select(all_of(expected_cols))
    
    return(dedup_X_mod)
  } else {
    return(dplyr::select(variant_descriptors, -ends_with("id")))
  }
}

environment(one_hot_encoding2) <- asNamespace('ideafix')
assignInNamespace("one_hot_encoding", one_hot_encoding2, ns = "ideafix")

# -----------------------------
# PREPARE INPUT FILE LIST
# -----------------------------

vcf_files <- list.files(
  INPUT_VCF_DIR,
  full.names = TRUE,
  pattern = VCF_PATTERN
)

vcf_names <- basename(vcf_files)

if (length(vcf_files) == 0) {
  stop("No VCF files found matching pattern.")
}

# -----------------------------
# MAIN PROCESSING LOOP
# -----------------------------

for (i in seq_along(vcf_files)) {
  
  tryCatch({
    
    message("Processing: ", vcf_names[i])
    
    sample_name <- sub(".Mutect2.filter.pass.vcf", "", vcf_names[i])
    output_vcf <- file.path(OUTPUT_DIR, paste0(sample_name, "_ideafix_filtered.vcf"))
    
    if (file.exists(output_vcf) && !OVERWRITE) {
      message("Output exists. Skipping: ", output_vcf)
      next
    }
    
    # -----------------------------
    # Extract variant descriptors
    # -----------------------------
    
    descriptors <- get_descriptors(
      vcf_filename = vcf_files[i],
      fasta_filename = REFERENCE_FASTA
    )
    
    descriptors <- descriptors[!is.infinite(descriptors$base.qual.frac), ]
    
    # -----------------------------
    # Classify variants (XGBoost)
    # -----------------------------
    
    predictions <- classify_variants(
      variant_descriptors = descriptors,
      algorithm = "XGBoost"
    )
    
    # Optional: write TSV annotations
    annotate_deaminations(
      predictions,
      format = "tsv",
      outfolder = ANNOTATION_OUTPUT_DIR,
      outname = paste0("ideafix_labels_", sample_name)
    )
    
    # -----------------------------
    # Load VCF and filter variants
    # -----------------------------
    
    vcf <- readVcf(vcf_files[i], genome = "hg38")
    
    to_remove <- predictions[predictions$DEAMINATION == "deamination",]
    
    vcf_chr <- as.character(seqnames(vcf))
    vcf_pos <- start(vcf)
    vcf_ref <- as.character(ref(vcf))
    vcf_alt <- sapply(alt(vcf), function(x) as.character(x[1]))
    
    mask <- !mapply(function(chr, pos, ref, alt) {
      any(
        to_remove$CHROM == chr &
          to_remove$POS == pos &
          to_remove$REF == ref &
          to_remove$ALT == alt
      )
    }, vcf_chr, vcf_pos, vcf_ref, vcf_alt)
    
    filtered_vcf <- vcf[mask]
    
    writeVcf(filtered_vcf, output_vcf)
    
    message("Finished: ", sample_name)
    
  }, error = function(e) {
    
    message("Error processing: ", vcf_files[i])
    message("Error message: ", e$message)
    
    write(
      paste(Sys.time(), "-", vcf_files[i], "-", e$message),
      file = file.path(OUTPUT_DIR, "error_log.txt"),
      append = TRUE
    )
  })
}

# -----------------------------
# DONE
# -----------------------------

message("All files processed.")

