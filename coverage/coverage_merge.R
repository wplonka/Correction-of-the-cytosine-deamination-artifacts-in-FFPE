# ============================================================
# Coverage Aggregation Pipeline
# ------------------------------------------------------------
# This script:
# 1. Reads merged variant TSVs generated from VCFs
# 2. Adds per-sample coverage from bedtools output
# 3. Optionally generates boxplots of DP and coverage
# 4. Saves final merged coverage tables
# ============================================================

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------
library(tidyr)
library(ggplot2)

# Paths (edit to your environment)
coverage_folder <- "path/to/folder_with_coverage_files/"
merged_folder   <- "path/to/folder_with_som.py_merged_tsvs/"
output_path     <- "path/to/final_coverage_output_dir/"

# Create output directory if it does not exist
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Option to generate plots (TRUE/FALSE)
GENERATE_PLOTS <- FALSE

# -----------------------------
# LIST FILES
# -----------------------------
coverage_files <- list.files(coverage_folder, pattern = ".txt")
merged_files   <- list.files(merged_folder)

# -----------------------------
# MAIN LOOP
# -----------------------------
for(file in merged_files){
  
  # Read merged variant TSV
  merged <- read.table(file.path(merged_folder, file), header = TRUE)
  
  # Match coverage files based on sample prefix
  cov_files <- coverage_files[grepl(substr(sub("^results_", "", sub("_FFpolish_labeled.tsv","",file)),1,12), coverage_files)] #edit here
  
  # Add BED coordinates
  merged$start <- merged$pozycja - 1
  merged$end   <- merged$pozycja - 1 + nchar(merged$ref)
  
  # Merge coverage columns
  if(length(cov_files) > 0){
    for(cov_file in cov_files){
      coverage <- read.table(file.path(coverage_folder, cov_file))
      coverage <- coverage[,1:4]  # columns: chrom, start, end, coverage
      colnames(coverage)[1:4] <- c("chromosome","start","end", sub("_coverage.txt","",cov_file))
      merged <- merge(merged, coverage, by=c("chromosome","start","end"), all=FALSE)
    }
  }
  
  # Get unique samples
  samples <- unique(merged$sample_name)
  
  # Optional: generate boxplots
  if(GENERATE_PLOTS){
    for(sam in samples){
      temp <- merged[merged$sample_name == sam,]
      
      # Select DP and coverage columns
      coverage_cols <- setdiff(colnames(temp), 
                               c("sample_name","chromosome","position","ref","alt","VAF",
                                 "ALT_COUNT","REF_COUNT","label","start","end"))      
      to_plot <- temp[, c("DP", coverage_cols)]
      
      # Pivot longer for plotting
      to_plot_long <- pivot_longer(to_plot,
                                   cols = everything(),
                                   names_to = "Category",
                                   values_to = "Value")
      
      # Create boxplot
      plot <- ggplot(to_plot_long, aes(x = Category, y = Value, fill = Category)) +
        geom_boxplot() +
        theme_minimal() +
        scale_y_log10() +
        labs(title=paste0("Boxplot of DP and coverage: ", sam),
             x="", y="Value")
      
      ggsave(filename = file.path(output_path, paste0(sam, "_", sub(".tsv","",file), "_coverage.png")),
             plot = plot, width = 8, height = 6)
    }
  }
      
  