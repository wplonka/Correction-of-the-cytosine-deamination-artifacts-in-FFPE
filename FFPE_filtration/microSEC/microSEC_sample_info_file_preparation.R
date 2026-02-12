# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# directories 
vcf_dir <- "path/to/vcf_directory"
bam_dir <- "path/to/bam_directory"
output_dir <- "path/to/output_info_files"
tsv_dir <- "path/to/microsec_mutation_information_tsv_directory"

# reference files
reference_fasta <- "path/to/reference_genome.fa"
simple_repeat <- "path/to/simpleRepeat.bed"

# fixed parameters
adaptor1_seq <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adaptor2_seq <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
organism_name <- "hg38"
analysis_type <- "TOP"
read_length <- 150 #or load file with this information


# -----------------------------
# CREATE OUTPUT DIR IF NEEDED
# -----------------------------

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# GET VCF FILES
# -----------------------------

files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = FALSE)

if (length(files) == 0) {
  stop("No VCF files found in: ", vcf_dir)
}

# Extract sample names
sample_names <- sub("\\.Mutect2\\.filter\\.pass\\.ffpe\\.vcf$", "", files)


# Construct full paths to input files
tsv_paths <- file.path(tsv_dir, paste0(sample_names, ".tsv"))
bam_paths <- file.path(bam_dir, paste0(sample_names, ".dedup.recal.mutect2.bam"))

read_lengths <- rep(read_length, length(sample_names))

# -----------------------------
# BUILD FINAL DATA FRAME
# -----------------------------

to_save <- data.frame(
  sample_names = sample_names,
  tsv_paths    = tsv_paths,
  bam_paths    = bam_paths,
  read_lengths = read_lengths,
  adaptor1     = rep(adaptor1_seq, length(sample_names)),
  adaptor2     = rep(adaptor2_seq, length(sample_names)),
  organism     = rep(organism_name, length(sample_names)),
  type         = rep(analysis_type, length(sample_names)),
  reference    = rep(reference_fasta, length(sample_names)),
  simpleRepeat = rep(simple_repeat, length(sample_names)),
  stringsAsFactors = FALSE
)

# -----------------------------
# WRITE INFO FILES
# -----------------------------

for (i in seq_len(nrow(to_save))) {
  
  sample_row <- to_save[i, ]
  
  output_file <- file.path(
    output_dir,
    paste0(sample_row$sample_names, "_info_file.tsv")
  )
  
  write.table(
    sample_row,
    file = output_file,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

}
