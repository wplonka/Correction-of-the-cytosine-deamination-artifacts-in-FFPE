# -----------------------------
# LOAD LIBRARIES
# -----------------------------

library(VariantAnnotation)
library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg38) # Or your reference genome
library(GenomicRanges)
library(rtracklayer)

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# Directories
vcf_dir    <- "path/to/vcf_directory"
output_dir <- "path/to/output_dir"
simple_repeat_bed <- "path/to/simpleRepeat.bed"

# Reference genome
reference_genome <- Hsapiens

# -----------------------------
# CREATE OUTPUT DIR IF NEEDED
# -----------------------------

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# FUNCTION: GET MUTATION TYPE
# -----------------------------

get_mutation_type <- function(vcf) {
  ref <- as.character(ref(vcf))
  alt <- as.character(unlist(alt(vcf)))
  
  mutation_type <- ifelse(nchar(ref) == 1 & nchar(alt) == 1, "snv",
                          ifelse(nchar(ref) < nchar(alt), "ins",
                                 ifelse(nchar(ref) > nchar(alt), "del", "Other")))
  
  return(mutation_type)
}

# -----------------------------
# GET VCF FILES
# -----------------------------

vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)

if (length(vcf_files) == 0) {
  stop("No VCF files found in: ", vcf_dir)
}

# -----------------------------
# CONFIGURATION (EDIT HERE)
# -----------------------------

# If the files are named like "sample1.Mutect2.filter.pass.ffpe.vcf",
# the sample name extracted from this file must match the sample name
# used in the corresponding TSV sample info file.
sample_names <- gsub(".Mutect2.filter.pass.ffpe.vcf", "", basename(vcf_files))

# -----------------------------
# IMPORT SIMPLE REPEATS
# -----------------------------

simple_repeats <- import(simple_repeat_bed)

# -----------------------------
# MAIN LOOP: PROCESS VCF FILES
# -----------------------------

for (i in seq_along(vcf_files)) {
  
  cat("Processing sample", i, ":", sample_names[i], "\n")
  
  # Read VCF
  vcf <- readVcf(vcf_files[i])
  
  # Extract row ranges
  vcf_temp <- as.data.frame(rowRanges(vcf))
  
  # Get mutation types
  Mut_type <- data.frame(mut = paste(vcf_temp$width, "-", get_mutation_type(vcf), sep=""))
  
  # Get neighborhood sequences
  region_starts <- GRanges(seqnames = as.character(vcf_temp$seqnames),
                           ranges = IRanges(start = vcf_temp$start - 20, end = vcf_temp$start - 1))
  sequences_starts <- getSeq(reference_genome, region_starts)
  
  region_ends <- GRanges(seqnames = as.character(vcf_temp$seqnames),
                         ranges = IRanges(start = vcf_temp$end + 1, end = vcf_temp$end + 20))
  sequences_ends <- getSeq(reference_genome, region_ends)
  
  final_sequences <- paste(as.character(sequences_starts),
                           as.character(unlist(alt(vcf))),
                           as.character(sequences_ends), sep = "")
  
  # Check overlap with simple repeats
  mutation_gr <- GRanges(seqnames = as.character(vcf_temp$seqnames),
                         ranges = IRanges(start = vcf_temp$start, end = vcf_temp$end))
  overlap <- findOverlaps(mutation_gr, simple_repeats)
  simple_repeat_flag <- ifelse(seq_along(mutation_gr) %in% queryHits(overlap), "Y", "N")
  
  # Build data frame for this sample
  data_to_save <- data.frame(
    Sample = rep(sample_names[i], nrow(vcf_temp)),
    Mut_type = Mut_type$mut,
    Chr = as.character(vcf_temp$seqnames),
    Pos = as.numeric(vcf_temp$start),
    Ref = vcf_temp$REF,
    Alt = as.character(unlist(alt(vcf))),
    SimpleRepeat_TRF = simple_repeat_flag,
    Neighborhood_sequence = final_sequences,
    stringsAsFactors = FALSE
  )
  
  # Write output to file
  write.table(
    data_to_save,
    file = file.path(output_dir, paste0(sample_names[i], ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
