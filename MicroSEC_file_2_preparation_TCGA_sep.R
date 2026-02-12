library(openxlsx)

setwd('C:\\Users\\wikto\\OneDrive - Politechnika Śląska\\Pulpit\\FFPE\\som.py\\TCGA_sep\\vcf')
files <- list.files(pattern= "*.vcf")

#BLGSP-71-06-00252_TumorFFPE.dedup.recal.mutect2.bam
sample_names <- sub(".Mutect2.filter.pass.ffpe.vcf", "", files)
tsv_paths <- paste0("/scratch/scratch-hdd/4all/wplonka/microSEC/TCGA_sep/dataframes_microSEC/",sub(".Mutect2.filter.pass.ffpe.vcf", "", files), ".tsv" )
BAM_paths <- paste("/scratch/scratch-hdd/4all/wplonka/microSEC/TCGA_sep/bam/", sub(".Mutect2.filter.pass.ffpe.vcf", "", files), ".dedup.recal.mutect2.bam", sep="")
#read_lengths <- rep(150, length(sample_names))
adaptor1 <- rep("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", length(sample_names))
adaptor2 <- rep("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", length(sample_names))
organism <- rep("hg38", length(sample_names))
type <- rep("TOP", length(sample_names))
ref <- rep("/scratch/scratch-hdd/4all/wplonka/ideafix/hg38.fa")
repead <- rep("/scratch/scratch-hdd/4all/wplonka/microSEC/simpleRepeat.bed")

# read lenghts
lengths <- read.table("C:\\Users\\wikto\\OneDrive - Politechnika Śląska\\Pulpit\\FFPE\\microSEC\\max_read_lengths_TCGA_sep.txt", header = T)
lengths$File <- gsub(".dedup.recal.mutect2.bam", "", lengths$File)

df_ordered <- lengths[match(sample_names, lengths$File), ]
read_lengths <- df_ordered$MaxReadLength



to_save <- data.frame(sample_names, tsv_paths, BAM_paths,
                      read_lengths, adaptor1, adaptor2, organism, type,
                      ref, repead)

setwd("C:\\Users\\wikto\\OneDrive - Politechnika Śląska\\Pulpit\\FFPE\\som.py\\TCGA_sep\\info_files")
for (i in 1:dim(to_save)[1]){
  tsv <- to_save[i,]
  write.table(tsv, paste(tsv$sample_names, "_info_file.tsv", sep=""), row.names = F, col.names = F, sep='\t', eol = "\n")
}
