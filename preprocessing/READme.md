   mutation_study.py samples_table pair_nr [--mc=mt2] [--ctrl=name] [--refgenome=GRCh37|GRCm38|GRCz11|hg38]  [--pr=N] [--novaseq]  [--smallbam] [-s]

   DESCRIPTION
     Identifies somatic mutations based on WGS/WES data starting from FASTQ file

   INPUT:
     samples_table - tab-delimited table of samples and associated files, acceptable formats:
                     patient_ID sample_ID  Lane  R1.fa  [R2.fa]  [Sex:M/F]
     pair_nr - sample number to be analyzed (to be used with SLURM job array or GNU parallel)

   OUTPUT:
     muTect output files: VCF, coverage

   PARAMETERS:
     --ctrl=name - name of the control sample for somatic calling [default:none]
     --mc=[mt2|none] - mutation caller type:  MuTect2 [default:mt2]
     --mcpon=[file_name|none] - file with the PoN (Panel of Normal) data to be used with Mutect2 [default:none]
     --refgenome=[hg38] - reference genome [default:hg38]
     --pr=N - number of processors to be used in each task [default:max_available]
     --smallbam - remove the original BAM file keeping only a new one that contains reads overlapping the variants detected using Mutect2
     --novaseq - enables an additional quality trimming in TrimGalore which is adapted to the 2 color chemistry used in Novaseq/Nextseq platforms (ignores quality scores of G bases)
     --adapter1/--adapter2=SEQUENCE - use custom 3'/5' adapter sequence (default: auto detect adapter sequences in TrimGalore)
     -s - simulation mode (print the commands only) [default:OFF]


   REQUIREMENTS:
     bedtools, samtools, bwa - HTS data processing software
     muTect - identification of SNVs based on HTS data in BAM format
     gatk - genome analysis toolkit


