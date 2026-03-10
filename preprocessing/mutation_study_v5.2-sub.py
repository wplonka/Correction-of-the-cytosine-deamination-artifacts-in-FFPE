#!/usr/bin/python3

###################################################
#
#  Application: mutation_study - subset of the original code
#       Author: Roman Jaksik
VERSION =       '5.2 - subset (05 Feb 2026)'
#               Silesian University of Technology
#               Gliwice, Poland
#      Contact: roman.jaksik@polsl.pl
#
###################################################


import os, re, math, sys, getopt, time, socket, multiprocessing, subprocess, signal
from subprocess import Popen



def usage():
   print('   ver. '+VERSION)
   print("""
   mutation_study.py samples_table pair_nr [--mc=mt2] [--ctrl=name] [--refgenome=GRCh37|GRCm38|GRCz11|hg38]  [--pr=N] [--novaseq]  [--smallbam] [-s] 

   DESCRIPTION
     Analyzes a set of WXS samples using various mutation callers
    
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
    """)





def main(argv):
  try:
    #processing options and parameters
    opts, args = getopt.gnu_getopt(argv, "sf", ["mc=","ctrl=","cnv=","annot=","gene_annot=","refgenome=","sv=","seqtype=","qc=","pr=","mcpon=","cnvponM=","cnvponF=","crosscheck","startbam","gdc=","ffpe=","sequenzaPatMatch=","phase=","adapter1=","adapter2=","novaseq","smallbam"])
  except:
    #incorrect parameters
    usage()
    printError("Incorrect parameters")    

  if (len(args)!=2):
    #application requires 1 parameter
    usage()
    sys.exit()
  else:

    

    #----starting the script----
    tstart=time.time()
     
    #parameters
    samples_table=args[0]
    pair_nr=int(args[1])

    #flags
    fSimulate=False
    fMutCaller='mt2'
    fCtrlSample=''
    fMutPoN=''
    fGenome='hg38'
    fNovaseq = False 
    fSmallbam = False
    fAdapterSeq1 = ''                    
    fAdapterSeq2 = ''    
    proc=multiprocessing.cpu_count() #number of processors to be used
    for opt,oarg in opts:
       if opt == "-s":
          fSimulate=True;
       if opt == "--mc":
          fMutCaller=oarg;
       if opt == "--ctrl":
          fCtrlSample=oarg;
       if opt == "--refgenome":
          fGenome=oarg;
       if opt == "--pr":
          proc=int(oarg);
       if opt == "--mcpon":
          fMutPoN=oarg;  
       if opt == "--novaseq":
          fNovaseq=True;          
       if opt == "--smallbam":
          fSmallbam=True;   
       if opt == "--adapter1":
          fAdapterSeq1=oarg; 
       if opt == "--adapter2":
          fAdapterSeq2=oarg;           
          

    #check if control sample ID is provided for Mutect2
    if (fMutCaller=='mt2') and (fCtrlSample=='' and fMutPoN==''):
       printError('You have to provide control sample ID through the --ctrl argument or PoN through --mcpon when running Mutect variant discovery')

    #check if selected mutation caller is valid
    if (fMutCaller not in ['mt2','none']):
       printError('The mutation caller that you selected is invalid')



    #################################### SOFTWARE ######################################
          
    #determine on which mashine the script is running
    hostname=socket.gethostname()
    print('\033[0;32mRunning mutation_study ver. %s on %s using %s processors\033[0;37m' %(VERSION,hostname,proc))
    print(sys.argv[1:])

    res_dir     = '/home/rjaksik/library/'
    gatk         = res_dir+'software/gatk-4.6.2.0/gatk4'
    bedtools     = res_dir+'software/bedtools-2.31.1/bin/bedtools'
    samtools     = res_dir+'software/samtools-1.23/samtools' 
    bcftools     = res_dir+'software/bcftools-1.16/bcftools'
    fastq_screen = res_dir+'software/FastQ-Screen-0.15.2/fastq_screen'
    bgzip        = res_dir+'software/htslib-1.15.1/bgzip'  ## from htslib
    bwa          = res_dir+'software/bwa-0.7.17/bwa'
    tabix        = 'tabix'
    trim_galore  = res_dir+'software/TrimGalore-0.6.10/trim_galore'    
    TMP_DIR      = '/scratch/scratch-ssd/rjaksik/tmp'   

 
    
    #################################### DATABASES ######################################   

	####### Human hg38
    if fGenome=='hg38':
       #Reference genome (BWA, BaseRecalibrator, IndelRealigner, MuTect, HaplotypeCaller)
       genome    = res_dir+'GENOMES/hg38/Homo_sapiens_assembly38.fasta'
       calling_reg = res_dir+'GENOMES/hg38/wgs_calling_regions.hg38.interval_list'        
       dbsnp       = res_dir+'SNP/Homo_sapiens_assembly38.dbsnp.vcf.gz'
       millsindels = res_dir+'SNP/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
       indels1kG   = res_dir+'SNP/Homo_sapiens_assembly38.known_indels.vcf.gz'    
       gnomad      = res_dir+'GATK/af-only-gnomad.hg38.vcf.gz'
       exac        = res_dir+'GATK/small_exac_common_3.hg38.vcf.gz'
       polish1k    = res_dir+'SNP/multisample_20210519.dv.bcfnorm.filtered.nogt.ACgt1.vcf.gz'
       phyloP100way= res_dir+'VEP_db/extras/hg38.phyloP100way.bw'
       
 	 
    ########################################################################################   



    #generate the dataset name based on the current folder
    dataset_name = os.path.basename(os.getcwd())


    #select apropriate sample pair
    samppairs=open(samples_table,'r')
    k=0
    pairExists=False
    prev_pat = ''
    cur_pat = ''
    patient = ''
    lFile_R1 =[]
    lFile_R2 =[]
    lFile_sample  =[]
    lFile_lane =[] 
    patient_sex = ''
    for line in samppairs:
       line=line.strip()
       if line and line[0] != '#':
          tline=line.split('\t')
          cur_pat=tline[0]

          #determine the number
          if prev_pat!=cur_pat:
             k=k+1

          #select patient with specific number
          if k==pair_nr:             
             patient = cur_pat

          #select all info for this patient
          if patient==cur_pat:
             lFile_sample.append(tline[1])
             lFile_lane.append(tline[2])
             lFile_R1.append(tline[3])             
             if len(tline)>4:
                lFile_R2.append(tline[4])
             else:
                lFile_R2.append('')
             if len(tline)>5:
                patient_sex=tline[5]
             pairExists=True

          prev_pat = cur_pat 
    samppairs.close()

    lLanes = list(set(lFile_lane))
    lSamples = list(set(lFile_sample))
    Nfiles = len(lFile_R1)
    


    #### used hierarchy:
    ####  - patient (only one is used in the execution of the app - multiple are parallelized on nodes)
    ####    - sample (multiple files can be processed for mutation discovery)
    ####      - lane (multiple lanes will be merged)


    if not(pairExists):
       printError("Sample pair number: "+str(pair_nr)+" does not exists")
    else:
       print('===== Analyzing patient %s out of %s (ID:%s,  Sex:%s,  samples:%s, lanes:%s)' %(pair_nr,k,patient,patient_sex,len(lSamples),len(lLanes)))
       
       ######### DATA PREPROCESSING #########
          
       #perform separatly for two sample files (if neccesary) 
       for i in range(0,len(lSamples)):
          print('== Pre-processing using sample: '+lSamples[i])

          tSample = patient+'_'+lSamples[i]
          logfile = tSample+'.log'
        

          #################### Standard processing from FASTQ level ########################### 
          #skip initial processing if bam files are available
          if (not bamCorrect(samtools,tSample+'.bam')) and (not bamCorrect(samtools,tSample+'.dedup.bam')) and (not bamCorrect(samtools,tSample+'.dedup.recal.bam')) and (not bamCorrect(samtools,tSample+'.dedup.recal.realign.bam')):          
             
             #=== trim adapters using trim_galore             
             #add Novaseq 2-color chemistry trimming                
             if fNovaseq:
                novaseq_command = " --nextseq 20 "                
             else:
                novaseq_command = ""                                
             for k in range(0,Nfiles):
                if lSamples[i]==lFile_sample[k]:
                   #PAIRED-END:
                   if lFile_R2[k]!='':
                      if not os.path.isfile(tSample+'_'+lFile_lane[k]+'_val_1.fq.gz') and not os.path.isfile(tSample+'_'+lFile_lane[k]+'_val_2.fq.gz'):                         
                         
                         if fAdapterSeq1!='':
                            adapter_options = '--adapter %s --adapter2 %s ' %(fAdapterSeq1,fAdapterSeq2)
                         else:
                            adapter_options = ''
                         command = "%s --cores %s --length 20 %s --trim-n --paired --basename %s_%s %s %s %s" % (trim_galore,proc,novaseq_command,tSample,lFile_lane[k],lFile_R1[k],lFile_R2[k],adapter_options) #--retain_unpaired 
                         runTask(command,fSimulate,'Adapter trimming (trim_galore)',logfile)                         

                         
                   #SINGLE-END
                   else:
                      if not os.path.isfile(tSample+'_'+lFile_lane[k]+'_trimmed.fq.gz'):
                         if fAdapterSeq1!='':
                            adapter_options = '--adapter %s' %(fAdapterSeq1)
                         else:
                            adapter_options = ''                                
                         command = "%s --cores %s --length 20 %s --trim-n %s %s" % (trim_galore,proc,novaseq_command,lFile_R1[k],adapter_options)
                         runTask(command,fSimulate,'Adapter trimming (trim_galore)',logfile) 


             #=== align reads using BWA
             for k in range(0,Nfiles):
                if lSamples[i]==lFile_sample[k]:
                   if not os.path.isfile(tSample+'_'+ lFile_lane[k]+'.sam') and not os.path.isfile(tSample+'.us.bam') and not os.path.isfile(tSample+'.bam'):                                      
					  #PAIRED-END:
                      if lFile_R2[k]!='':
                        command_FQ = '%s_%s_val_1.fq.gz %s_%s_val_2.fq.gz' %(tSample,lFile_lane[k],tSample,lFile_lane[k])
                      #SINGLE-END:
                      else:
                        command_FQ = '%s_%s_trimmed.fq.gz' %(tSample,lFile_lane[k])
                      command = "%s mem -M -Y -t %s -R \'@RG\\tID:%s_%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s\'  %s  %s > %s_%s.sam" %(bwa,proc,tSample,lFile_lane[k],tSample,tSample,genome,command_FQ,tSample,lFile_lane[k])
                      runTask(command,fSimulate,'BWA mem alignment',logfile)
                
             #=== combine SAM files into one unsorted BAM
             if not bamCorrect(samtools,tSample+'.us.bam') and not bamCorrect(samtools,tSample+'.bam'):
                if Nfiles>1:
                   command = '%s merge --output-fmt bam %s.us.bam' %(samtools,tSample)
                   for k in range(0,Nfiles):
                      if lSamples[i]==lFile_sample[k]:
                         if bamCorrect(samtools,tSample+'_'+lFile_lane[k]+'.sam') or fSimulate:
                            command = command+' '+tSample+'_'+lFile_lane[k]+'.sam'
                         else:
                            printError('File: %s_%s.sam is empty or does not exist.' %(tSample,lFile_lane[k]))
                else:
                   command = '%s view -S -b %s_%s.sam > %s.us.bam' %(samtools,tSample,lFile_lane[0],tSample)
                runTask(command,fSimulate,'Samtools combine SAM and BAM conversion',logfile)
                     

             #=== check if the previous steps were successful and if yes remove the additional output files
             if bamCorrect(samtools,tSample+'.us.bam'):
                #remove the SAM files
                for k in range(0,Nfiles):
                   if lSamples[i]==lFile_sample[k]:               
                      if os.path.isfile(tSample+'_'+lFile_lane[k]+'.sam'):
                         runTask('rm '+tSample+'_'+lFile_lane[k]+'.sam',fSimulate,'',logfile)

             #=== sort the BAM
             if not bamCorrect(samtools,tSample+'.bam'):    
                command = '%s sort --threads %s %s.us.bam -o %s.bam' %(samtools,proc,tSample,tSample)
                runTask(command,fSimulate,'Samtools BAM sorting',logfile)

             #=== check if the previous steps were successful and if yes remove the additional output files
             if bamCorrect(samtools,tSample+'.bam') and os.path.isfile(tSample+'.us.bam'):
                #remove unsorted bam
                runTask('rm '+tSample+'.us.bam',fSimulate,'',logfile)


              

          if (not bamCorrect(samtools,tSample+'.dedup.bam')) and (not bamCorrect(samtools,tSample+'.dedup.recal.bam')) and (not bamCorrect(samtools,tSample+'.dedup.recal.realign.bam')):          
          
			 #=== GATK:Picard tools MarkDuplicates
             if (not bamCorrect(samtools,tSample+'.dedup.bam')):
                command = '%s MarkDuplicates --INPUT %s.bam  --OUTPUT %s.dedup.bam  --METRICS_FILE %s.dedup_metrics.txt --TMP_DIR %s' %(gatk,tSample,tSample,tSample,TMP_DIR)  #VALIDATION_STRINGENCY=LENIENT
                runTask(command,fSimulate,'GATK:Picard MarkDuplicates',logfile)
                if not os.path.exists('QC/DedupMetrics'):
                   runTask('mkdir QC/DedupMetrics',fSimulate,'',logfile)
                command = 'mv %s.dedup_metrics.txt QC/DedupMetrics' %tSample
                runTask(command,fSimulate,'',logfile)
                
             #=== check if the previous steps were successful and if yes remove the additional output files
             if bamCorrect(samtools,tSample+'.dedup.bam'):
                #remove sorted bam
                runTask('rm '+tSample+'.bam',fSimulate,'',logfile)
                
             #=== index BAM
             if bamCorrect(samtools,tSample+'.dedup.bam'):
                if not os.path.isfile(tSample+'.dedup.bam.bai'):
                   command = '%s index %s.dedup.bam -@ %s' %(samtools,tSample,proc)
                   runTask(command,fSimulate,'Samtools BAM indexing (dedup)',logfile)
             

          if (not bamCorrect(samtools,tSample+'.dedup.recal.bam')) and (not bamCorrect(samtools,tSample+'.dedup.recal.realign.bam')):  

             #=== Recalibrate Bases: https://www.broadinstitute.org/gatk/guide/article?id=2801
             if not os.path.isfile(tSample+'.dedup.recal.bam'):
                command = '%s BaseRecalibrator  -R %s  -I %s.dedup.bam  --known-sites %s  --known-sites %s --known-sites %s -O %s.recal_data.table' %(gatk,genome,tSample,dbsnp,millsindels,indels1kG,tSample)
                runTask(command,fSimulate,'GATK BaseRecalibrator',logfile)
                command = '%s ApplyBQSR  -R %s -I %s.dedup.bam --bqsr-recal-file %s.recal_data.table  -O %s.dedup.recal.bam' %(gatk,genome,tSample,tSample,tSample)                                
                runTask(command,fSimulate,'GATK ApplyBQSR',logfile)
             
             #remove files from previous step
             if bamCorrect(samtools,tSample+'.dedup.recal.bam'):
                runTask('rm '+tSample+'.dedup.bam',fSimulate,'',logfile)
                runTask('rm '+tSample+'.dedup.bam.bai',fSimulate,'',logfile)
                
             #=== index BAM
             if not os.path.isfile(tSample+'.dedup.recal.bam.bai'):
                command = '%s index %s.dedup.recal.bam -@ %s' %(samtools,tSample,proc)
                runTask(command,fSimulate,'Samtools BAM indexing (dedup+recal)',logfile)
          

 
                         

       ######### MUTATION ANALYSIS #########              
       # Option1: paired tumor-normal (somatic)
       #     run with --mt=mt2 and --ctrl=name parameters
       # Option2: tumor vs multisample PoN and optionally normal (somatic)
       #     run with --mt=mt2pon to create *.Mutect2.vcf files for all control samples
       #     create PoN outside of the app using all ctrl samples combined with GATK CreateSomaticPanelOfNormals -vcfs file_list -O out.vcf.gz 
       #     run the app again with --mt=mt2 --mcpon=out.vcf.gz --ctrl=name
       print('== Running mutation analysis: ')

   
       #=== Mutect2 analysis
       if fMutCaller=='mt2':
          mcName = 'Mutect2'
          if not os.path.isfile(patient+'.Mutect2.vcf'):       
             # main command
             command = '%s Mutect2  -R %s  --germline-resource %s  --genotype-germline-sites true --genotype-pon-sites true ' %(gatk,genome,gnomad)
             for samp in lSamples: 
                if samp != fCtrlSample:
                   command += '-I %s_%s.dedup.recal.bam ' %(patient,samp)
             #add PoN
             if fMutPoN!='':
                command += ' --panel-of-normals '+fMutPoN
             #add matching control sample
             if fCtrlSample!='':
                command +=' -I %s_%s.dedup.recal.bam -normal %s_%s' %(patient,fCtrlSample,patient,fCtrlSample) 

             #parallel support
             if proc==1:
                command +=' -L %s -O %s.Mutect2.vcf --f1r2-tar-gz %s_f1r2.tar.gz ' %(calling_reg,patient,patient)
                runTask(command,fSimulate,'GATK Mutect2',logfile)				   
             else:
                   #====Parallell run===
                      #1: split intervals
                      command_split = '%s SplitIntervals -R %s -L %s --scatter-count %s -O %s_INTERVALS' %(gatk,genome,calling_reg,proc,patient)
                      runTask(command_split,fSimulate,'GATK SplitIntervals',logfile)

                      #2: create Mutect2 commands for each I/O file					 
                      mutectTasks = []
                      for task in range(0,int(proc)):
                         interval_file = '%s_INTERVALS/%04d-scattered.interval_list' %(patient,task)                      
                         temp_output = '%s.Mutect2_%04d.vcf' %(patient,task)
                         temp_f1r2 = '%s.Mutect2_f1r2_%04d.tar.gz' %(patient,task)
                         command_task = command+' -L %s -O %s --f1r2-tar-gz %s' %(interval_file,temp_output,temp_f1r2)    
                         if not os.path.isfile(temp_output) or not os.path.isfile(temp_output+'.stats'):
                            mutectTasks.append(command_task)                           

					  #3: run in parallel multiple mutect instances for each interval file					 
                      print('--- Parallel Mutect2 tasks ---\033[0;35m')
                      print(mutectTasks)
                      print('\033[0;37m------------------------------')
                      if not fSimulate:
                         mtstart = time.time()	 
                         processes = [Popen(program, shell=True) for program in mutectTasks]
                         for process in processes:
                            process.wait()    
                         print('\033[0;36mRUNTIME: Parallel Mutect2 - '+runTime(mtstart)+'\033[0;37m')

					  #5: check if all files are provided, if not do not proceed!
                      nFiles = 0
                      for task in range(0,int(proc)):
                         temp_output1 = '%s.Mutect2_%04d.vcf' %(patient,task)   
                         temp_output2 = '%s.Mutect2_%04d.vcf.stats' %(patient,task) 
                         if os.path.isfile(temp_output1) and os.path.isfile(temp_output2):
                            nFiles += 1
                      if nFiles!=int(proc) and not fSimulate:
                         printError('Missing partial result files from the Mutect2 study')

					  #6: combine result files
                      command = '%s MergeVcfs -O %s.Mutect2.vcf' %(gatk,patient)
                      for task in range(0,int(proc)):
                         command += ' -I %s.Mutect2_%04d.vcf' %(patient,task)
                      runTask(command,fSimulate,'',logfile)

					  #7: remove temporary files
                      #Intervals:
                      command = 'rm -r %s_INTERVALS' %(patient)
                      runTask(command,fSimulate,'',logfile)
                      #Partial mutation files
                      if os.path.isfile(patient+'.Mutect2.vcf'):
                         command = 'rm '
                         for task in range(0,int(proc)):
                            command += ' %s.Mutect2_%04d.vcf' %(patient,task)
                            command += ' %s.Mutect2_%04d.vcf.idx' %(patient,task)
                         runTask(command,fSimulate,'',logfile)
                
             #Create read orientation model
             if not os.path.isfile(patient+'.Mutect2_read-orientation-model.tar.gz'):
                command = '%s LearnReadOrientationModel -O %s.Mutect2_read-orientation-model.tar.gz' % (gatk,patient)
                if proc==1:
                   command += '-I %s.Mutect2_f1r2.tar.gz' % (patient)
                else:					                        
                   for task in range(0,int(proc)):
                      command += ' -I %s.Mutect2_f1r2_%04d.tar.gz' %(patient,task)
                runTask(command,fSimulate,'GATK LearnReadOrientationModel',logfile)
             #remove f1r2 files				   
             if os.path.isfile(patient+'.Mutect2_read-orientation-model.tar.gz'): 
                command = ''
                for task in range(0,int(proc)):
                   fname = '%s.Mutect2_f1r2_%04d.tar.gz' %(patient,task)
                   if os.path.isfile(fname):
                      command += ' '+fname
                if command!='':
                   runTask('rm '+command,fSimulate,'',logfile)  

             #Merge the .stats files
             if proc>1:
			 #combine files
                if not os.path.isfile(patient+'.Mutect2.vcf.stats'):       
                   command = '%s MergeMutectStats -O %s.Mutect2.vcf.stats' % (gatk,patient)
                   for task in range(0,int(proc)):
                      command += ' -stats %s.Mutect2_%04d.vcf.stats' %(patient,task)
                   runTask(command,fSimulate,'GATK MergeMutectStats',logfile)
                #remove partial stats files				   
                if os.path.isfile(patient+'.Mutect2.vcf.stats'):
                   command = ''
                   for task in range(0,int(proc)):
                      fname = '%s.Mutect2_%04d.vcf.stats' %(patient,task)
                      if os.path.isfile(fname):
                         command += ' '+fname
                   if command!='':
                      runTask('rm '+command,fSimulate,'',logfile)                      


          #GetPileupSummaries for CalculateContamination          
          # IN: BAM     OUT: pileup.table          
          for samp in lSamples: 
             if not os.path.isfile(patient+'_'+samp+'.pileup.table'):                    
                command = '%s GetPileupSummaries -V %s -L %s -O %s_%s.pileup.table -I %s_%s.dedup.recal.bam' %(gatk,exac,exac,patient,samp,patient,samp)                    
                runTask(command,fSimulate,'GATK GetPileupSummaries: '+samp,logfile) 
          #GetPileupSummaries for CalculateContamination          
          # IN: 2x pileup.table (tumor+normal)     OUT:  contamination.table, segments.table
          for samp in lSamples:
             if samp!=fCtrlSample:   
                if not os.path.isfile(patient+'_'+samp+'.contamination.table'):
                    command = '%s CalculateContamination -I %s_%s.pileup.table -tumor-segmentation %s_%s.segments.table -O %s_%s.contamination.table -matched %s_%s.pileup.table' %(gatk,patient,samp,patient,samp,patient,samp,patient,fCtrlSample)
                    runTask(command,fSimulate,'GATK CalculateContamination',logfile)

          #Filter Mutect2 calls   
          # IN: VCF, read-orientation-model, contamination.table(s), segments.table(s)           OUT: VCF, filterstats
          if not os.path.isfile(patient+'.Mutect2.filter.vcf'):       
             command = '%s FilterMutectCalls -R %s -V %s.Mutect2.vcf -O %s.Mutect2.filter.vcf  --orientation-bias-artifact-priors %s.Mutect2_read-orientation-model.tar.gz --filtering-stats %s.Mutect2.filterstats' %(gatk,genome,patient,patient,patient,patient)
             for samp in lSamples:
                if samp!=fCtrlSample: 
                   command += ' --tumor-segmentation %s_%s.segments.table --contamination-table %s_%s.contamination.table' %(patient,samp,patient,samp)             
             runTask(command,fSimulate,'GATK FilterMutectCalls',logfile)


       #=== filter mutations
       if fMutCaller=='mt2' or fMutCaller=='mt':
          if (not os.path.exists(patient+'.'+mcName+'.filter.pass.vcf')):
             runTask('egrep "^#|PASS" %s.%s.filter.vcf > %s.%s.filter.pass.vcf' %(patient,mcName,patient,mcName), fSimulate, '',logfile)

					
    
    
    
       
       #=== Remove the original BAM file keeping only a new one that contains reads overlapping the variants detected using Mutect2
       if fSmallbam and fMutCaller=="mt2":
          #check if we have all of the input files
          cRemove = True
          for samp in lSamples: 
             if samp != fCtrlSample:
                if not bamCorrect(samtools,patient+'_'+samp+'.dedup.recal.mutect2.bam'):
                   if not bamCorrect(samtools,patient+'_'+samp+'.dedup.recal.bam'):
                      cRemove = False
                      printError('Not all bam files found - removal of BAM files halted') 
                   if lineNumber(patient+'.Mutect2.filter.pass.vcf') == 0:
                      cRemove = False
                      printError('Mutect2 VCF file is empty - removal of BAM files halted') 

          #make the new small BAM files
          if cRemove:        
             if not os.path.isfile(patient+'.Mutect2.filter.pass.bed'):
                command = "%s query -f '%%CHROM\t%%POS0\t%%END\n' %s.Mutect2.filter.pass.vcf > %s.Mutect2.filter.pass.bed" %(bcftools,patient,patient)
                runTask(command,fSimulate,'bcftools: convert Mutect2 VCF to BED',logfile)               
             for samp in lSamples: 
               if samp != fCtrlSample:               
                  if not bamCorrect(samtools,patient+'_'+samp+'.dedup.recal.mutect2.bam'):               
                     command = "%s intersect -abam %s_%s.dedup.recal.bam -b %s.Mutect2.filter.pass.bed > %s_%s.dedup.recal.mutect2.bam" %(bedtools,patient,samp,patient,patient,samp)
                     runTask(command,fSimulate,'bedtools: Extract variants containing reads from BAM file',logfile)
                     #index the new BAM file
                     command = "%s index %s_%s.dedup.recal.mutect2.bam -@ %s" %(samtools,patient,samp,proc)
                     runTask(command,fSimulate,'Samtools BAM indexing',logfile)
          
          #proceed to the BAM removal
          if cRemove:               
             for samp in lSamples:
                if bamCorrect(samtools,patient+'_'+samp+'.dedup.recal.mutect2.bam') or samp == fCtrlSample:                    
                   if bamCorrect(samtools,patient+'_'+samp+'.dedup.recal.bam'):
                      #remove BAM files                
                      #runTask('rm '+patient+'_'+samp+'.dedup.recal.bam',fSimulate,'',logfile)
                      #runTask('rm '+patient+'_'+samp+'.dedup.recal.bam.bai',fSimulate,'',logfile)
                      True
                else:
                    printError('Small BAM was not created - removal of BAM files halted') 
                
                    

       #report run time
       total_time = runTime(tstart)
       print('\033[0;36mTotal running time: '+total_time+'\033[0;37m')
       

       
       
def runTask(command,simulate,desc,logfile):
  if simulate:
     print('\033[0;35m'+command+'\033[0;37m') 
  else:
     ststart = time.time()    
     print('\033[0;35mRUNNING: '+command+'\033[0;37m')
     os.system(command)
     os.system('echo "%s" >> %s' %('\033[0;35m'+time.strftime('%Y-%m-%d %H:%M')+' COMPLETED: '+command+'\033[0;37m',logfile))
     if desc!='':
        text = '\033[0;36mRUNTIME: '+desc+' - '+runTime(ststart)+'\033[0;37m'
        print(text)
        os.system('echo "   %s" >> %s' % (text,logfile))


def runTime(tstart):    
   tend = time.time()
   calc_time = tend - tstart;
   hrs = math.floor(calc_time/60/60)
   mins = math.floor(calc_time/60) - hrs*60
   secs = round(calc_time - hrs*60*60 - mins*60)

   if calc_time/60>60:
      rtime =  '%.d hr %.d min %.d sec' % (hrs,mins,secs)
   elif calc_time>60:
      rtime =  '%.d min %.d sec' % (mins,secs)
   elif calc_time>=1:
      rtime =  '%s sec' % (secs)
   else:
      rtime =  '<1 sec'
   return rtime 
       

def bamCorrect(samtools,bamfile):
   command = '%s quickcheck %s' %(samtools,bamfile)
   return subprocess.getstatusoutput(command)[0]==0
     
  
def lineNumber(filename):
   if not os.path.isfile(filename):
      lines=0
   else:
      command = 'grep -v "^#" '+filename+' | wc -l'
      lines = int(subprocess.check_output(command, shell=True, preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL)))
   return lines


def printError(msg):
   print('\033[0;31mERROR: '+msg+'\033[0;37m')
   sys.exit()


def printWarning(msg):
   print('\033[0;33mWARNING: '+msg+'\033[0;37m')



if __name__ == "__main__":
    main(sys.argv[1:])
