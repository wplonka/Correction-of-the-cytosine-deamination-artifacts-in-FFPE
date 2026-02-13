1. Run microSEC_mutation_info_file_preparation.R
2. Run microSEC_sample_info_file_preparation.R (If you are working with chromosomes separetly edit this code)
3. Make a list of info sample files and save eg. info_files_list.txt
4. Run microSEC_run.sh
   
In case of errors related to incorrect data types in the chromStart and chromEnd columns (e.g., coercion errors, non-numeric values, or issues converting to integer), the following lines should be added to the original microSEC code to explicitly enforce numeric-to-integer conversion:

download_region$chromStart <- as.integer(as.numeric(download_region$chromStart))
download_region$chromEnd   <- as.integer(as.numeric(download_region$chromEnd))

5. (Optional) merge results to one tsv if you worked with chromosomes separetly.
6. Make a list of vcf and output tsv from MicroSEC files eg. MicroSEC_paired_files.txt
7. Run microSEC_filtration.sh

