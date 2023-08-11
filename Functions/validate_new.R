suppressPackageStartupMessages(library(tibble)) #add_colum()
suppressPackageStartupMessages(library(pracma)) # isempty()
suppressPackageStartupMessages(library(stringr)) #str_extract()
suppressPackageStartupMessages(library(Biostrings)) # DNAString, reverseComplement
suppressPackageStartupMessages(library(stats)) # na.omit()
source(paste0(pung_path,"/Functions/remove_snp.R"))
source(paste0(pung_path,"/Functions/separate_ambiguities.R"))
source(paste0(pung_path,"/Functions/read_PING_files.R"))
source(paste0(pung_path,"/Functions/update_snp_tables.R"))
source(paste0(pung_path,"/Functions/QC_snps.R"))
source(paste0(pung_path,"/Functions/kir_blastn.R"))
source(paste0(pung_path,"/Functions/kir_recall.R"))
source(paste0(pung_path,"/Functions/get_snps.R"))

validate_new <- function(data_directory, results_directory, pung_path, kir_filter="KIR3DL2"){
  
  # Record log of run
  log_file <- paste0(results_directory, "PUNG_log.txt")
  if(file.exists(log_file)){file.remove(log_file)}
  sink(paste0(results_directory, "PUNG_log.txt"), append=FALSE, split=TRUE)
  start.time <- Sys.time()
  
  ### PREPARE DATA
  
  # Read PING files and make our own dataframe
  if(!grepl("/$",data_directory)){
    data_directory<-paste0(data_directory,"/")
  }
  cat(paste0("Running PUNG on ",data_directory,"\n"), sep="\n")
  unresolved_data <- read_PING_files(data_directory, kir_filter=kir_filter)
  if(nrow(unresolved_data) == 0){
    return()
  }
  
  # Update reference snp tables from PING/Resources/genotype_files/SNP_files 
  update_snp_tables(pung_path)
  
  
  ### VALIDATE NEW ALLELE
  ### Get sample sequence, QC, check if its been recently deposited on IPD, if not then filter FASTQ files and BLAST reads
  
  # Iterate samples
  for(i in 1:nrow(unresolved_data)){
    
    # Get gene name
    gene_name <- unresolved_data[i,"locus"]
    
    # Get sample name
    sample_name <- unresolved_data[i,"ind"]
    
    ## Get list of SNPs and validate each one 
    snp_list <- get_snps(unresolved_data, i)
    
    # Copy snp_output files to results folder
    snp_new_folder <- paste0(results_directory, "/snp_output/")
    if(!dir.exists(snp_new_folder)){
      dir.create(snp_new_folder)
    }
    snp_file <- paste0(data_directory,"alignmentFiles/",
                        sample_name,"/iterAlign/final_", gene_name,"_",sample_name,"_SNP.csv")
    snp_new_file <- paste0(snp_new_folder, basename(snp_file))
    file.copy(to=snp_new_file, from=snp_file, overwrite = TRUE, copy.mode = FALSE, copy.date = FALSE)
    
    
    ## 1. Running Quality Control (QC) of new SNPs found. 
    cat(paste0("(",i,"/",nrow(unresolved_data),") Checking sequence of ", unresolved_data[i,"genotype_string"], "..."),sep="\n")
    unresolved_data <- QC_snps(unresolved_data, index=i, snp_list)
    # End validation if all snps have been removed
    snp_list <- get_snps(unresolved_data, i)
    if(isempty(snp_list)){next}
    
    
    ## 2. Create bioinformatic probes for each of the novel single nucleotide 
    ## polymorphisms (SNPs) found in the potential new alleles.
    ## 3. The probes created will be used to filter relevant sequencing reads in the FASTQ files.
    ## 4. Because KIR genes have high sequence identity, it will be necessary to confirm that the
    ## reads containing novel SNPs are specific to the investigated gene. For this purpose, we
    ## will use the sequencing reads to search the human dataset with BLAST. If the reads are 
    ## specific to the gene being investigated, we will move forward with the validation.
    unresolved_data <- kir_blastn(unresolved_data, index=i, snp_list, pung_path)
    # End validation if all snps have been removed
    snp_list <- get_snps(unresolved_data, i)
    if(isempty(snp_list)){next}
    
    ## 5. Determining the sequences of potential new KIR alleles from unresolved genotypes
    ## 6. Given that IPD-KIR is a constantly updated public sourced database, 
    ## we will check if the potential new allele found has been recently deposited on IPD-KIR. 
    genotype <- unlist(strsplit(unresolved_data[i,"genotype_string"],  "[+]|\\s+"))
    for(allele in genotype){
      updated_genotype <- unlist(strsplit(unresolved_data[i,"genotype_string"],  "[+]|\\s+"))
      if(allele %in% updated_genotype){
        if(!grepl("[$]",allele)){next}
        unresolved_data <- kir_recall(unresolved_data, index=i, allele)
      }
      if(!grepl("[$]",unresolved_data[i,"genotype_string"])){break} # if new alleles have been removed
    }
    
    
    ## 7. Based on remaining SNPs that are confirmed, confirm the phase, if possible 
    # TBD
    
  } # END Iterate samples
  

  # Summarize results from unresolved_data, and make report on NewAllelesList
  unresolved_data$snp_list <- NA
  # add list of SNPs to be sequenced 
  for(i in 1:nrow(unresolved_data)){
    snp_list <- get_snps(unresolved_data, i)
    snp_list <- paste(snp_list, collapse = "^")
    snp_list <- paste(unresolved_data$locus,snp_list, sep="*")
    unresolved_data$snp_list[i] <- snp_list
  }
  write.csv(unresolved_data, file=paste0(results_directory,"newAllelesList.csv"), quote=F, row.names=F)
  
  # Save new genotypes to the main table and report it on finalAlleleCalls.csv
  final_call <- read.csv("finalAlleleCalls.csv", colClasses = "character")
  for(i in 1:nrow(unresolved_data)){
    sample_name <- unresolved_data[i, "ind"]
    gene_name <- unresolved_data[i, "locus"]
    new_genotype <- unresolved_data[i, "genotype_string"]
    final_call[which(final_call[1]== sample_name), gene_name] <- new_genotype
  }
  write.csv(final_call, file=paste0(results_directory, "finalAlleleCalls-PUNG.csv"), quote=F, row.names=F)
  
  # Close time
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste0("\n",time.taken,"\n"), sep="\n")
  
  # Final messages
  resolved_samples <- length(grep("[$]|_", unresolved_data$genotype_string, invert = T))
  cat(paste0("PUNG completely resolved ", resolved_samples," out of ", nrow(unresolved_data), " genotypes"), sep="\n")
  cat("\nEND OF RUN", sep="\n")
  
  # Close log 
  sink()

}
