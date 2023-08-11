suppressPackageStartupMessages(library(operators)) #%!in%
suppressPackageStartupMessages(library(crayon)) #color in cat

check_files <- function(data_directory, results_directory, pung_path){
  ## Check for finalAlleleCalls.csv, iterAlleleCalls.csv, Variant Calling,  Filtered FastQ files
  setwd(data_directory)
  files <- c("finalAlleleCalls.csv", "iterAlleleCalls.csv", "alignmentFiles/", "extractedFastq/")
  
  for(each_file in files){
    if(!file.exists(each_file)){
      cat(red(paste0("Your data directory is missing ",each_file," output from PING. Provide it before running PUNG.")), sep="\n")
    }
    ## Check for Variant Calling data - per sample
    if(each_file == "alignmentFiles/"){
      if(file.exists("finalAlleleCalls.csv")){
        # Get list of samples
        sample_names <- read.csv("finalAlleleCalls.csv", colClasses = "character")[,1]
        # Get list of snp output files
        snp_files <- list.files(paste0(data_directory,"alignmentFiles/"), recursive = T, pattern = "final.*_DP")
        snp_files <- unlist(unique(gsub("_KIR[0-9A-Z]*_","_",snp_files)))
        # Check if they match
        if(length(sample_names) < length(snp_files)){
          missing_samples <- sample_names [ sample_names %!in% gsub("final_|_DP\\.csv","",basename(snp_files)) ]
          cat(red("Some Variant Calling data for your samples are missing from 'alignmentFiles/sample/iterAlign/'. Provide it before running PUNG."), sep="\n")
          cat(silver(paste(c("\tMissing: ", missing_samples))), sep="\n\t" )
        }
      }
    }
    
  }
  
  # Check reading permission on data folder, (0 for success and -1 for failure)
  read_perm <- ifelse(file.access(data_directory, mode =4) == 0,T,F)[[1]]
  if(!read_perm){
    cat(red(paste0("You do not have permission to read files in ",data_directory,". Provide it before running PUNG.")), sep="\n")
    
  }
  
  # Check writing permission on results folder, (0 for success and -1 for failure)
  write_perm <- ifelse(file.access(dirname(results_directory), mode =2) == 0,T,F)[[1]]
  if(write_perm == FALSE){
    cat(red(paste0("You do not have permission to create a folder in ", dirname(results_directory))), sep="\n")
  }

  # Create an output folder. 
  if(!dir.exists(results_directory)){
    dir.create(results_directory)
  }
  
  cat("Files ok.", sep="\n")
  
}

  
  