read_PING_files <- function(data_directory, kir_filter="KIR3DL2"){
  
  setwd(data_directory)
  
  # Get unresolved genotypes from finalAlleleCalls.csv
  final_call <- read.csv("finalAlleleCalls.csv", colClasses = "character")
  names(final_call)[1] <- "sample_id"
  unresolved_data <- data.frame(ind=character(0), locus=character(0), genotype=character(0), stringsAsFactors = FALSE)
  for(col in names(final_call)[-1]){
    for(i in 1:nrow(final_call)){
      if(grepl("unresolved",final_call[i,col])){
        sample_name <- final_call[i,1]
        kir_genotype <- final_call[i,col]
        add<- data.frame(ind=sample_name, locus=col, genotype=kir_genotype, stringsAsFactors=FALSE)
        unresolved_data <- rbind(unresolved_data, add)
      }
    }
  }
  
  ## Only process KIR3DL2 - TEMPORARY
  if("all" %in% kir_filter){
    kir_filter <- names(final_call)[-1]
  }
  unresolved_data <- subset(unresolved_data, locus %in% kir_filter)
  if(nrow(unresolved_data) == 0){
    cat("No unresolved data found in the genes selected.", sep="\n")
    return(unresolved_data)
  }
  row.names(unresolved_data) <- 1:nrow(unresolved_data)
  
  # Add copy number info
  cn_df <- read.csv(file= "manualCopyNumberFrame.csv", row.names = 1)
  unresolved_data$copy_number <- NA 
  for(i in 1:nrow(unresolved_data)){
    sample_name <- unresolved_data[i, "ind"] 
    gene_name <- unresolved_data[i, "locus"] 
    cn_add <- cn_df[sample_name,gene_name]
    unresolved_data[i, "copy_number"] <- cn_add
  }
  
  
  cat(paste0("Found ", nrow(unresolved_data), " unresolved genotypes in the data provided..."), sep="\n")
  if(nrow(unresolved_data) == 0){
    return(unresolved_data)
  }
  
  # Get Detailed genotype description from iterAlleleCalls.csv
  unresolved_data$genotype_string <- NA 
  genotype_description <- read.csv("iterAlleleCalls.csv", colClasses = "character", row.names = NULL)
  for(i in 1:nrow(unresolved_data)){
    sample_name <- unresolved_data[i, "ind"]
    locus_name <- unresolved_data[i, "genotype"]
    locus_name <- unresolved_data[i, "locus"]
    if(locus_name == "KIR2DS35"){
      locus_name <- c("KIR2DS3","KIR2DS5")
    } else if(locus_name == "KIR2DL23"){
      locus_name <- c("KIR2DL2","KIR2DL3")
    } else if(locus_name == "KIR3DL1S1"){
      locus_name <- c("KIR3DL1","KIR3DS1")
    }
    index <- which(genotype_description[,1] == sample_name)
    genotype_string <- genotype_description[index, locus_name]
    # remove all 'new' snps with a *, it only indicates the allele in question was unsequenced
    genotype_string <- gsub("[EI][0-9]{1,2}_[0-9]{1,4}[.][*][\\^]{0,1}","",genotype_string)
    genotype_string <- gsub("[$]\\s|[$][+]|[$]$|[\\^]$","",genotype_string)
    # join genotype strings of alleles coding different proteins for the same gene (2DS35, 2DL23, 3DL1S1)
    if(length(locus_name)>1){
      genotype_string <- paste(genotype_string, collapse=" ") # FIX THIS, KNOW HOW TO REPORT AMBIGUOUS ALLELES
    }
    unresolved_data[i, "genotype_string"] <-genotype_string
  }
  #write.csv(unresolved_data, file=paste0(data_directory,"unresolved_data.csv"), row.names = F, quote = F)
  
  # # Separate each genotype possibility found into a new column in the unresolved_data
  # max_ambiguities <- c()
  # for(i in 1:nrow(unresolved_data)){
  #   genotype <- unlist(strsplit(unresolved_data[i,"genotype_string"], "[+]"))
  #   allele1 <- unlist(strsplit(genotype[1], "\\s+"))
  #   allele2 <- unlist(strsplit(genotype[2], "\\s+"))
  #   add <- length(allele1) * length(allele2)
  #   max_ambiguities <- c(max_ambiguities, add)
  # }
  # max_ambiguities <- max(max_ambiguities)
  # # add new columns
  # new_columns <- paste0("genotype_string", seq(1,max_ambiguities))
  # for(new_column in new_columns){
  #   unresolved_data[[new_column]] <- NA
  # }
  # #separate ambiguous genotypes
  # for(i in 1:nrow(unresolved_data)){
  #   unresolved_data <- separate_ambiguities(unresolved_data, i)
  # }
  # 
  return(unresolved_data)
}