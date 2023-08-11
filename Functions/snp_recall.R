snp_recall <- function(unresolved_data, index=i, each_snp, snp_output1){
  
  sample_name <- unresolved_data[index,"ind"]
  gene_name <- unresolved_data[index, "locus"]
  position_snp <- unlist(strsplit(each_snp, "\\."))[1]
  
  # check if there is enough depth in that position
  total <- sum(as.numeric(snp_output1))
  if(total < 20){
    cat(paste0("\tNot enough reads for ", position_snp, " in ", sample_name," (",total,"). Removing SNP..."), sep="\n")
    unresolved_data <- remove_snp(unresolved_data, index, sample_name, each_snp)
    return(unresolved_data)
  }
  # check depth per each base/indel to make call
  ratios <- as.numeric(snp_output1)/total
  # if below 20%, discard
  ratios[which(ratios < 0.20)] <- 0
  # call SNP
  snp_call <- names(snp_output1) [ ratios != 0]
  
  # Add new call to snp_final
  null_to_add <- 3-length(snp_call)
  null_to_add <- rep("N",null_to_add)
  snp_call_final <- c(snp_call, null_to_add)
  # Import file with the Typing of SNPs in a given sample 
  snp_file_type <- paste0(results_directory,"/snp_output/", "final_", gene_name,"_",sample_name,"_SNP.csv")
  snp_final <- read.csv(snp_file_type, stringsAsFactors = F, row.names = 1, check.names = F)
  # Mutate snp_final
  snp_final[, position_snp] <- snp_call_final
  # Save updated typing table
  write.csv(snp_final, file=snp_file_type, quote = F, row.names = T)
  
  # if it is indeed an indel, pass QC and return
  if(any(snp_call %in% c(".","INS"))){
    cat(paste0("\t",each_snp," passed QC."), sep="\n")
    return(unresolved_data)
  }
  
  # If bases are called
  if(length(snp_call) ==1){
    
    # check if the new snp call indeed different from the sequence of the allele its in
    # import consensus
    gene_snp_path <- list.files(paste0(pung_path,"/Resources/SNP_files"), pattern = paste0("^",gene_name), full.names = T)
    snp_df <- read.csv(gene_snp_path, colClasses = "character", row.names = 1)
    exon_indexes <- grep("^X{0,1}E", names(snp_df))
    snp_df <- snp_df[,exon_indexes]
    genotype <- unlist(strsplit(unresolved_data[index,"genotype_string"],  "[+]|\\s+"))
    for(each_allele in genotype){
      if(grepl(gsub("\\.","\\\\.", each_snp),each_allele)){
        base_allele <- unlist(lapply(each_allele, function(x) unlist(strsplit(x, "\\$"))[1])) #function(x) x ^ 2
        position_snp <- unlist(strsplit(each_snp, "\\.{1}"))[1]
        nt_consensus <- snp_df[base_allele,position_snp]
        # If the base allele is unsequenced, lower resolution and get closest allele
        is_unsequenced <- is.na(nt_consensus)
        if(is_unsequenced){
          resolution <- nchar(unlist(strsplit(base_allele,"[*]"))[2])
          if(resolution == 7){
            base_allele <- str_extract(base_allele, "KIR.*[*][0-9]{5}")
          } else if(resolution == 5){
            base_allele <- str_extract(base_allele, "KIR.*[*][0-9]{3}")
          } 
          base_allele <- grep(gsub("[*]","[*]",base_allele), row.names(snp_df), value=T)
          nt_consensus <- snp_df[base_allele,position_snp]
          nt_consensus <- unique(nt_consensus)
          nt_consensus <- paste(nt_consensus, collapse = "")
        }
        if(snp_call == nt_consensus){
          old_call <- gsub("[*]","[*]",
                           gsub("[$]","[$]",
                                gsub("\\.","\\\\.", 
                                     gsub("[\\^]","[\\\\^]",each_allele))))
          new_call <- gsub(gsub("\\.","\\\\.", each_snp),"",each_allele)
          new_call <- gsub("[$][\\^]","$",new_call)
          new_call <- gsub("[\\^]{2}","^",new_call)
          new_call <- gsub("[\\^]$","",new_call)
          new_call <- gsub("[$]$","",new_call)
          # Change in the main data
          new_string <- unresolved_data[index,"genotype_string"] 
          new_string <- gsub(old_call, new_call, new_string)
          unresolved_data[index,"genotype_string"]  <- new_string
          next
        } else{
          # First mutate new snp
          each_snp1 <- paste(position_snp, snp_call, sep = ".")
          old_call <- gsub("[*]","[*]", gsub("[$]","[$]", gsub("\\.","\\\\.", each_allele)))
          new_call <- gsub(gsub("\\.","\\\\.", each_snp),each_snp1,each_allele)
          # Change in the main data
          new_string <- unresolved_data[index,"genotype_string"] 
          new_string <- gsub(old_call, new_call, new_string)
          unresolved_data[index,"genotype_string"]  <- new_string
          next
        }
      }
    }
    cat(paste0("\t",each_snp," passed QC."), sep="\n")
    
  } else if(length(snp_call) ==2){
    
    # check if the new snp call indeed different from the sequence of the allele its in
    # import consensus
    gene_snp_path <- list.files(paste0(pung_path,"/Resources/SNP_files"), pattern = paste0("^",gene_name), full.names = T)
    snp_df <- read.csv(gene_snp_path, colClasses = "character", row.names = 1)
    exon_indexes <- grep("^X{0,1}E", names(snp_df))
    snp_df <- snp_df[,exon_indexes]
    genotype <- unlist(strsplit(unresolved_data[index,"genotype_string"],  "[+]|\\s+"))
    # check if any base on the new call is 'new'
    base_alleles <- unlist(lapply(genotype, function(x) unlist(strsplit(x, "\\$"))[1])) #function(x) x ^ 2
    position_snp <- unlist(strsplit(each_snp, "\\.{1}"))[1]
    nt_consensus <- snp_df[base_alleles,position_snp]
    # add new SNP call into string
    if(any(snp_call %!in% nt_consensus)){
      snp_call <- snp_call[snp_call %!in% nt_consensus]
      snp_call <- unique(snp_call)
      gsub("\\.\\.", paste0(".",snp_call))
      
      #get old call and new call
      old_call <- gsub("[*]","[*]",
                       gsub("[$]","[$]",
                            gsub("\\.","\\\\.", 
                                 gsub("[\\^]","[\\\\^]",each_allele))))
      new_call <- gsub(gsub("\\.","\\\\.", each_snp),
                       gsub("\\.\\.", paste0(".",snp_call), each_snp),
                       each_allele)
      new_call <- gsub("[$]$","",new_call)
      # Change in the main data
      new_string <- unresolved_data[index,"genotype_string"] 
      new_string <- gsub(old_call, new_call, new_string)
      unresolved_data[index,"genotype_string"]  <- new_string
      next
    }
    
    cat(paste0("\t",each_snp," passed QC."), sep="\n")
    
  } else{
    cat(paste0("\t",each_snp," failed QC."), sep="\n")
    unresolved_data <- remove_snp(unresolved_data, sample_name, each_snp)
    next
  }
  return(unresolved_data)
}