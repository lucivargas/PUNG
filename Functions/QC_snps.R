source(paste0(pung_path,"/Functions/snp_recall.R"))

QC_snps <- function(unresolved_data, index=i, snp_list){
  
  #### QC - Check quality 
  # Prep
  sample_name <- unresolved_data[index,"ind"]
  gene_name <- unresolved_data[index, "locus"]
  
  # import file with Depth of SNP in the given sample
  snp_file_depth <- paste0(data_directory,"alignmentFiles/",
                     sample_name,"/iterAlign/final_", gene_name,"_",sample_name,"_DP.csv")
  snp_output <- read.csv(snp_file_depth, stringsAsFactors = F)
  snp_output <- data.frame(t(snp_output), stringsAsFactors = F)
  names(snp_output) <- snp_output[1,]
  snp_output <- snp_output[-1,]
  row.names(snp_output) <- gsub("^X","",row.names(snp_output))
  
  # Iterate 'new' snps found
  for(each_snp in snp_list){
    
    # Define SNP variables
    position_snp <- unlist(strsplit(each_snp, "\\."))[1]
    new_nt <- unlist(strsplit(each_snp, "\\."))[2]
    snp_index <- which(row.names(snp_output) == position_snp)
    snp_output1 <- snp_output[snp_index,]
    
    # Some snp calls are indels or are not made, indicated by (.), try to solve it.
    if(grepl("\\.\\.", each_snp)){
      unresolved_data <- snp_recall(unresolved_data, index=index, each_snp, snp_output1)
      next
    }
    
    # Get depth of reads
    depth_all <- sum(as.numeric(snp_output1[]))
    depth_nt <- as.numeric(snp_output1[,new_nt])
    depth_ratio <- depth_nt/depth_all
    
    # If less than 20 reads, reject it, or less than 25% of total reads
    if(depth_ratio <0.25 | depth_nt<20){
      cat(paste0("\tNot enough reads for ", each_snp, " in ", sample_name," (",depth_nt,"). Removing SNP..."), sep="\n")
      unresolved_data <- remove_snp(unresolved_data, index, sample_name, each_snp)
    }else{
      cat(paste0("\t",each_snp," passed QC."), sep="\n")
    }
  }
  
  return(unresolved_data)
    
}