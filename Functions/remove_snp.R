# Remove snp from unresolved_data
remove_snp <- function(unresolved_data, index, sample_name, each_snp){
  
  # get gene name
  gene_name <- unresolved_data[index, "locus"]
  
  # get original string
  old_string <- unresolved_data[index, "genotype_string"]
  
  # get position
  position_snp <- unlist(strsplit(each_snp, "\\."))[1]
  
  # remove snp from string
  match_string <- gsub("\\.","[.]",each_snp)
  new_string <- gsub(paste0("[\\^]{0,1}",match_string),"", old_string)
  new_string <- gsub("[$]\\s+"," ",new_string)
  new_string <- gsub("[$][+]","+",new_string)
  new_string <- gsub("[$]$","",new_string)
  new_string <- gsub("[$][\\^]","$",new_string)
  unresolved_data[index, "genotype_string"] <-  new_string
  
  # if any of ambiguities the are the same, join
  genotype1 <- unlist(strsplit(unresolved_data[index,"genotype_string"], "[+]"))
  list_genotype <- c()
  for(allele in genotype1){
    add <- unlist(strsplit(allele, "\\s+"))
    add <- unique(add)
    add <- paste(add, collapse=" ")
    list_genotype <- c(list_genotype, add)
  }
  unresolved_data[index,"genotype_string"] <- paste(list_genotype, collapse="+")
  
  # Mutate the snp_output table 
  # Import file with the Typing of SNPs in a given sample 
  snp_file_type <- paste0(results_directory,"/snp_output/", "final_", gene_name,"_",sample_name,"_SNP.csv")
  snp_final <- read.csv(snp_file_type, stringsAsFactors = F, row.names = 1)
  names(snp_final) <- gsub("^X","",names(snp_final))
  # Generate new call
  nt_to_remove <- unlist(strsplit(each_snp, "\\."))[2]
  new_call <- snp_final[, position_snp]
  new_call <- gsub(nt_to_remove, "N", new_call)
  # Mutate snp_final
  if(!all(unique(new_call) == "N")){
    snp_final[, position_snp] <- new_call
  }
  # Save updated typing table
  write.csv(snp_final, file=snp_file_type, quote = F, row.names = T)
  
  return(unresolved_data)
}
