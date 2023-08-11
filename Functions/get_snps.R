get_snps <- function(unresolved_data, i){
  # Get genotype
  genotype <- unresolved_data[i,"genotype_string"]
  if(is.na(genotype)){next}
  genotype <- unlist(strsplit(unresolved_data[i,"genotype_string"],  "[+]|\\s+"))
  # Get SNPs
  snp_list <- unlist(lapply(genotype, function(x) unlist(strsplit(x, "\\$"))[2])) #function(x) x ^ 2
  snp_list <- unlist(strsplit(snp_list,"\\^"))
  snp_list <- unique(snp_list)
  snp_list <- as.character(na.omit(snp_list))
  return(snp_list)
}