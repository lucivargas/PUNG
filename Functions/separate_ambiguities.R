separate_ambiguities <- function(unresolved_data, i){
  genotype <- unlist(strsplit(unresolved_data[i,"genotype_string"], "[+]"))
  allele1 <- unlist(strsplit(genotype[1], "\\s+"))
  allele2 <- unlist(strsplit(genotype[2], "\\s+"))
  if(length(genotype)>1){
    add <- c()
    for(each_allele1 in allele1){
      add <- c(add, paste(each_allele1, allele2, sep = "+"))
    }
  } else{
    add <- genotype
  }
  for(j in 1:length(add)){
    index <- paste0("genotype_string",j)
    unresolved_data[i,index] <- add[j]
  }
  return(unresolved_data)
}