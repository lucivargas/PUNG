kir_recall <- function(unresolved_data, index=i, allele){
  
  ## 1.Determining the sequences of potential new KIR alleles from unresolved genotypes
  
  # Import sequences of alleles deposited on IPD-KIR
  gene_name <- unresolved_data[index,"locus"]
  gene_snp_path <- list.files(paste0(pung_path,"/Resources/SNP_files"), pattern = paste0("^",gene_name), full.names = T)
  snp_df <- read.csv(gene_snp_path, colClasses = "character", row.names = 1)
  exon_indexes <- grep("^X{0,1}E", names(snp_df))
  snp_df <- snp_df[,exon_indexes]
  
  # Filter sequence of the most similar allele
  base_allele <- gsub("\\$.*$","",allele)
  allele_index <- grep(gsub("[*]","[*]",base_allele), row.names(snp_df))[1]
  if(is.na(allele_index)){
    cat(red(paste0("\tWarning: Allele ",base_allele," not found in IPD-KIR reference tables. Could not verify ",allele,".")), sep="\n")
    return(unresolved_data)
  }
  new_sequence <- snp_df[allele_index,]
  #if many alleles are found, keep first one
  if(nrow(new_sequence)>1){
    new_sequence <- new_sequence[1,]
  }
  
  # Create list of 'new' SNPs 
  snp_list <- unlist(strsplit(allele, "\\$"))[-1] 
  snp_list <- unlist(strsplit(snp_list,"\\^"))
  # remove all 'new' snps with a *, it only indicates the position in question was unsequenced
  snp_list <- snp_list[grep("[*]",snp_list, invert=T)]
  
  # Mutate in the consensus sequence the 'new' SNPs observed
  for(each_snp in snp_list){
    position_snp <- unlist(strsplit(each_snp, "\\."))[1]
    new_nt <- unlist(strsplit(each_snp, "\\."))[2]
    if(new_sequence[,position_snp] == new_nt){
      cat(paste0("\tThe new SNP ", each_snp, " is already found in the consensus sequence of ", base_allele), sep="\n")
    }
    new_sequence[,position_snp] <- new_nt
  }
  
  
  ## 2. Given that IPD-KIR is a constantly updated public sourced database, 
  ## we will check if the potential new allele found has been recently deposited on IPD-KIR. 
  ## If the potential new allele sequence matches an existing allele, we will return the 
  ## updated genotype and end the validation.
  
  # Check if this sequence matches existing alleles deposited on IPD-KIR
  query <- paste(new_sequence, collapse = "")
  query <- gsub("[*]",".", query)
  match_list <- c()
  for(each_allele in row.names(snp_df)){
    target <- paste(snp_df[each_allele,], collapse = "")
    target <- gsub("[*]",".", target)
    is_match <- grepl(query, target) | grepl(target, query) 
    #which( new_sequence != snp_df[each_allele,]) # view differences
    if(is_match){
      match_list <- c(match_list, each_allele)
    }
  }
  
  
  ## Format and report results
  # If "new" allele found has been deposited on IPD-KIR, end validation, 
  if(!isempty(match_list)){
    cat(paste0("\tThe new allele ",allele," has already been deposited on IPD-KIR! Updating results table..."), sep="\n")
    # Trim down list if there are many high resolution alleles found
    if(length(match_list) > 1){
      for(k in c(7,5,3)){
        match_list <- unique(str_extract(match_list, paste0("^KIR.*[*][0-9]{3,",k,"}"))) #fields digits 3-2-2
        if(length(match_list)==1){break}
      }
    }
    # if there is ambiguity, format string 
    if(length(match_list)>1){
      match_list <- paste(match_list, collapse = " ")
    }
    
    # Change allele AND all the ambiguous alleles in genotyping list
    old_allele <- gsub("[*]","[*]",allele)
    old_allele <- gsub("[$]","[$]",old_allele)
    old_allele <- gsub("\\.","\\\\.",old_allele)
    
    genotype1 <- unlist(strsplit(unresolved_data[index,"genotype_string"],  "[+]"))
    index_genotype <- grep(old_allele, genotype1)
    genotype1[index_genotype] <- match_list
    unresolved_data[index,"genotype_string"] <- paste(genotype1, collapse="+")
    # 
    # break
    # # NEW GENOTYPING
    # #import sample name
    # seq_path <- list.files(paste0(data_directory,"alignmentFiles/"), recursive=T, pattern=paste("final",gene_name,sample_name,"SNP", sep = ".*"), full.names = T)
    # seq <- read.csv(seq_path, colClasses = "character")
    # #keep exons
    # seq <- seq[,grep("^E", names(seq))]
    # # Get sequence from the other allele
    # seq_known <- snp_df[match_list,] 
    # seq_other <- c()
    # for(var in 1:length(seq)){
    #   geno <- seq[c(1:2),var]
    #   if(geno[2] == "N"){
    #     geno[2] <- geno[1]
    #   }
    #   if(all(geno[2]=="N")){
    #     geno <- c("*","*") 
    #   }
    #   query <- seq_known[,var]
    #   delete_var <- which(geno == query)[1]
    #   if(isempty(delete)){
    #     cat(paste0("Could not find ",query," in ",names(seq)[var]), sep="\n")
    #   }
    #   geno <- geno[-delete_var]
    #   seq_other <- c(seq_other,geno)
    # }
    # seq_other <- as.data.frame(rbind(names(seq_known), seq_other), stringsAsFactor = FALSE)
    # names(seq_other) <- seq_other[1,]
    # seq_other <- seq_other[-1,]
    # 
    # # Format and Save
    # # add new genotype to unresolved data
    # genotype1 <- gsub("[*]","[*]",allele)
    # genotype1 <- gsub("[$]","[$]",genotype1)      
    # genotype1 <- gsub("\\.","[.]",genotype1) 
    # genotype1 <- gsub("\\^","\\\\^",genotype1) 
    # confirmed_snps <- unlist(strsplit(allele,"[$]"))[2]
    # confirmed_snps <- unlist(strsplit(confirmed_snps, "\\^"))
    # 
    # # change the old annotaded allele string for the new allele found on IPD
    # # unresolved_data[index,"genotype_string"] <- gsub(genotype1, match_list, unresolved_data[i,"genotype_string"])
    # old_string <- unlist(strsplit(unresolved_data[index,"genotype_string"] , "[+]"))
    # amb_index <- grep(genotype1, old_string)
    # old_string[amb_index] <- match_list
    # unresolved_data[index,"genotype_string"] <- paste(old_string, collapse="+")
    # unresolved_data <- separate_ambiguities(unresolved_data, index)
    # 
    # # Remove found snps in our data
    # sample_name <- unresolved_data[index,"ind"]
    # copy_number <- unresolved_data[index, "copy_number"]
    # if(copy_number >1){
    #   for (confirmed_snp in confirmed_snps){
    #     position_snp <- unlist(strsplit(confirmed_snp,"\\."))[1]
    #     base_snp <- unlist(strsplit(confirmed_snp,"\\."))[2]
    #     #import sample sequence
    #     seq_path <- list.files(paste0(data_directory,"alignmentFiles/"), recursive=T, pattern=paste("final",gene_name,sample_name,"SNP", sep = ".*"), full.names = T)
    #     seq <- read.csv(seq_path, colClasses = "character")
    #     #keep exons
    #     seq <- seq[,grep("^E", names(seq))]
    #     
    #     
    #     geno <- seq[,position_snp]
    #     geno <- geno[geno != "N"]
    #     if(!all(geno==position_snp)){
    #       unresolved_data <- remove_snp(unresolved_data, sample_name, confirmed_snp)
    #     }
    #   }
    # }
    # 
    # # -----END VALIDATION ----
    # next
  }
  
  return(unresolved_data)
  
}
