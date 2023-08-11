update_snp_tables <- function(pung_path){
  
  # Update reference snp tables from PING/Resources/genotype_files/SNP_files 
  # if PING is more up to date than PUNG, copy from - allele sequences are updated everytime PING is run
  ping_master_snp_path <- grep("PING/Resources/genotype_resources/SNP_files$", list.dirs(path="~"), value = TRUE)
  ping_new_snp_path <- paste0(pung_path, "/Resources/SNP_files")
  if(exists(ping_master_snp_path)){
    list_files <- list.files(ping_master_snp_path, full.names = TRUE)
    date_ping_files <- file.mtime(list_files)[1]
    date_ping_new_files <- file.mtime(ping_new_snp_path)[1]
    if(date_ping_files > date_ping_new_files){
      file.copy(from = list.files(ping_master_snp_path, full.names = TRUE),
                to = ping_new_snp_path, overwrite = TRUE, recursive = FALSE,
                copy.mode = FALSE, copy.date = TRUE)
    }
  }
}