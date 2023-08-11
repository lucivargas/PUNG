kir_blastn <- function(unresolved_data, index=i, snp_list, pung_path){
  
  
  ### BLAST
  
  # Iterate 'new' snps found
  for(each_snp in snp_list){
    
    # skip blast if position is a potential deletion (vcf has reads only for .)
    if(grepl("\\.\\.$", each_snp)){
      next
    }
    ### PREPARE BLAST ###
    
    # Only run blast if results do not already exist
    gene_name <- unresolved_data[index, "locus"]
    sample_name <- unresolved_data[index, "ind"]
    blast_results <- paste0(results_directory,"/filteredFastq/",gene_name, ".", each_snp,"/", sample_name, "_kir_blastn.txt")
    do_run <- !file.exists(blast_results)
    if(do_run){
    
      ## Get probes for this SNP - from genomic sequence (w/ introns) 
      # get the sample sequence for this gene
      seq_file <- paste0(data_directory,"alignmentFiles/",
                         sample_name,"/iterAlign/final_", gene_name,"_",sample_name,"_SNP.csv")
      snp_df_probe <- read.csv(seq_file, colClasses = "character")
      # Format
      row.names(snp_df_probe) <- snp_df_probe[,1]
      snp_df_probe <- snp_df_probe[,-1]
      snp_df_probe <- snp_df_probe[c(1:2),] # delete 3rd row
      # remove pseudoexon prefix P
      names(snp_df_probe) <- gsub("^P","",names(snp_df_probe))
      
      # Get genotyping for the region surrounding snp of interest (where the probe will be)
      position_snp <- unlist(strsplit(each_snp, "\\."))[1]
      position_index <- which(names(snp_df_probe) %in% position_snp)
      position_index <- seq(from=position_index-20, to=position_index+20) # get 41 bp
      snp_df_probe <-  snp_df_probe[,position_index]
      
      # Convert genotype calling table to probe sequences
      for(position in 1:length(snp_df_probe)){
        allele1 <- snp_df_probe[1,position]
        allele2 <- snp_df_probe[2,position]
        if(allele2 =="N" | allele2 =="."){
          snp_df_probe[2,position] <- allele1
        }
      }
      
      # Check if there are other snps in LD with our target snp
      is_there_ld <- any(c(snp_df_probe[1,] != snp_df_probe[2,])[- which(names(snp_df_probe) == position_snp)])
      if(is_there_ld){
        # create all combinations of sequences
        ld_snps <- which(c(snp_df_probe[1,] != snp_df_probe[2,]))
        for(ld_snp in ld_snps){
          # create new sequence to shuffle positions
          add <- snp_df_probe
          # invert
          add[1,ld_snp] <- snp_df_probe[2,ld_snp]
          add[2,ld_snp] <- snp_df_probe[1,ld_snp]
          # append
          snp_df_probe <-rbind(snp_df_probe, add)
        }
        # keep unique combinations
        snp_df_probe <- unique(snp_df_probe)
      }
      
      # select mutation
      new_nt <- unlist(strsplit(each_snp, "\\."))[2]
      index_seq <- which(snp_df_probe[,position_snp]  == new_nt) 
      if(isempty(index_seq)){
        snp_df_probe[, position_snp] <- new_nt
        index_seq <- which(snp_df_probe[,position_snp]  == new_nt) 
      }
      snp_df_probe <- snp_df_probe[index_seq,] 
      
      # if homozygous for mutation, keep only unique rows
      if (nrow(snp_df_probe) >1){
        # keep unique combinations
        snp_df_probe <- unique(snp_df_probe)
      }
      
      # make forward and reverse complements
      probes <- c()
      for(p in 1:nrow(snp_df_probe)){
        probe1 <- paste(snp_df_probe[p,], collapse = "")
        probe2 <- reverseComplement(DNAString(probe1))
        probe2 <- as.character(probe2)
        probes <- c(probes, probe1, probe2)
      }
      
      
      ## Filter relevant sequencing reads in the FASTQ alignment files.
      
      # Find FASTQ files
      sample_name <- unresolved_data[index,"ind"] 
      fastq_files <- list.files(path = paste0(data_directory,"extractedFastq/"), pattern = sample_name, full.names = TRUE)
      
      # Search probes in each FASTQ file
      reads <- c()
      for (fastq_file in fastq_files){
        add <- readLines(fastq_file)
        indexes <- grep(paste(probes, collapse="|"), add)
        # get header line (-1) and quality lines (+2)
        indexes <- c(indexes-1, indexes, indexes+1, indexes+2)
        indexes <- indexes[order(indexes)]
        add <- add[indexes]
        reads <- c(reads, add)
        rm(add)
      }
      
      if(isempty(reads)){
        cat(paste0("\tNo matches for ",each_snp,
                   " in sequencing reads of sample id ", unresolved_data[index,"ind"],". Removing SNP..."), sep="\n")
        unresolved_data <- remove_snp(unresolved_data, index, sample_name, each_snp)
        next()
      } else if(length(reads)<15){
        cat(paste0("\tNot sufficient depth (",length(reads)," reads) for 'new' SNP ",each_snp,
                   " in sequencing reads for sample id ", unresolved_data[index,"ind"]), sep="\n")
        unresolved_data <- remove_snp(unresolved_data, index, sample_name, each_snp)
        next
      }
      
      
      # Create complement reverse of reads
      reads_cr <- reads
      index_seqs <- grep("^[ATCG]*$", reads)
      for(seq in index_seqs){
        convert_seq <- reads_cr[seq]
        convert_seq <- reverseComplement(DNAString(convert_seq))
        convert_seq <- as.character(convert_seq)
        reads_cr[seq] <- convert_seq
        # add reverse complement (-REV) label 
        reads_cr[seq-1] <- sub(" ","-REV ", reads_cr[seq-1])
      }
      reads <- c(reads,reads_cr)
      
      # Save to filteredFastq in PING-NEW resources
      if(!dir.exists(paste0(results_directory,"/filteredFastq/"))){
        dir.create(paste0(results_directory,"/filteredFastq/"))
      }
      filteredFastq_path <- paste0(results_directory,"/filteredFastq/", gene_name,".",each_snp,"/")
      filteredFastq_path <- gsub("[$]",".", filteredFastq_path)
      filteredFastq_path <- gsub("[*]",".", filteredFastq_path)
      if(!dir.exists(filteredFastq_path)){
        dir.create(filteredFastq_path)
      }
      filteredFastq_file <- gsub("KIR_[12]","KIR",gsub("\\.gz","",basename(fastq_file)))
      filteredFastq_file <- paste0(filteredFastq_path, filteredFastq_file)
      #filteredFastq_file <- paste0(filteredFastq_file, ".gz")
      if(file.exists(filteredFastq_file)){
        system(paste0("rm ", filteredFastq_file))
      }
      if(file.exists(paste0(filteredFastq_file,".gz"))){
        system(paste0("rm ", paste0(filteredFastq_file,".gz")))
      }
      writeLines(reads, filteredFastq_file)
      system(paste0("gzip ", filteredFastq_file))
      # save probes sequences too
      #probes_path <- gsub("\\.fastq\\.gz$","_probes.txt",filteredFastq_path)
      #write.table(c(probe1, probe2), file=probes_path, row.names = F, col.names = F)
      
      
      ## Perform BLAST of these sequences to KIR alignment file
      # Make FASTA input
      reads_fasta <- gsub("^@","> @", reads)
      # keep only header and sequences
      reads_fasta <- reads_fasta[-c(grep("F[^@]",reads_fasta), # remove quality line
                                    which(reads_fasta == "+"))] #remove separator
      #save fasta
      reads_fasta_path <- gsub("\\.fastq",".fsa",filteredFastq_file)
      writeLines(reads_fasta, reads_fasta_path)
      
      # Input files for BLAST 
      # To create db of all KIR fasta # makeblastdb -in '~/PING-NEW/Resources/AllKIRGen_noGapsAlleleNames_gen.fas' -dbtype nucl
      target <- reads_fasta_path
      db_blast <- paste0(pung_path,"/Resources/AllKIRGen_noGapsAlleleNames_gen.fas")
      blast_results <- paste0(dirname(reads_fasta_path),"/", sample_name, "_kir_blastn.txt")
      
      
      ### RUN BLAST ###
      # Filters: 
      # -outfmt 6 Tabular output, no comments
      # -perc_identity 100
      # significance cutoff -evalue 1e-30
      system(sprintf("blastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-30 -perc_identity 90", # 
                     target, db_blast, blast_results))
      
      # Add header to result
      blastn <- readLines(blast_results)
      header <- "Query_id	Subject_id	identity	alignment_length	mismatches	gap_openings	query_start	query_end	subject_start	subject_end	e-value	bit_score_"
      blastn <- c(header,blastn)
      writeLines(blastn, con=blast_results, sep="\n")
    }
    
    
    # Read results
    blastn <- readLines(blast_results)
    blastn <- as.data.frame(do.call(rbind, strsplit(blastn, split="\t")), stringsAsFactors=FALSE)
    names(blastn) <- blastn[1,]
    blastn <- blastn[-1,]
    
    # Filter which alleles match exactly the reads found, with max of one mismatch
    blastn <- subset(blastn, mismatches <=1)
    if(nrow(blastn) == 0){
      cat(paste0("\t",each_snp," failed BLAST - did not align to any gene."), sep="\n")
      #unresolved_data <- remove_snp(unresolved_data, index=i, sample_name, each_snp)
      next
    }
    
    # Filter out 3DL1_059, is it a hybrid allele?? # CONFIRM THIS
    # Import probes for 3DL1_059
    # If probes match fastq files, remove those reads
    if(gene_name == "KIR3DL2"){
      fusion_probes <- read.delim(paste0(pung_path,"/Resources/probes/Fusion_Probes.txt"), header=F, sep="\t", colClasses = "character")
      names(fusion_probes) <- c("name", "probe")
      index_fusion <- grep("3DL12|3DL2",fusion_probes$name )
      fusion_probes <- fusion_probes[index_fusion,]
      # Search probes in each FASTQ file
      fastq_files <- list.files(path = paste0(data_directory,"extractedFastq/"), pattern = sample_name, full.names = TRUE)
      reads <- c()
      for (fastq_file in fastq_files){
        add <- readLines(fastq_file)
        indexes <- grep(paste(fusion_probes$probe, collapse="|"), add)
        # get header line (-1) and quality lines (+2)
        indexes <- c(indexes-1, indexes, indexes+1, indexes+2)
        indexes <- indexes[order(indexes)]
        add <- add[indexes]
        reads <- c(reads, add)
        rm(add)
      }
      # If the fusion gene is not present, remove its blast matches from main gene
      if(isempty(reads)){
        blastn <- subset(blastn, Subject_id != "3DL1_059")
        # Working on reporting this fusion gene result 
      }
    } 
    
    # Table results
    matches <- table(blastn$Subject_id)
    nonspecific_matches <- matches[grep(gsub("^KIR","", gene_name),names(matches), invert = T)]
    
    # save list of blast hits
    blast_df <- data.frame(matches)
    names(blast_df) <- c("Subject_id", "hits")
    filteredFastq_path <- dirname(blast_results)
    write.csv(blast_df, file=paste0(filteredFastq_path, "/",sample_name, "_blast_list.csv"), quote = F, row.names = F)
    
    # If there are misaligned reads, remove snp. If not, go forward
    if (!isempty(nonspecific_matches)){
      misaligned_gene <- unique(gsub("_.*$","",names(nonspecific_matches)))
      cat(paste0("\tThe 'new' variant ", each_snp, " in ",gene_name," matches with sequences in ", paste(misaligned_gene, collapse=" "),". Removing SNP..."), sep="\n")
      unresolved_data <- remove_snp(unresolved_data, index, sample_name, each_snp)
    } else{
      cat(paste0("\t",each_snp," passed BLAST."), sep="\n")
      
    }
  }
  
  
  
  return(unresolved_data)
  
  
}