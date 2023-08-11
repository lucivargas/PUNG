## PUNG
## Description: Bioinformatic algorithm that, using PING output files, will determine whether the 
## unresolved genotypes from PING may be potential new alleles.
## Author: Vargas, L. B.
## Year: 2021

# Location of PING results, change accordingly
data_directory <- "/home/projects/normanlab_shared/temp/IHW_J/"
results_directory <- "~/PUNG_outputs/IHW_J-PUNG/"

##############################
#### START OF RUN
##############################

# Setting work location
setwd(data_directory)

# Loading libraries 
library("rstudioapi") #getActiveDocumentContext

# Get location of PUNG folder
pung_path <- dirname(getActiveDocumentContext()$path)

# Check if files and permissions needed are available, create results folder
source(paste0(pung_path,"/Functions/check_files.R"))
check_files(data_directory, results_directory, pung_path)

############
## ANALYSES
############

## Validate new alles
# QC, BLAST and check if its been recently deposited on IPD
source(paste0(pung_path,"/Functions/validate_new.R"))
validate_new(data_directory, results_directory, pung_path, 
             kir_filter= c("KIR3DL3", "KIR2DS2", "KIR2DL1", "KIR2DP1", "KIR3DP1",
                           "KIR2DL4", "KIR2DS1", "KIR2DS4", "KIR3DL2"))


