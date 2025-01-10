#!/usr/bin/env Rscript

#### PREPARATION STEPS ####
## You can run this code as it is and reproduce Ag2 (TcCLB.511671.60) Alanine Scan results, or you can follow these next steps for the entire dataset.

## 1. In chagastope_data/inputs/11_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 2. In chagastope_data/inputs/12_individual_serums_raw_data place each of the raw data files from CHAGASTOPE-v2 (such as AR_P1_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### WARNINGS ####
## This code uses large amounts of RAM


#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- "."
testing <- TRUE


for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--main_folder") {
    main_folder <- args[i + 1]
  }
  if (args[i] == "--testing") {
    testing <- as.logical(args[i + 1])
  }
  if (args[i] == "--selected_protein") {
    selected_protein <- paste0("\"", args[i + 1], "\"")
  }
}

#### OTHER PARAMETERS ####  


#### PATH CONFIG ####

if (testing == TRUE) {
  #For testing
  project_folder <- sprintf("%s/test_data", main_folder)
  design_alanine_file <- sprintf("%s/extra_files/alanine_scan_design_test_subset.tsv", main_folder)  
} else {
  #For running the actual data
  project_folder <- sprintf("%s/chagastope_data", main_folder)
  design_alanine_file <- sprintf("%s/extra_files/alanine_scan_design.tsv", main_folder)
}

design_data_file <- sprintf("%s/inputs/11_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)


#### Aux function path ####
functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/07_alanine_scan_analysis_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

alanine_scan(main_folder, testing, selected_protein)
