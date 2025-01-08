#!/usr/bin/env Rscript

#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

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
sources <- c("AR_P1", "AR_P2", "AR_P3", "AR_P4", "AR_P5", "AR_P6", "BO_P1", "BO_P2", "BO_P3", "BO_P4", "BO_P5", "BO_P6", "BR_P1", "BR_P2", "BR_P3", "BR_P4", "BR_P5", "CO_P1", "CO_P2",
             "CO_P3", "CO_P4", "MX_P1", "MX_P2", "MX_P3", "MX_P4", "MX_P5", "MX_P6", "US_P1", "US_P2", "US_P3", "US_P4", "US_P5", "US_P6", "AR_E1", "AR_E2", "AR_E3", "AR_E4", "AR_E5",
             "AR_E6", "BO_E1", "BO_E2", "BO_E3", "BO_E4", "BO_E5", "BO_E6", "BR_E1", "BR_E2", "BR_E3", "BR_E4", "BR_E5", "BR_E6", "BR_E7", "CO_E1", "CO_E2", "CO_E3", "CO_E4", "CO_E5",
             "CO_E6", "CO_E7", "MX_E1", "MX_E2", "MX_E3", "MX_E4", "MX_E5", "MX_E6", "US_E1", "US_E2", "US_E3", "US_E4", "US_E5", "US_E6")



for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--main_folder") {
    main_folder <- args[i + 1]
  }
  if (args[i] == "--testing") {
    testing <- as.logical(args[i + 1])
  }
  if (args[i] == "--sources") {
    sources <- unlist(strsplit(args[i + 1], ","))
  }
}

#### OTHER PARAMETERS ####  
types <- c("PO") #Positive

output_signal_decimals <- 2
output_statistics_mode_decimals <- 2
output_statistics_sd_decimals <- 2


#### PATH CONFIG ####
if (testing == TRUE) {
  #For testing
  project_folder <- sprintf("%s/test_data", main_folder)
} else {
  #For running the actual data
  project_folder <- sprintf("%s/chagastope_data", main_folder)
}

design_data_file <- sprintf("%s/inputs/11_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)

raw_data_folder <- sprintf("%s/inputs/12_individual_serums_raw_data", project_folder)

output_folder <- sprintf("%s/outputs/11_individual_serums_normalized_data", project_folder)
output_suffix <- "_processed.tsv"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#### Aux function path ####
functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/05_individual_serums_normalize_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

normalize_serums(raw_data_folder, output_folder, testing, sources, types, output_suffix, output_signal_decimals, output_statistics_mode_decimals, output_statistics_sd_decimals)
