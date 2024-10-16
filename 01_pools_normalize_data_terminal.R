#!/usr/bin/env Rscript
# Shebang, interpreter for executable

#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. In chagastope_data/inputs/02_pools_raw_data place each of the raw data files from CHAGASTOPE-v1 (such as AR_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11651/)
## 2. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 3. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### WARNINGS ####
## This code uses large amounts of RAM


#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  main_folder <- "." #Default value
  testing <- TRUE #Default value
  sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE") #Default value
} else {
  main_folder <- args[1] # first argument
  testing <- as.logical(args[2]) # second argument (TRUE o FALSE)
  sources <- unlist(strsplit(args[3], ",")) # third argument (comma-separated list of sources)
}

setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/test_data", main_folder)
} else {
  sprintf("%s/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
output_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/01_pools_normalize_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

normalize_and_process_data(raw_data_folder, output_folder, sources)