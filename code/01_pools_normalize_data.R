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

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")

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

setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/data/test_data", main_folder)
} else {
  sprintf("%s/data/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
output_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

if (!testing && length(list.files(raw_data_folder)) == 0) {
  stop("Download the CHAGASTOPE Assay Pool Data and CHAGASTOPE Assay Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}


if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/01_pools_normalize_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

normalize_and_process_data(raw_data_folder, output_folder, sources)