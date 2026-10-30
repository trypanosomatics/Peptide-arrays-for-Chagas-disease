#!/usr/bin/env Rscript

#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. Make sure you have run all previous codes
## 2. Set the "testing" variable in the config below to FALSE, or run this code with the "--testing FALSE" argument
## 3. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code


#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v1"
sources <- c("AR", "BO", "BR", "CO", "MX", "US")
min_num_of_peptides_in_peak <- 2
# sd_multiplier_for_cutoff determines cutoff as: mode + sd_multiplier_for_cutoff * sd
sd_multiplier_for_cutoff <- 4
# smoothed data to use, as printed by 02_pools_smooth_data.R
smoothing <- "med3_mean0_zeros"

if (length(args) == 0) {
  cat("No arguments provided. Using default values.\n")
} else {
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
    if (args[i] == "--min_num_of_peptides_in_peak") {
      min_num_of_peptides_in_peak <- as.numeric(args[i + 1])
    }
    if (args[i] == "--sd_multiplier_for_cutoff") {
      sd_multiplier_for_cutoff <- as.numeric(args[i + 1])
    }
    if (args[i] == "--smoothing_label") {
      smoothing <- args[i + 1]
    }
    if (args[i] == "--array") {
      array_id <- args[i + 1]
    }
  }
}

source(sprintf("%s/code/functions/00_common.R", main_folder))
array_config <- get_array(array_id, supported = "v1")

# Print values 
cat("Main folder:", main_folder, "\n")
cat("Array:", array_config$id, "\n")
cat("Testing:", testing, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Min num of peptides in peak:", min_num_of_peptides_in_peak, "\n")
cat("Sd multiplier for cutoff:", sd_multiplier_for_cutoff, "\n")
cat("Smoothing label:", smoothing, "\n")
cat("Peptide length / offset:", array_config$peptide_length, "/", array_config$peptide_offset, "\n")


#### PATH CONFIG ####
if (testing == TRUE) {
  project_folder <- sprintf("%s/data/test_data", main_folder)
} else {
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
}

output_folder <- sprintf("%s/outputs/03_pools_antigenic_peaks", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

profile_data_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)

if (!dir.exists(profile_data_folder)) {
  dir.create(profile_data_folder, recursive = TRUE)
}

#### Aux function path ####
functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/03_calculate_peaks_function.R", functions_folder)

#### OTHER PARAMETERS ####  
# Types and peptide length/offset come from the array config (functions/00_common.R)
other_type_proportion_decimals <- 2
combined_mean_signal_decimals <- 2


#### CALL MAIN FUNCTION ####
source(functions_file)

determine_peaks(array_config, project_folder, profile_data_folder, output_folder, sources, min_num_of_peptides_in_peak, sd_multiplier_for_cutoff, smoothing,
                other_type_proportion_decimals, combined_mean_signal_decimals)

cat("Output saved in:", output_folder, "\n")