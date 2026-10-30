#!/usr/bin/env Rscript

#### CONFIG for running in terminal ####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v1"
# peaks to use: same smoothing label and cutoff parameters as in 03_calculate_peaks.R
smoothing <- "med3_mean0_zeros"
sd_multiplier_for_cutoff <- 4
min_num_of_peptides_in_peak <- 2

if (length(args) == 0) {
  cat("No arguments provided. Using default values.\n")
} else {
# User values
  for (i in seq(1, length(args), by = 2)) {
    if (args[i] == "--main_folder") {
      main_folder <- args[i + 1]
    }
    if (args[i] == "--testing") {
      testing <- as.logical(args[i + 1])
    }
    if (args[i] == "--smoothing_label") {
      smoothing <- args[i + 1]
    }
    if (args[i] == "--sd_multiplier_for_cutoff") {
      sd_multiplier_for_cutoff <- as.numeric(args[i + 1])
    }
    if (args[i] == "--min_num_of_peptides_in_peak") {
      min_num_of_peptides_in_peak <- as.numeric(args[i + 1])
    }
    if (args[i] == "--array") {
      array_id <- args[i + 1]
    }
  }
}

source(sprintf("%s/code/functions/00_common.R", main_folder))
array_config <- get_array(array_id, supported = "v1")
peaks_tag <- cutoff_label(sd_multiplier_for_cutoff, min_num_of_peptides_in_peak)

#### OTHER PARAMETERS#### 
# Peptide length and offset come from the array config (functions/00_common.R)
chronic_peaks_max_other_type_signal_ratio <- 0.2
chronic_peaks_min_peptide_length <- 12
chronic_peaks_border_peptide_amount <- 16
chronic_peaks_types <- c("PO")

# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Array:", array_config$id, "\n")
cat("Smoothing label:", smoothing, "\n")
cat("Chronic peaks tag:", peaks_tag, "\n")
cat("Chronic peaks types:", paste(chronic_peaks_types, collapse = ", "), "\n")
cat("Chronic peaks max other type signal ratio:", chronic_peaks_max_other_type_signal_ratio, "\n")
cat("Peptide length / offset:", array_config$peptide_length, "/", array_config$peptide_offset, "\n")
cat("Chronic peaks min peptide length:", chronic_peaks_min_peptide_length, "\n")
cat("Chronic peaks border peptide amount:", chronic_peaks_border_peptide_amount, "\n")


####PATH CONFIG####
if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
}

output_folder <- sprintf("%s/outputs/04_antigenic_regions", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

detailed_region_output_file <- output_file_name(output_folder, array_config$id, "regions_peptides", smoothing = smoothing, cutoff = peaks_tag)
summary_region_output_file <- output_file_name(output_folder, array_config$id, "regions_summary", smoothing = smoothing, cutoff = peaks_tag)


#### Aux function path ####
functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/04_calculate_regions_function.R", functions_folder)

#### CALL MAIN FUNCTION####
source(functions_file)


#### CALCULATE REGIONS ####
region_data <- processAntigenicRegions(
  array_config,
  project_folder, 
  extra_files_folder,
  smoothing,
  peaks_tag,
  chronic_peaks_types, 
  chronic_peaks_max_other_type_signal_ratio, 
  chronic_peaks_min_peptide_length, 
  chronic_peaks_border_peptide_amount, 
  detailed_region_output_file, 
  summary_region_output_file
)

cat("Output saved in:", output_folder, "\n")