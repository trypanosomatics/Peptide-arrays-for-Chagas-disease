#!/usr/bin/env Rscript

#### CONFIG for running in terminal ####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- "."
testing <- FALSE
peaks_tag <- "4SD_2pep"

# User values
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--main_folder") {
    main_folder <- args[i + 1]
  }
  if (args[i] == "--testing") {
    testing <- as.logical(args[i + 1])
  }
  if (args[i] == "--peaks_tag") {
    peaks_tag <- args[i + 1]
  }
}

output_suffix <- peaks_tag

#### OTHER PARAMETERS#### 
#In this design, the length of the peptides is 16 residues and the overlap is 12, it cannot be changed
# If you use another microarray design, you can change it
chronic_peaks_max_other_type_signal_ratio <- 0.2
chronic_peaks_peptide_length <- 16
chronic_peaks_peptide_offset <- 4
chronic_peaks_min_peptide_length <- 12
chronic_peaks_border_peptide_amount <- 16
chronic_peaks_types <- c("PO")

# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Chronic peaks tag:", peaks_tag, "\n")
cat("Chronic peaks types:", paste(chronic_peaks_types, collapse = ", "), "\n")
cat("Chronic peaks max other type signal ratio:", chronic_peaks_max_other_type_signal_ratio, "\n")
cat("Chronic peaks peptide length:", chronic_peaks_peptide_length, "\n")
cat("Chronic peaks peptide offset:", chronic_peaks_peptide_offset, "\n")
cat("Chronic peaks min peptide length:", chronic_peaks_min_peptide_length, "\n")
cat("Chronic peaks border peptide amount:", chronic_peaks_border_peptide_amount, "\n")
cat("Suffix file:", output_suffix, "\n")


####PATH CONFIG####
if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/test_data", main_folder)
  extra_files_folder <- sprintf("%s/extra_files", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/chagastope_data", main_folder)
  extra_files_folder <- sprintf("%s/extra_files", main_folder)
}

output_folder <- sprintf("%s/outputs/04_antigenic_regions", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

detailed_region_output_file <- sprintf("%s/outputs/04_antigenic_regions/%s_peptide_region_data.tsv", project_folder, output_suffix)
summary_region_output_file <- sprintf("%s/outputs/04_antigenic_regions/%s_summary_region_data.tsv", project_folder, output_suffix)


#### Aux function path ####
functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/04_calculate_regions_function.R", functions_folder)

#### CALL MAIN FUNCTION####
source(functions_file)


#### CALCULATE REGIONS ####
region_data <- processAntigenicRegions(
  project_folder, 
  extra_files_folder,
  peaks_tag,
  chronic_peaks_types, 
  chronic_peaks_max_other_type_signal_ratio, 
  chronic_peaks_peptide_length, 
  chronic_peaks_peptide_offset, 
  chronic_peaks_min_peptide_length, 
  chronic_peaks_border_peptide_amount, 
  detailed_region_output_file, 
  summary_region_output_file
)
