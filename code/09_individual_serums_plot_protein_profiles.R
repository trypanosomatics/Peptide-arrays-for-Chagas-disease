#!/usr/bin/env Rscript

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
sources <- c("AR_P1", "AR_P2", "AR_P3", "AR_P4", "AR_P5", "AR_P6", "BO_P1", "BO_P2", "BO_P3", "BO_P4", "BO_P5", "BO_P6", "BR_P1", "BR_P2", "BR_P3", "BR_P4", "BR_P5", "CO_P1", "CO_P2",
             "CO_P3", "CO_P4", "MX_P1", "MX_P2", "MX_P3", "MX_P4", "MX_P5", "MX_P6", "US_P1", "US_P2", "US_P3", "US_P4", "US_P5", "US_P6", "AR_E1", "AR_E2", "AR_E3", "AR_E4", "AR_E5",
             "AR_E6", "BO_E1", "BO_E2", "BO_E3", "BO_E4", "BO_E5", "BO_E6", "BR_E1", "BR_E2", "BR_E3", "BR_E4", "BR_E5", "BR_E6", "BR_E7", "CO_E1", "CO_E2", "CO_E3", "CO_E4", "CO_E5",
             "CO_E6", "CO_E7", "MX_E1", "MX_E2", "MX_E3", "MX_E4", "MX_E5", "MX_E6", "US_E1", "US_E2", "US_E3", "US_E4", "US_E5", "US_E6")
profile_data_suffix <- "smoothed_signals.tsv"
protein <- NULL  
sd_multiplier_for_cutoff <- 4
only_proteins_above <- 0
only_proteins_below <- 0
output_suffix <- NULL
pdf_height <- 11.69
pdf_width <- 8.27

# User values
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
    if (args[i] == "--sources") {
      sources <- unlist(strsplit(args[i + 1], ","))
    }
    if (args[i] == "--profile_data_suffix") {
      profile_data_suffix <- args[i + 1]
      output_suffix <- args[i + 1]
    }
    if (args[i] == "--protein") {
      protein <- args[i + 1]
    }
    if (args[i] == "--sd_multiplier_for_cutoff") {
      sd_multiplier_for_cutoff <- as.numeric(args[i + 1])
    }
    if (args[i] == "--pdf_height") {
      pdf_height <- as.numeric(args[i + 1])
    }
    if (args[i] == "--pdf_width") {
      pdf_width <- as.numeric(args[i + 1])
    }
  }
}

# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Profile data suffix:", profile_data_suffix, "\n")
cat("Protein (optional):", protein, "\n")
cat("SD multiplier for cutoff>", sd_multiplier_for_cutoff, "\n")

#### PATH CONFIG ####
if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
}

input_folder <- sprintf("%s/outputs/06_individual_serums_smoothed_data", project_folder)
design_file <- sprintf("%s/inputs/03_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)
output_folder <- sprintf("%s/outputs/08_09_10_plots", project_folder)

design_data_folder <- sprintf("%s/inputs/03_individual_serums_array_design", project_folder)

if (!testing && length(setdiff(list.files(design_data_folder), ".gitkeep")) == 0) {
  stop("Download the CHAGASTOPE Assay Individual Serums Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)

#### CALL MAIN FUNCTION ####
source(sprintf("%s/09_individual_serums_plot_protein_profiles_function.R", functions_folder))

plot_proteins_ridge(project_folder = project_folder,
              input_folder = input_folder,
              design_file = design_file,
              sources = sources,
              profile_data_suffix = profile_data_suffix,
              output_folder = output_folder,
              protein = protein,
              only_proteins_above = only_proteins_above,
              only_proteins_below = only_proteins_below,
              sd_multiplier_for_cutoff = sd_multiplier_for_cutoff,
              output_suffix = output_suffix,
              pdf_height = pdf_height,
              pdf_width = pdf_width)

cat("Output saved in:", output_folder, "\n")