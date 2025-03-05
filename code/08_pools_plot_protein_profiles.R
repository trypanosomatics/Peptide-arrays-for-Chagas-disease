#!/usr/bin/env Rscript

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")
profile_data_suffix <- "smoothed_signals.tsv"
protein <- NULL  
sd_multiplier_for_cutoff <- 4
only_proteins_above <- 0
only_proteins_below <- 0
output_suffix <- NULL
fixed_scale <- TRUE

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
    if (args[i] == "--only_proteins_above") {
      only_proteins_above <- as.numeric(args[i + 1])
    }
    if (args[i] == "--only_proteins_below") {
      only_proteins_below <- as.numeric(args[i + 1])
    }
    if (args[i] == "--output_suffix") {
      output_suffix <- args[i + 1]
    }
    if (args[i] == "--fixed_scale") {
      fixed_scale <- args[i + 1]
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

input_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)
design_file <- sprintf("%s/inputs/01_pools_array_design/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv", project_folder)
output_folder <- sprintf("%s/outputs/08_plots", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)

#### CALL MAIN FUNCTION ####
source(sprintf("%s/08_pools_plot_protein_profiles_function.R", functions_folder))

plot_proteins(project_folder = project_folder,
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
              fixed_scale = fixed_scale)

