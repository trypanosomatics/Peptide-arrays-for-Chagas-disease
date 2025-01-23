#!/usr/bin/env Rscript
# Shebang, interpreter for executable

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- "."
testing <- FALSE
sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")
smoothing_median_window_size <- 3
smoothing_mean_window_size <- 0
smooth_borders_option <- "zeros"
output_suffix <- ""


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
      sources <- unlist(strsplit(args[i + 1], ",")) #
    }
    if (args[i] == "--smoothing_median_window_size") {
      smoothing_median_window_size <- as.numeric(args[i + 1])
    }
    if (args[i] == "--smoothing_mean_window_size") {
      smoothing_mean_window_size <- as.numeric(args[i + 1])
    }
    if (args[i] == "--smooth_borders_option") {
      smooth_borders_option <- args[i + 1]
    }
    if (args[i] == "--output_suffix") {
      output_suffix <- paste0(args[i + 1], "_")
    }
  }
}

# Print values 
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Smoothing median window size:", smoothing_median_window_size, "\n")
cat("Smoothing mean window size:", smoothing_mean_window_size, "\n")
cat("Smooth borders option:", smooth_borders_option, "\n")
cat("Output suffix:", output_suffix, "\n")

# CONFIG
setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/test_data", main_folder)
} else {
  sprintf("%s/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
design_data_folder <- sprintf("%s/inputs/01_pools_array_design", project_folder)
output_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)

functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/02_pools_smooth_data_function.R", functions_folder)

design_data_file <- sprintf("%s/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv", design_data_folder)
design_group <- "Trypanosoma_cruzi"

normalized_data_folder <- sprintf("%s/outputs/01_pools_normalized_data/", project_folder)
normalized_data_suffix <- "_processed.tsv"

#### CALL MAIN FUNCTION ####
source(functions_file)

output_type_order <- c("PO", "NE")
output_signal_mean_decimals <- 2
output_signal_sd_decimals <- 2

smooth_data(design_data_file, design_group, normalized_data_folder, normalized_data_suffix, sources, smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option, output_type_order, output_signal_mean_decimals, output_signal_sd_decimals, output_suffix)
