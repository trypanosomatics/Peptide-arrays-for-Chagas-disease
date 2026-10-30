#!/usr/bin/env Rscript

#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v1" # v1 (pools) or v2 (individual sera)
# smoothed data to plot, as printed by 02_pools_smooth_data.R / 06_individual_serums_smooth_data.R
smoothing <- "med3_mean0_zeros"
y_factors <- c(10000, 2, 0.5)

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
    if (args[i] == "--array") {
      array_id <- args[i + 1]
    }
    if (args[i] == "--smoothing_label") {
      smoothing <- args[i + 1]
    }
    if (args[i] == "--y_factors") {
      y_factors <- as.numeric(unlist(strsplit(args[i + 1], ",")))
    }
  }
}

source(sprintf("%s/code/functions/00_common.R", main_folder))
array_config <- get_array(array_id, supported = c("v1", "v2"))

# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Array:", array_config$id, "\n")
cat("Smoothing label:", smoothing, "\n")
cat("Y factors:", y_factors, "\n")

#### PATH CONFIG ####
project_folder <- if (testing) {
  sprintf("%s/data/test_data", main_folder)
} else {
  sprintf("%s/data/chagastope_data", main_folder)
}

output_folder <- sprintf("%s/outputs/08_09_10_plots", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/10_plot_signal_variations_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

result <- plot_density(array_config, project_folder, output_folder, smoothing, y_factors)

cat("Output saved in:", output_folder, "\n")
