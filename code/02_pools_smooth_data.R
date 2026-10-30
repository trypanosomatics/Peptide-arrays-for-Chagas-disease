#!/usr/bin/env Rscript
# Shebang, interpreter for executable

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v1"
sources <- NULL # default: all sources of the array (functions/00_common.R)
smoothing_median_window_size <- 3
smoothing_mean_window_size <- 0
# "zeros" fill the sides with 0, "repeat" fill the side with signals equal to the last one
smooth_borders_option <- "zeros"


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
    if (args[i] == "--array") {
      array_id <- args[i + 1]
    }
  }
}

source(sprintf("%s/code/functions/00_common.R", main_folder))
array_config <- get_array(array_id, supported = "v1")
if (is.null(sources)) {
  sources <- array_config$default_sources
}

# Print values 
cat("Main folder:", main_folder, "\n")
cat("Array:", array_config$id, "\n")
cat("Testing:", testing, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Smoothing median window size:", smoothing_median_window_size, "\n")
cat("Smoothing mean window size:", smoothing_mean_window_size, "\n")
cat("Smooth borders option:", smooth_borders_option, "\n")
cat("Smoothing label (use it in later steps):", smoothing_label(smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option), "\n")

# CONFIG
setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/data/test_data", main_folder)
} else {
  sprintf("%s/data/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
design_data_file <- sprintf("%s/%s", project_folder, array_config$design_file)
design_data_folder <- dirname(design_data_file)
output_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)

if (!testing && (length(setdiff(list.files(raw_data_folder), ".gitkeep")) == 0 || length(setdiff(list.files(design_data_folder), ".gitkeep")) == 0)) {
  stop("Download the CHAGASTOPE Assay Pool Data and CHAGASTOPE Assay Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}

functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/02_pools_smooth_data_function.R", functions_folder)

normalized_data_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#### CALL MAIN FUNCTION ####
source(functions_file)

output_signal_mean_decimals <- 2
output_signal_sd_decimals <- 2

smooth_data(array_config, design_data_file, normalized_data_folder, output_folder, sources, smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option, output_signal_mean_decimals, output_signal_sd_decimals)

cat("Output saved in:", output_folder, "\n")