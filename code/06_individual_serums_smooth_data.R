#!/usr/bin/env Rscript

#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. Make sure you have run all previous codes
## 2. In chagastope_data/inputs/03_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "--testing FALSE" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v2"
sources <- NULL # default: all sources of the array (functions/00_common.R)
smoothing_median_window_size <- 3
smoothing_mean_window_size <- 0
# smooth_borders_option is used to have the same amount of data after the smoothing.
# "zeros" fill the sides with 0, "repeat" fill the side with signals equal to the last one
smooth_borders_option <- "zeros" 

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
array_config <- get_array(array_id, supported = "v2")
if (is.null(sources)) {
  sources <- array_config$default_sources
}

#### OTHER PARAMETERS ####  

output_signal_mean_decimals <- 2
output_signal_sd_decimals <- 2


# Print values 
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Array:", array_config$id, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Smoothing median window size:", smoothing_median_window_size, "\n")
cat("Smoothing mean window size:", smoothing_mean_window_size, "\n")
cat("Smooth borders option:", smooth_borders_option, "\n")
cat("Types:", paste(array_config$types, collapse = ", "), "\n")
cat("Smoothing label (use it in later steps):", smoothing_label(smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option), "\n")


#### PATH CONFIG ####
if (testing == TRUE) {
  #For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
} else {
  #For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
}

design_data_file <- sprintf("%s/%s", project_folder, array_config$design_file)
design_data_folder <- dirname(design_data_file)

if (!testing && length(setdiff(list.files(design_data_folder), ".gitkeep")) == 0) {
  stop("Download the CHAGASTOPE Assay Individual Serums Data and CHAGASTOPE Assay Individual Serums Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}

normalized_data_folder <- sprintf("%s/outputs/05_individual_serums_normalized_data", project_folder)

output_folder <- sprintf("%s/outputs/06_individual_serums_smoothed_data", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


#### Aux function path ####
functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/06_individual_serums_smooth_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

smooth_serums(array_config, design_data_file, normalized_data_folder, output_folder, sources, smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option,
              output_signal_mean_decimals, output_signal_sd_decimals)

cat("Output saved in:", output_folder, "\n")