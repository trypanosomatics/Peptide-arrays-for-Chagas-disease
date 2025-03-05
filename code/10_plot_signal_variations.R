#!/usr/bin/env Rscript

#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- TRUE
type_data <- "pools" #or "individual_serums"
raw_data_suffix <- "raw\\.tsv$"
normalized_data_suffix <- "processed\\.tsv$"
smoothed_data_suffix <- "_smoothed_signals\\.tsv$"
y_factors <- c(10000, 2, 0.5)

for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--main_folder") {
    main_folder <- args[i + 1]
  }
  if (args[i] == "--testing") {
    testing <- as.logical(args[i + 1])
  }
  if (args[i] == "--type_data") {
    type_data <- args[i + 1]
  }
  if (args[i] == "--raw_data_suffix") {
    raw_data_suffix <- args[i + 1]
  }
  if (args[i] == "--normalized_data_suffix") {
    normalized_data_suffix <- args[i + 1]
  }
  if (args[i] == "--smoothed_data_suffix") {
    smoothed_data_suffix <- args[i + 1]
  }
  if (args[i] == "--y_transformation") {
    y_transformation <- as.logical(args[i + 1])
  }
  if (args[i] == "--y_factors") {
    y_factors <- as.numeric(unlist(strsplit(args[i + 1], ",")))
  }
}

#### PATH CONFIG ####
setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/data/test_data", main_folder)
} else {
  sprintf("%s/data/chagastope_data", main_folder)
}

output_folder <- sprintf("%s/outputs/08_plots", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/10_plot_signal_variations_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

result <- plot_density(project_folder, type_data, raw_data_suffix, 
                       normalized_data_suffix, smoothed_data_suffix, y_factors, plot_data)
