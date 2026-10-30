#!/usr/bin/env Rscript

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
array_id <- "v1"
sources <- NULL # default: all sources of the array (functions/00_common.R)
# smoothed data to plot, as printed by 02_pools_smooth_data.R
smoothing <- "med3_mean0_zeros"
protein <- NULL
sd_multiplier_for_cutoff <- 4
only_proteins_above <- 0
only_proteins_below <- 0
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
    if (args[i] == "--array") {
      array_id <- args[i + 1]
    }
    if (args[i] == "--sources") {
      sources <- unlist(strsplit(args[i + 1], ","))
    }
    if (args[i] == "--smoothing_label") {
      smoothing <- args[i + 1]
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
    if (args[i] == "--fixed_scale") {
      fixed_scale <- as.logical(args[i + 1])
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
cat("Testing:", testing, "\n")
cat("Array:", array_config$id, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Smoothing label:", smoothing, "\n")
cat("Protein (optional):", protein, "\n")
cat("SD multiplier for cutoff:", sd_multiplier_for_cutoff, "\n")

#### PATH CONFIG ####
if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
}

input_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)
output_folder <- sprintf("%s/outputs/08_09_10_plots", project_folder)

design_data_folder <- dirname(sprintf("%s/%s", project_folder, array_config$design_file))

if (!testing && length(setdiff(list.files(design_data_folder), ".gitkeep")) == 0) {
  stop("Download the CHAGASTOPE Assay Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)

#### CALL MAIN FUNCTION ####
source(sprintf("%s/08_pools_plot_protein_profiles_function.R", functions_folder))

plot_proteins(array_config = array_config,
              project_folder = project_folder,
              input_folder = input_folder,
              sources = sources,
              smoothing = smoothing,
              output_folder = output_folder,
              protein = protein,
              only_proteins_above = only_proteins_above,
              only_proteins_below = only_proteins_below,
              sd_multiplier_for_cutoff = sd_multiplier_for_cutoff,
              fixed_scale = fixed_scale)

cat("Output saved in:", output_folder, "\n")
