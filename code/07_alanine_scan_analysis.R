#!/usr/bin/env Rscript

#### PREPARATION STEPS ####
## You can run this code as it is and reproduce Ag2 (TcCLB.511671.60) Alanine Scan results, or you can follow these next steps for the entire dataset.

## 1. In chagastope_data/inputs/03_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 2. In chagastope_data/inputs/04_individual_serums_raw_data place each of the raw data files from CHAGASTOPE-v2 (such as AR_P1_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### WARNINGS ####
## This code uses large amounts of RAM


#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- dirname(getwd())
testing <- FALSE
selected_protein <- NULL
heatmap_plot <- TRUE
sequence_logo_plot <- TRUE
sequence_logo_source <- "all"

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
    if (args[i] == "--selected_protein") {
      selected_protein <- unlist(strsplit(args[i + 1], ","))
    }
    if (args[i] == "--heatmap_plot") {
      heatmap_plot <- as.logical(args[i + 1])
    }
    if (args[i] == "--sequence_logo_plot") {
      sequence_logo_plot <- as.logical(args[i + 1])
    }
    if (args[i] == "--sequence_logo_source") {
      if (args[i + 1] == "all") {
        sequence_logo_source <- "all"
      } else {
        sequence_logo_source <- unlist(strsplit(args[i + 1], ","))
      }
    }
  }
}
  
# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Selected_protein:", paste(selected_protein, collapse = ", "), "\n")
cat("Heatmap_plot:", heatmap_plot, "\n")
cat("Sequence_logo_plot:", sequence_logo_plot, "\n")
cat("Sequence_logo_source", paste(sequence_logo_source, collapse = ", "), "\n")


#### PATH CONFIG ####

if (testing == TRUE) {
  #For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
  design_alanine_file <- sprintf("%s/data/extra_files/alanine_scan_design_test_subset.tsv", main_folder)  
} else {
  #For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
  design_alanine_file <- sprintf("%s/data/extra_files/alanine_scan_design.tsv", main_folder)
}

design_data_folder <- sprintf("%s/inputs/03_individual_serums_array_design", project_folder)
design_data_file <- sprintf("%s/inputs/03_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)

if (!testing && length(setdiff(list.files(design_data_folder), ".gitkeep")) == 0) {
  stop("Download the CHAGASTOPE CHAGASTOPE Assay Individual Serums Design Data to perform this operation. Or use the test subset data downloaded with this repository using --testing TRUE")
}

#### Aux function path ####
functions_folder <- sprintf("%s/code/functions", main_folder)
functions_file <- sprintf("%s/07_alanine_scan_analysis_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

alanine_scan(main_folder, testing, selected_protein, heatmap_plot, sequence_logo_plot, sequence_logo_source)

cat("Output saved in:", output_folder, "\n")