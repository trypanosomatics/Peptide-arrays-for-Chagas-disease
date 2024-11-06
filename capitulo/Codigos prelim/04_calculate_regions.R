#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. Make sure you have run all previous codes
## 2. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 3. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### CONFIG ####
main_folder <- "." #When running in Rstudio, set this to the absolute path of the folder containing this code
main_folder <- "/Users/romer/Documentos/Peptide-arrays-for-Chagas-disease/Peptide-arrays-for-Chagas-disease"

functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/04_calculate_regions_function.R", functions_folder)

testing <- TRUE #set this to FALSE when running the actual data


#### READ ARGUMENTS AND GET PATH (DO NOT CHANGE) ####
args <- commandArgs(TRUE)

if (length(args == 2)) {
  if (args[1] == "-test") {
    testing <- as.logical(args[2])
  }
}

if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/test_data", main_folder)
  extra_files_folder <- sprintf("%s/extra_files", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/chagastope_data", main_folder)
  extra_files_folder <- sprintf("%s/extra_files", main_folder)
}

#### INTERNAL CONFIG (DO NOT CHANGE) ####
source(functions_file) # Load the auxiliary functions

# File paths
#chronic_peaks_data_file <- sprintf("%s/outputs/03_pools_antigenic_peaks/pools_peaks_cutoff4SD_2pep.tsv", project_folder)
#sylvio_fasta_data_file <- sprintf("%s/sylvio-x10-tritrypdb-5.Dec.2016-no-pseudo_SIMPLIFIED_TABBED.tsv", extra_files_folder)
#brener_fasta_data_file <- sprintf("%s/tcruzi-clbrener-tritrypdb-5.Dec.2016-no-pseudo_SIMPLIFIED_TABBED.tsv", extra_files_folder)

# Parameters
chronic_peaks_tag <- "analyzed_proteome_cutoff_4_SD_2_pep"
chronic_peaks_types <- c("PO")
chronic_peaks_max_other_type_signal_ratio <- 0.2
chronic_peaks_peptide_length <- 16
chronic_peaks_peptide_offset <- 1
chronic_peaks_min_peptide_length <- 12
chronic_peaks_border_peptide_amount <- 16

# Outputs
detailed_region_output_file <- sprintf("%s/outputs/04_antigenic_regions/peptide_region_data.tsv", project_folder)
summary_region_output_file <- sprintf("%s/outputs/04_antigenic_regions/summary_region_data.tsv", project_folder)

#### LOAD DATA ####
#chronic_peaks_data <- fread(chronic_peaks_data_file, header = TRUE, sep = "\t", na.strings = NULL)
#sylvio_fasta_data <- fread(sylvio_fasta_data_file, header = TRUE, sep = "\t", na.strings = NULL)
#brener_fasta_data <- fread(brener_fasta_data_file, header = TRUE, sep = "\t", na.strings = NULL)

# Combine both fasta data
#chronic_all_fasta <- rbindlist(list(sylvio_fasta_data, brener_fasta_data))

#### CALCULATE REGIONS ####
region_data <- processAntigenicRegions(
  project_folder, 
  extra_files_folder,
  chronic_peaks_tag,
  chronic_peaks_types, 
  chronic_peaks_max_other_type_signal_ratio, 
  chronic_peaks_peptide_length, 
  chronic_peaks_peptide_offset, 
  chronic_peaks_min_peptide_length, 
  chronic_peaks_border_peptide_amount, 
  detailed_region_output_file, 
  summary_region_output_file
)
