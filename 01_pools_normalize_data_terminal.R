#### CONFIG for running in terminal####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  main_folder <- "C:/Users/Ramiro/Documents/GitHub/Peptide-arrays-for-Chagas-disease" #valor por defecto
  testing <- TRUE #valor por defecto
  sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE") #valor por defecto
} else {
  main_folder <- args[1] # primer argumento
  testing <- as.logical(args[2]) # segundo argumento (TRUE o FALSE)
  sources <- unlist(strsplit(args[3], ",")) # tercer argumento (lista de sources separados por comas)
}

setwd(main_folder)

project_folder <- if (testing) {
  sprintf("%s/test_data", main_folder)
} else {
  sprintf("%s/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
output_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/01_pools_normalize_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

normalize_and_process_data(raw_data_folder, output_folder, sources)