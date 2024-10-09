#### CONFIG ####
main_folder <- "." #When running in Rstudio, set this to the absolute path of the folder containing this code
setwd(main_folder)

testing <- TRUE #set this to TRUE for testing purposes
sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")

project_folder <- if (testing) {
  sprintf("%s/test_data", main_folder)
} else {
  sprintf("%s/chagastope_data", main_folder)
}

#### PATH CONFIG ####
raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)
output_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

functions_folder <- sprintf("%s/functions", main_folder)
functions_file <- sprintf("%s/01_pools_normalize_data_function.R", functions_folder)

#### CALL MAIN FUNCTION ####
source(functions_file)

normalize_and_process_data(raw_data_folder, output_folder, sources)