#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. In chagastope_data/inputs/11_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 2. In chagastope_data/inputs/12_individual_serums_raw_data place each of the raw data files from CHAGASTOPE-v2 (such as AR_P1_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### WARNINGS ####
## This code uses large amounts of RAM

#### CONFIG ####
main_folder <- "." #When running in Rstudio, set this to the absolute path of the folder containing this code
testing <- TRUE #set this to FALSE when running the actual data

#### READ ARGUMENTS AND GET PATH (DO NOT CHANGE) ####
args <- commandArgs(TRUE)

if (length(args == 2)) {
    if (args[1] == "-test") {
        testing <- as.logical(args[2])
    }
}

if (testing == TRUE) {
    #For testing
    project_folder <- sprintf("%s/test_data", main_folder)
} else {
    #For running the actual data
    project_folder <- sprintf("%s/chagastope_data", main_folder)
}

#### INTERNAL CONFIG (DO NOT CHANGE) ####
library(data.table)
library(preprocessCore) #quantile.normalization

design_data_file <- sprintf("%s/inputs/11_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)

raw_data_folder <- sprintf("%s/inputs/12_individual_serums_raw_data", project_folder)

sources <- c("AR_P1", "AR_P2", "AR_P3", "AR_P4", "AR_P5", "AR_P6", "BO_P1", "BO_P2", "BO_P3", "BO_P4", "BO_P5", "BO_P6", "BR_P1", "BR_P2", "BR_P3", "BR_P4", "BR_P5", "CO_P1", "CO_P2",
             "CO_P3", "CO_P4", "MX_P1", "MX_P2", "MX_P3", "MX_P4", "MX_P5", "MX_P6", "US_P1", "US_P2", "US_P3", "US_P4", "US_P5", "US_P6", "AR_E1", "AR_E2", "AR_E3", "AR_E4", "AR_E5",
             "AR_E6", "BO_E1", "BO_E2", "BO_E3", "BO_E4", "BO_E5", "BO_E6", "BR_E1", "BR_E2", "BR_E3", "BR_E4", "BR_E5", "BR_E6", "BR_E7", "CO_E1", "CO_E2", "CO_E3", "CO_E4", "CO_E5",
             "CO_E6", "CO_E7", "MX_E1", "MX_E2", "MX_E3", "MX_E4", "MX_E5", "MX_E6", "US_E1", "US_E2", "US_E3", "US_E4", "US_E5", "US_E6")

types <- c("PO")

output_signal_decimals <- 2
output_statistics_mode_decimals <- 2
output_statistics_sd_decimals <- 2

output_folder <- sprintf("%s/outputs/11_individual_serums_normalized_data", project_folder)
output_suffix <- "_processed.tsv"

#### AUXILIAR FUNCTIONS ####
calculateMode <- function(x, decimals = 0) {
    x <- round(x, decimals)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#### NORMALIZE ALL DATA ####
#Get data
output_initialized <- F
for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
        # type_for <- types[1]
        raw_data_file <- sprintf("%s/%s_%s_raw.tsv", raw_data_folder, source_for, type_for)  
        raw_data <- fread(raw_data_file, header = T, sep = "\t", na.strings = NULL)
     
        setnames(raw_data, "allData.replica1", sprintf("%s_%s_rawData_allData_r1", source_for, type_for))
        setnames(raw_data, "allData.replica2", sprintf("%s_%s_rawData_allData_r2", source_for, type_for))
        
        #Check if I need to initialize the tables to normalize
        if (output_initialized == F) {
            full_raw_PO_data <- raw_data[, .(`Reporter Name`)]
            output_initialized <- T
        }
        
        #Add the data to the corresponding table
        full_raw_PO_data <- merge(full_raw_PO_data,
                                  raw_data,
                                  by = "Reporter Name")
    }
}

#Normalize the positive data
matrix_PO_aux <- normalize.quantiles(as.matrix(full_raw_PO_data[, -c("Reporter Name")]))
matrix_PO_aux <- round(matrix_PO_aux, digits = output_signal_decimals)
normalized_PO_data <- as.data.table(matrix_PO_aux)
colnames(normalized_PO_data) <- colnames(full_raw_PO_data[, -c("Reporter Name")])
normalized_PO_data$`Reporter Name` <- full_raw_PO_data$`Reporter Name`
rm(full_raw_PO_data)
rm(matrix_PO_aux)
gc()

#Combine data
colnames(normalized_PO_data) <- gsub("rawData", "normalizedData", colnames(normalized_PO_data))

full_normalized_data <- normalized_PO_data
rm(normalized_PO_data)
gc()

#### NORMALIZE ONLY REGION DATA ####
design_data <- fread(design_data_file, header = T, sep = "\t", na.strings = NULL)
ids_to_keep <- unique(design_data[group == "antigenic_regions_Tcruzi"]$array_express_id)

#Get data
output_initialized <- F
for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
        # type_for <- types[1]
        raw_data_file <- sprintf("%s/%s_%s_raw.tsv", raw_data_folder, source_for, type_for)  
        raw_data <- fread(raw_data_file, header = T, sep = "\t", na.strings = NULL)
        
        raw_data <- raw_data[`Reporter Name` %in% ids_to_keep]
                 
        setnames(raw_data, "allData.replica1", sprintf("%s_%s_rawData_onlyRegionsData_r1", source_for, type_for))
        setnames(raw_data, "allData.replica2", sprintf("%s_%s_rawData_onlyRegionsData_r2", source_for, type_for))
        
        #Check if I need to initialize the tables to normalize
        if (output_initialized == F) {
            full_raw_PO_data <- raw_data[, .(`Reporter Name`)]
            output_initialized <- T
        }
        
        #Add the data to the corresponding table
        full_raw_PO_data <- merge(full_raw_PO_data,
                                  raw_data,
                                  by = "Reporter Name")
    }
}

#Normalize the positive data
matrix_PO_aux <- normalize.quantiles(as.matrix(full_raw_PO_data[, -c("Reporter Name")]))
matrix_PO_aux <- round(matrix_PO_aux, digits = output_signal_decimals)
normalized_PO_data <- as.data.table(matrix_PO_aux)
colnames(normalized_PO_data) <- colnames(full_raw_PO_data[, -c("Reporter Name")])
normalized_PO_data$`Reporter Name` <- full_raw_PO_data$`Reporter Name`
rm(full_raw_PO_data)
rm(matrix_PO_aux)
gc()

#Combine data
colnames(normalized_PO_data) <- gsub("rawData", "normalizedData", colnames(normalized_PO_data))

full_normalized_data <- merge(full_normalized_data,
                              normalized_PO_data,
                              by = "Reporter Name",
                              all.x = T)

#### CALCULATE GLOBAL STATISTICS (FROM THE ONLY REGION NORMALIZATION) ####
cols_to_select <- setdiff(colnames(normalized_PO_data), "Reporter Name")
global_signals <- as.vector(as.matrix(normalized_PO_data[, cols_to_select, with = F]))

mode_aux <- calculateMode(global_signals, decimals = output_statistics_mode_decimals)
sd_aux <- round(sd(global_signals), output_statistics_sd_decimals)

global_statistics <- data.table(mode = mode_aux,
                                sd = sd_aux)

rm(normalized_PO_data)
rm(global_signals)
gc()

#### SAVE DATA ####
global_statistics_output_file <- sprintf("%s/global_statistics.tsv", output_folder)
write.table(global_statistics, file = global_statistics_output_file, col.names = T, row.names = F, sep = "\t", quote = T)

for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
        # type_for <- types[1]    
        cols_to_select <- colnames(full_normalized_data)
        cols_to_select <- cols_to_select[grepl(sprintf("%s_%s", source_for, type_for), cols_to_select)]
        cols_to_select <- c("Reporter Name", cols_to_select)
        
        normalized_data_for <- full_normalized_data[, cols_to_select, with = F]
        
        setnames(normalized_data_for, sprintf("%s_%s_normalizedData_allData_r1", source_for, type_for), "allData.replica1")
        setnames(normalized_data_for, sprintf("%s_%s_normalizedData_allData_r2", source_for, type_for), "allData.replica2")
        setnames(normalized_data_for, sprintf("%s_%s_normalizedData_onlyRegionsData_r1", source_for, type_for), "onlyRegionsData.replica1")
        setnames(normalized_data_for, sprintf("%s_%s_normalizedData_onlyRegionsData_r2", source_for, type_for), "onlyRegionsData.replica2")
        
        normalized_data_output_file <- sprintf("%s/%s_%s%s", output_folder, source_for, type_for, output_suffix)
        write.table(normalized_data_for, file = normalized_data_output_file, col.names = T, row.names = F, sep = "\t", quote = T, na = "")
    }
}
