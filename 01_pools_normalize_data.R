#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. In chagastope_data/inputs/02_pools_raw_data place each of the raw data files from CHAGASTOPE-v1 (such as AR_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11651/)
## 2. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 3. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

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

raw_data_folder <- sprintf("%s/inputs/02_pools_raw_data", project_folder)

sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")
sources_to_normalize_as_negatives <- c("LE")
sources_to_exclude_from_global_statistics <- c("LE")

types <- c("PO", "NE")

output_signal_decimals <- 2
output_statistics_mode_decimals <- 2
output_statistics_sd_decimals <- 2

output_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)
output_suffix <- "_processed.tsv"

#### AUXILIAR FUNCTIONS ####
calculateMode <- function(x, decimals = 0) {
    x <- round(x, decimals)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#### NORMALIZE DATA ####
#Get data
output_initialized <- F
for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
        # type_for <- types[1]
        raw_data_file <- sprintf("%s/%s_%s_raw.tsv", raw_data_folder, source_for, type_for)  
        raw_data <- fread(raw_data_file, header = T, sep = "\t", na.strings = NULL)
     
        #Check if I need to initialize the tables to normalize
        if (output_initialized == F) {
            full_raw_PO_data <- raw_data[, .(`Reporter Name`)]
            full_raw_NE_data <- raw_data[, .(`Reporter Name`)]
            output_initialized <- T
        }
        
        #Add the data to the corresponding table
        if ((type_for == "NE") | (source_for %in% sources_to_normalize_as_negatives)) {
            full_raw_NE_data <- merge(full_raw_NE_data,
                                      raw_data,
                                      by = "Reporter Name")
        } else {
            full_raw_PO_data <- merge(full_raw_PO_data,
                                      raw_data,
                                      by = "Reporter Name")
        }
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

#Normalize the negative data
matrix_NE_aux <- normalize.quantiles(as.matrix(full_raw_NE_data[, -c("Reporter Name")]))
matrix_NE_aux <- round(matrix_NE_aux, digits = output_signal_decimals)
normalized_NE_data <- as.data.table(matrix_NE_aux)
colnames(normalized_NE_data) <- colnames(full_raw_NE_data[, -c("Reporter Name")])
normalized_NE_data$`Reporter Name` <- full_raw_NE_data$`Reporter Name`
rm(full_raw_NE_data)
rm(matrix_NE_aux)
gc()

#Combine data
full_normalized_data <- merge(normalized_PO_data,
                              normalized_NE_data,
                              by = "Reporter Name")
rm(normalized_PO_data)
rm(normalized_NE_data)
gc()

colnames(full_normalized_data) <- gsub("rawData", "normalizedData", colnames(full_normalized_data))

#### CALCULATE GLOBAL STATISTICS ####
cols_to_select <- setdiff(colnames(full_normalized_data), "Reporter Name")
if (length(sources_to_exclude_from_global_statistics) > 0) {
    for (source_for in sources_to_exclude_from_global_statistics) {
        # source_for <- sources_to_exclude_from_global_statistics[1]
        cols_to_select <- cols_to_select[!grepl(sprintf("%s_", source_for), cols_to_select)]
    }
}

global_signals <- as.vector(as.matrix(full_normalized_data[, cols_to_select, with = F]))

mode_aux <- calculateMode(global_signals, decimals = output_statistics_mode_decimals)
sd_aux <- round(sd(global_signals), output_statistics_sd_decimals)

global_statistics <- data.table(mode = mode_aux,
                                sd = sd_aux)
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
        
        normalized_data_output_file <- sprintf("%s/%s_%s%s", output_folder, source_for, type_for, output_suffix)
        write.table(normalized_data_for, file = normalized_data_output_file, col.names = T, row.names = F, sep = "\t", quote = T) 
    }
}
