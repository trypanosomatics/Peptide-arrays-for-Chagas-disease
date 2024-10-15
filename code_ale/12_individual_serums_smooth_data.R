#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. Make sure you have run all previous codes
## 2. In chagastope_data/inputs/11_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

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
library(zoo) #rollmean rollmedian

design_data_file <- sprintf("%s/inputs/11_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)
design_groups <- c("antigenic_regions_Tcruzi")

normalized_data_folder <- sprintf("%s/outputs/11_individual_serums_normalized_data", project_folder)

sources <- c("AR_P1", "AR_P2", "AR_P3", "AR_P4", "AR_P5", "AR_P6", "BO_P1", "BO_P2", "BO_P3", "BO_P4", "BO_P5", "BO_P6", "BR_P1", "BR_P2", "BR_P3", "BR_P4", "BR_P5", "CO_P1", "CO_P2",
             "CO_P3", "CO_P4", "MX_P1", "MX_P2", "MX_P3", "MX_P4", "MX_P5", "MX_P6", "US_P1", "US_P2", "US_P3", "US_P4", "US_P5", "US_P6", "AR_E1", "AR_E2", "AR_E3", "AR_E4", "AR_E5",
             "AR_E6", "BO_E1", "BO_E2", "BO_E3", "BO_E4", "BO_E5", "BO_E6", "BR_E1", "BR_E2", "BR_E3", "BR_E4", "BR_E5", "BR_E6", "BR_E7", "CO_E1", "CO_E2", "CO_E3", "CO_E4", "CO_E5",
             "CO_E6", "CO_E7", "MX_E1", "MX_E2", "MX_E3", "MX_E4", "MX_E5", "MX_E6", "US_E1", "US_E2", "US_E3", "US_E4", "US_E5", "US_E6")

types <- c("PO")

smoothing_median_window_size <- 3
smoothing_mean_window_size <- 0
smooth_borders_option <- "zeros"

output_signal_mean_decimals <- 2
output_signal_sd_decimals <- 2

output_folder <- sprintf("%s/outputs/12_individual_serums_smoothed_data", project_folder)
output_suffix <- "_smoothed_signals.tsv"

#### AUXILIAR FUNCTIONS ####
smoothVector <- function(vector, median_window_size = 5, mean_window_size = 7, borders = "zeros") {
    # borders can be "repeat" or "zeros"
    
    if (median_window_size > 0) {
        #Fill the borders for median
        if (borders == "repeat") {
            #Fill the sides with the first and last number to have the same amount of data after the smoothing
            prefix <- rep(vector[1], floor((median_window_size - 1) / 2))
            suffix <- rep(vector[length(vector)], ceiling((median_window_size - 1) / 2))        
        } else if (borders == "zeros") {
            #Fill the sides with 0 to have the same amount of data after the smoothing
            ### This flattens the borders a bit
            prefix <- rep(0, floor((median_window_size - 1) / 2))
            suffix <- rep(0, ceiling((median_window_size - 1) / 2))
        } else {
            writeLines("WARNING: Incorrect border option.")
        }
        vector_aux <- c(prefix, vector, suffix)
        
        #Calculate the rolling median
        smoothed_vector <- round(rollmedian(vector_aux, median_window_size), 3)    
    } else {
        smoothed_vector <- vector #this is because the name change
    }
    
    if (mean_window_size > 0) {
        #Fill the borders for mean
        if (borders == "repeat") {
            #Fill the sides with the first and last number to have the same amount of data after the smoothing
            prefix <- rep(smoothed_vector[1], floor((mean_window_size - 1) / 2))
            suffix <- rep(smoothed_vector[length(smoothed_vector)], ceiling((mean_window_size - 1) / 2))        
        } else if (borders == "zeros") {
            #Fill the sides with 0 to have the same amount of data after the smoothing
            ### This flattens the borders a bit
            prefix <- rep(0, floor((mean_window_size - 1) / 2))
            suffix <- rep(0, ceiling((mean_window_size - 1) / 2))
        } else {
            writeLines("WARNING: Incorrect border option.")
        }
        vector_aux <- c(prefix, smoothed_vector, suffix)
        
        #Calculate the rolling mean
        smoothed_vector <- round(rollmean(vector_aux, mean_window_size), 3)        
    }
    
    smoothed_vector
}

#### SMOOTH DATA ####
design_data <- fread(design_data_file, header = T, sep = "\t", na.strings = NULL)
design_data <- design_data[group %in% design_groups]

ids_in_design <- design_data$array_express_id

#Get data
for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
        # type_for <- types[1]
        normalized_data_file <- sprintf("%s/%s_%s_processed.tsv", normalized_data_folder, source_for, type_for)  
        normalized_data <- fread(normalized_data_file, header = T, sep = "\t", na.strings = "")
        
        setnames(normalized_data, "Reporter Name", "array_express_id")
        setnames(normalized_data, "onlyRegionsData.replica1", "r1")
        setnames(normalized_data, "onlyRegionsData.replica2", "r2")
        
        normalized_data <- normalized_data[, .(array_express_id, r1, r2)]
        
        #Keep only sequences in design
        normalized_data <- normalized_data[array_express_id %in% ids_in_design]
        
        #Combine the replicas in the same column
        normalized_data <- rbindlist(list(normalized_data[, .(array_express_id, source = source_for, type = type_for, replica = 1, signal = r1)],
                                          normalized_data[, .(array_express_id, source = source_for, type = type_for, replica = 2, signal = r2)]))
        
        #Add the design data
        normalized_data <- merge(normalized_data,
                                 design_data[, .(protein, region, start, array_express_id)],
                                 by = "array_express_id",
                                 allow.cartesian = T)
        
        normalized_data <- normalized_data[order(source, type, replica, protein, region, start)]
        
        #Smooth the signal
        smoothed_normalized_data_aux <- normalized_data[, .(smoothed_signal = smoothVector(vector = signal,
                                                                                           median_window_size = smoothing_median_window_size, 
                                                                                           mean_window_size = smoothing_mean_window_size, 
                                                                                           borders = smooth_borders_option)),
                                                        by = .(source, type, replica, protein, region)]
        normalized_data$smoothed_signal <- smoothed_normalized_data_aux$smoothed_signal
        
        #Combine both replicas
        normalized_data <- normalized_data[, .(mean_smoothed_signal = round(mean(smoothed_signal), output_signal_mean_decimals),
                                               sd_smoothed_signal = round(sd(smoothed_signal), output_signal_sd_decimals)),
                                           by = .(source, type, protein, region, start)]
        
        #Add sequence data
        normalized_data <- merge(normalized_data,
                                 design_data[, .(protein, start, sequence, truncated)],
                                 by = c("protein", "start"))
        
        #Sort the columns
        normalized_data <- normalized_data[order(source, type, protein, region, start)]
        setcolorder(normalized_data, c("source", "type", "protein", "region", "start",
                                       "mean_smoothed_signal", "sd_smoothed_signal",
                                       "sequence", "truncated"))
        
        #Write data
        output_file <- sprintf("%s/%s_%s%s", output_folder, source_for, type_for, output_suffix)
        write.table(normalized_data, file = output_file, col.names = T, row.names = F, sep = "\t", quote = T) 
    }
}
