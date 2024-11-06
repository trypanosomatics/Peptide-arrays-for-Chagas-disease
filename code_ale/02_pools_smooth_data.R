#### PREPARATION STEPS ####
## You can run this code as it is to process a small subset of proteins, or you can follow these next steps to analyze the entire dataset.

## 1. Make sure you have run all previous codes
## 2. In chagastope_data/inputs/01_pools_array_design place the "Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv" file (links in the paper)
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

design_data_file <- sprintf("%s/inputs/01_pools_array_design/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv", project_folder)
design_groups <- c("Trypanosoma_cruzi")

normalized_data_folder <- sprintf("%s/outputs/01_pools_normalized_data", project_folder)

sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")

types <- c("PO", "NE")

smooth_median_window_size <- 3
smooth_mean_window_size <- 0
smooth_borders_option <- "zeros"

output_signal_mean_decimals <- 2
output_signal_sd_decimals <- 2

output_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)
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
        normalized_data <- fread(normalized_data_file, header = T, sep = "\t", na.strings = NULL)
        
        setnames(normalized_data, "Reporter Name", "array_express_id")
        setnames(normalized_data, sprintf("%s_%s_normalizedData_r1", source_for, type_for), "r1")
        setnames(normalized_data, sprintf("%s_%s_normalizedData_r2", source_for, type_for), "r2")
        
        #Keep only sequences in design
        normalized_data <- normalized_data[array_express_id %in% ids_in_design]
        
        #Combine the replicas in the same column
        normalized_data <- rbindlist(list(normalized_data[, .(array_express_id, source = source_for, type = type_for, replica = 1, signal = r1)],
                                          normalized_data[, .(array_express_id, source = source_for, type = type_for, replica = 2, signal = r2)]))
        
        #Add the design data
        normalized_data <- merge(normalized_data,
                                 design_data[, .(protein, start, array_express_id)],
                                 by = "array_express_id",
                                 allow.cartesian = T)
        
        normalized_data <- normalized_data[order(source, type, replica, protein, start)]
        
        #Smooth the signal
        smoothed_normalized_data_aux <- normalized_data[, .(smoothed_signal = smoothVector(vector = signal,
                                                                                           median_window_size = smooth_median_window_size, 
                                                                                           mean_window_size = smooth_mean_window_size, 
                                                                                           borders = smooth_borders_option)),
                                                        by = .(source, type, replica, protein)]
        normalized_data$smoothed_signal <- smoothed_normalized_data_aux$smoothed_signal
        
        #Combine both replicas
        normalized_data <- normalized_data[, .(mean_smoothed_signal = round(mean(smoothed_signal), output_signal_mean_decimals),
                                               sd_smoothed_signal = round(sd(smoothed_signal), output_signal_sd_decimals)),
                                           by = .(source, type, protein, start)]
        
        #Add sequence data
        normalized_data <- merge(normalized_data,
                                 design_data[, .(protein, start, sequence, truncated)],
                                 by = c("protein", "start"))
        
        #Sort the columns
        normalized_data <- normalized_data[order(source, type, protein, start)]
        setcolorder(normalized_data, c("source", "type", "protein", "start",
                                       "mean_smoothed_signal", "sd_smoothed_signal",
                                       "sequence", "truncated"))
        
        #Write data
        output_file <- sprintf("%s/%s_%s%s", output_folder, source_for, type_for, output_suffix)
        write.table(normalized_data, file = output_file, col.names = T, row.names = F, sep = "\t", quote = T) 
    }
}
