#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}
if (!require(zoo, quietly = TRUE)) {
  writeLines("Installing library 'zoo' for R")
  install.packages("zoo", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(zoo)
}

#### AUXILIARY FUNCTIONS ####
smoothVector <- function(vector, median_window_size, mean_window_size, borders) {
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

#### MAIN FUNCTION ####
#### SMOOTH DATA ####

smooth_data <- function(design_data_file, design_group, normalized_data_folder, normalized_data_suffix, sources, smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option, output_type_order, output_signal_mean_decimals, output_signal_sd_decimals, output_suffix) {
  # Get data
  design_data <- fread(design_data_file, header = TRUE, sep = "\t", na.strings = NULL)
  design_data <- design_data[group == design_group]
  types <- c("PO", "NE")
  
  # Create an extra design
  design_data_aux <- design_data[, .(protein, start, sequence = final_sequence, truncated, array_express_id)] #remove extra columns
  sequences_in_design <- unique(design_data_aux$array_express_id)
  
  for (source_for in sources) {
    #source_for = sources[1]
    print(source_for)
    for (type_for in types) {
      normalized_data_file <- sprintf("%s%s_%s%s", normalized_data_folder, source_for, type_for, normalized_data_suffix)
      normalized_data <- fread(normalized_data_file, header = TRUE, sep = "\t", na.strings = NULL)
      
      normalized_data <- normalized_data[`Reporter Name` %in% sequences_in_design]

      normalized_data <- merge(normalized_data, design_data_aux, by.x = "Reporter Name", by.y = "array_express_id", allow.cartesian = TRUE)
      
      normalized_data <- normalized_data[order(protein, start)]
      
      # Smooth the signal
      smoothed_normalized_data_aux <- normalized_data[, .(
        smoothed_signal_r1 = smoothVector(get(sprintf("%s_%s_%s", source_for, type_for, "normalizedData_r1")), smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option),
        smoothed_signal_r2 = smoothVector(get(sprintf("%s_%s_%s", source_for, type_for, "normalizedData_r2")), smoothing_median_window_size, smoothing_mean_window_size, smooth_borders_option)
      ), by = .(protein)]
      
      normalized_data$smoothed_signal_r1 <- smoothed_normalized_data_aux$smoothed_signal_r1
      normalized_data$smoothed_signal_r2 <- smoothed_normalized_data_aux$smoothed_signal_r2
      
      # Combine both replicas
      normalized_data <- normalized_data[, .(
        mean_smoothed_signal = round(mean(c(smoothed_signal_r1, smoothed_signal_r2), na.rm = TRUE), output_signal_mean_decimals),
        sd_smoothed_signal = round(sd(c(smoothed_signal_r1, smoothed_signal_r2), na.rm = TRUE), output_signal_sd_decimals)
      ), by = .(protein, start)]
      
      # Add sequence data
      normalized_data <- merge(normalized_data, design_data[, .(protein, start, sequence, truncated)], by = c("protein", "start"))
      
      # Sort columns
      normalized_data$type <- type_for
      normalized_data$source <- source_for
      normalized_data <- normalized_data[order(protein, start)]
      setcolorder(normalized_data, c("source", "type", "protein", "start",
                                     "mean_smoothed_signal", "sd_smoothed_signal",
                                     "sequence", "truncated"))
      
      # Save data
      output_file <- sprintf("%s/%s_%s_smoothed.tsv", output_folder, source_for, type_for)
      write.table(normalized_data, file = output_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = TRUE)
    }
  }
}
