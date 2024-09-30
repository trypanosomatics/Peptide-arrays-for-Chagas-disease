#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}
if (!require(preprocessCore, quietly = TRUE)) {
  writeLines("Installing library 'preprocessCore' for R")
  install.packages("preprocessCore", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(preprocessCore) #to modify pheatmap
}

#### AUXILIARY FUNCTIONS ####
calculateMode <- function(x, decimals = 0) {
  x <- round(x, decimals)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#### MAIN FUNCTION ####
#### NORMALIZE DATA ####
normalize_and_process_data <- function(raw_data_folder, output_folder, sources) {
  #Get data
  output_initialized <- FALSE
  sources_to_normalize_as_negatives <- c("LE")
  sources_to_exclude_from_global_statistics <- c("LE")
  types <- c("PO", "NE")
  
  for (source_for in sources) {
    # source_for <- sources[1]
    for (type_for in types) {
      # type_for <- types[1]
      raw_data_file <- sprintf("%s/%s_%s_raw.tsv", raw_data_folder, source_for, type_for)
      raw_data <- fread(raw_data_file, header = TRUE, sep = "\t", na.strings = NULL)
      
      #Check if I need to initialize the tables to normalize
      if (!output_initialized) {
        full_raw_PO_data <- raw_data[, .(`Reporter Name`)]
        full_raw_NE_data <- raw_data[, .(`Reporter Name`)]
        output_initialized <- TRUE
      }
      
      #Add the data to the corresponding table
      if ((type_for == "NE") | (source_for %in% sources_to_normalize_as_negatives)) {
        full_raw_NE_data <- merge(full_raw_NE_data, raw_data, by = "Reporter Name")
      } else {
        full_raw_PO_data <- merge(full_raw_PO_data, raw_data, by = "Reporter Name")
      }
    }
  }
  
  #Normalize the positive data
  matrix_PO_aux <- normalize.quantiles(as.matrix(full_raw_PO_data[, -c("Reporter Name")]))
  matrix_PO_aux <- round(matrix_PO_aux, digits = 2)
  normalized_PO_data <- as.data.table(matrix_PO_aux)
  colnames(normalized_PO_data) <- colnames(full_raw_PO_data[, -c("Reporter Name")])
  normalized_PO_data$`Reporter Name` <- full_raw_PO_data$`Reporter Name`
  rm(full_raw_PO_data)
  rm(matrix_PO_aux)
  gc()
  
  #Normalize the negative data
  matrix_NE_aux <- normalize.quantiles(as.matrix(full_raw_NE_data[, -c("Reporter Name")]))
  matrix_NE_aux <- round(matrix_NE_aux, digits = 2)
  normalized_NE_data <- as.data.table(matrix_NE_aux)
  colnames(normalized_NE_data) <- colnames(full_raw_NE_data[, -c("Reporter Name")])
  normalized_NE_data$`Reporter Name` <- full_raw_NE_data$`Reporter Name`
  rm(full_raw_NE_data)
  rm(matrix_NE_aux)
  gc()
  
  #Combine data
  full_normalized_data <- merge(normalized_PO_data, normalized_NE_data, by = "Reporter Name")
  rm(normalized_PO_data)
  rm(normalized_NE_data)
  gc()
  
  colnames(full_normalized_data) <- gsub("rawData", "normalizedData", colnames(full_normalized_data))
  
  #### CALCULATE GLOBAL STATISTICS ####
  cols_to_select <- setdiff(colnames(full_normalized_data), "Reporter Name")
  if (length(sources_to_exclude_from_global_statistics) > 0) {
    for (source_for in sources_to_exclude_from_global_statistics) {
      cols_to_select <- cols_to_select[!grepl(sprintf("%s_", source_for), cols_to_select)]
    }
  }
  
  global_signals <- as.vector(as.matrix(full_normalized_data[, cols_to_select, with = FALSE]))
  mode_aux <- calculateMode(global_signals, decimals = 2)
  sd_aux <- round(sd(global_signals), 2)
  
  global_statistics <- data.table(mode = mode_aux, sd = sd_aux)
  
  #### SAVE DATA ####
  global_statistics_output_file <- sprintf("%s/global_statistics.tsv", output_folder)
  write.table(global_statistics, file = global_statistics_output_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = TRUE)
  
  for (source_for in sources) {
    for (type_for in types) {
      cols_to_select <- colnames(full_normalized_data)
      cols_to_select <- cols_to_select[grepl(sprintf("%s_%s", source_for, type_for), cols_to_select)]
      cols_to_select <- c("Reporter Name", cols_to_select)
      
      normalized_data_for <- full_normalized_data[, cols_to_select, with = FALSE]
      normalized_data_output_file <- sprintf("%s/%s_%s_processed.tsv", output_folder, source_for, type_for)
      write.table(normalized_data_for, file = normalized_data_output_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = TRUE)
    }
  }
}