#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}

if (!require(preprocessCore, quietly = TRUE)) {
  writeLines("Installing library 'preprocessCore' for R")
  install.packages("preprocessCore", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(preprocessCore) #quantile.normalization
}

#### AUXILIAR FUNCTIONS ####
calculateMode <- function(x, decimals = 0) {
  x <- round(x, decimals)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#### MAIN FUNCTION ####
normalize_serums <- function(main_folder, testing, sources) {
  
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
}