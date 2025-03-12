#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}

determine_peaks <- function(main_folder, testing, sources, types, project_folder, output_folder, min_num_of_peptides_in_peak, sd_multiplier_for_cutoff, profile_data_suffix) {
  
  global_statistics_file <- sprintf("%s/outputs/01_pools_normalized_data/global_statistics.tsv", project_folder)
  global_statistics <- fread(global_statistics_file, header = TRUE, sep = "\t", na.strings = NULL)
  cutoff <- global_statistics$mode + sd_multiplier_for_cutoff * global_statistics$sd
  
  #### AUXILIAR FUNCTIONS ####
  calculateIntervals <- function(vector, threshold_for_ignoring_negative_interval = 0) {
    output <- c()
    if (length(vector) >= 2) {
      #This function receives a list of 1/0 and returns the intervals where the 1s are found
      vector_aux <- vector[-length(vector)]
      vector_next_pos <- vector[-1]
      
      #In comparison the 1s are the positions of the original vector where an interval ended
      comparison <- as.numeric(vector_aux != vector_next_pos)
      
      #This gives me the positions of each 1 in comparison (where the intervals ended)
      changing_positions <- which(comparison %in% 1)
      
      #Check for every interval if it's of 0s or of 1s
      start <- 1
      if (length(changing_positions) > 0) {
        for (i in 1:length(changing_positions)) {
          #i <- 1
          position_for <- changing_positions[i]
          
          if (vector[position_for] == 1) {
            end <- position_for
            
            output <- c(output, start, end)
          }
          
          start <- position_for + 1
        }
      }
      #Add the last case by hand
      if (vector[length(vector)] == 1) {
        end <- length(vector)
        output <- c(output, start, end)
      }
      
      #Now, merge together peaks that are closer than the threshold
      if (threshold_for_ignoring_negative_interval > 0) {
        if (length(output) >= 4) {
          parsed_output <- c()
          
          first_index <- 1
          second_index <- 3
          first_start <- output[first_index]
          first_end <- output[first_index + 1]
          second_start <- output[second_index]
          second_end <- output[second_index + 1]
          finished <- 0
          while(!finished) {
            if (((second_start - first_end) - 1) <= threshold_for_ignoring_negative_interval) {
              first_end <- second_end
              
              second_index <- second_index + 2
              if (second_index < length(output)) {
                second_start <- output[second_index]
                second_end <- output[second_index + 1]                            
              } else {
                parsed_output <- c(parsed_output, first_start, first_end)
                finished <- 1
              }
            } else {
              parsed_output <- c(parsed_output, first_start, first_end)
              
              first_index <- second_index
              second_index <- first_index + 2
              first_start <- output[first_index]
              first_end <- output[first_index + 1]
              if (second_index < length(output)) {
                second_start <- output[second_index]
                second_end <- output[second_index + 1]                            
              } else {
                parsed_output <- c(parsed_output, first_start, first_end)
                finished <- 1
              }
            }
          }
          
          output <- parsed_output
        }
      }
    }
    
    output
  }
  
  #### CALCULATE PEAKS ####
  sequence_offset <- sequence_length - sequence_overlap
  
  output_initialized <- F
  for (source_for in sources) {
    # source_for <- "AR"
    
    for (type_for in types) {
      # type_for <- types[1]
      
      profile_data_file <- sprintf("%s/%s_%s_%ssmoothed.tsv", profile_data_folder, source_for, type_for, profile_data_suffix)
      profile_data <- fread(profile_data_file, header = TRUE, sep = "\t", na.strings = NULL)
      
      if (output_initialized == F) {
        all_profile_data <- profile_data
        output_initialized <- T
      } else {
        all_profile_data <- rbindlist(list(all_profile_data, profile_data))
      }
    }
  }
  rm(profile_data)
  gc()
  ## Find the peaks
  #Set a flag about passing or not the cutoff
  all_profile_data$above_cutoff <- 0
  all_profile_data[mean_smoothed_signal >= cutoff, above_cutoff := 1]
  #Find the INDEX of the starting and ending position of peaks (named fake start)
  peak_data <- all_profile_data[, calculateIntervals(above_cutoff), by = .(source, type, protein)]
  setnames(peak_data, "V1", "fake_peak_start") #the real starts will be calculated later
  #Divide the starts and ends and add them to different columns
  peak_data$is_start <- rep(c(1, 0), peak_data[, .N] / 2)
  peak_data_aux <- peak_data[is_start == 1]
  peak_data_aux$fake_peak_last_start <- peak_data[is_start == 0]$fake_peak_start
  peak_data <- peak_data_aux[, -c("is_start")]
  peak_data[, peptide_amount := fake_peak_last_start - fake_peak_start + 1]
  #Sort the data
  peak_data <- peak_data[order(protein, source, fake_peak_start)]
  
  ## Filter peaks by width
  peak_data <- peak_data[peptide_amount >= min_num_of_peptides_in_peak]
  
  ## Prepare the data to calculate the real starts later on
  real_starts <- unique(all_profile_data[, .(protein, start)])
  real_starts[, fake_peak_start := c(1:.N), by = protein]
  setnames(real_starts, "start", "peak_start")
  real_last_starts <- real_starts
  real_last_starts[1] <- real_last_starts[1] #unlink
  setnames(real_last_starts, "peak_start", "peak_last_start")
  setnames(real_last_starts, "fake_peak_start", "fake_peak_last_start")
  
  ## In each protein, find matching peaks and get data for the other peaks
  output_initialized <- 0
  unique_proteins <- unique(peak_data$protein)
  unique_types <- unique(all_profile_data$type)
  progress_i <- 1
  progress_total <- length(unique_proteins)
  for (protein_for in unique_proteins) {
    # protein_for <- unique_proteins[4]
    writeLines(sprintf("%s (%s/%s)", protein_for, progress_i, progress_total))
    progress_i <- progress_i + 1
    
    sub_peak_data <- peak_data[protein == protein_for]
    sub_all_profile_data <- all_profile_data[protein == protein_for]
    
    for (type_for in unique_types) {
      # type_for <- "PO"
      # type_for <- "NE"
      peaks_to_parse <- sub_peak_data[type == type_for]
      
      if (peaks_to_parse[, .N] > 0) {
        overlapping_peaks_to_parse <- peaks_to_parse[0,] #save format for later
        
        ### SINGLE SOURCE PEAKS
        ## Parse all peaks, fixing the starts and fetching signal information
        #Get the REAL starts and ends for the peptides
        #1, 5, 9...
        peaks_to_parse <- merge(peaks_to_parse, real_starts, by = c("protein", "fake_peak_start"), all.x = T)
        peaks_to_parse <- merge(peaks_to_parse, real_last_starts, by = c("protein", "fake_peak_last_start"), all.x = T)
        peaks_to_parse[, peak_last_end := peak_last_start + sequence_length - 1]
        
        #Sort to put similar sources together (to reduce the amount of subsets later on)
        peaks_to_parse <- peaks_to_parse[order(source, fake_peak_start, fake_peak_last_start)] 
        
        #For each peak, fetch the data
        last_source <- ""
        for (peak_to_parse_i in 1:peaks_to_parse[, .N]) {
          # peak_to_parse_i <- 3
          # peak_to_parse_i <- 1
          peak_to_parse_for <- peaks_to_parse[peak_to_parse_i]
          
          #See if you need to update the data because a new source combination
          if (peak_to_parse_for$source != last_source) {
            peak_to_parse_sources <- unlist(strsplit(peak_to_parse_for$source, ", ")) #this isn't needed here, but this way the code matches the one below
            sub_sub_all_profile_data <- sub_all_profile_data[source %in% peak_to_parse_sources]
            
            last_source <- peak_to_parse_for$source
          }
          
          #Extract the data for this start range
          peak_profile_data <- sub_sub_all_profile_data[(start >= peak_to_parse_for$peak_start) & (start <= peak_to_parse_for$peak_last_start)]
          
          #Fetch the information for this peak
          best_peptide <- peak_profile_data[type == type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
          best_peak_by_source <- peak_profile_data[type == type_for][, .(peak_signal = sum(mean_smoothed_signal)), by = source]
          peak_sequence_data <- peak_profile_data[type == type_for, .(protein, start, sequence)]
          peak_sequence_data[, sequence_aux := substring(sequence, sequence_length - sequence_offset + 1, sequence_length)]
          peak_sequence <- paste(c(peak_sequence_data[1]$sequence, peak_sequence_data[-1]$sequence_aux), collapse = "")
          
          #Fetch information for the other type
          other_type_best_peptide <- peak_profile_data[type != type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
          other_type_best_peak_by_source <- peak_profile_data[type != type_for][, .(peak_signal = sum(mean_smoothed_signal)), by = source]
          other_type_best_peak_by_source <- merge(other_type_best_peak_by_source,
                                                  best_peak_by_source[, .(source, main_type_peak_signal = peak_signal)],
                                                  by = "source")
          other_type_best_peak_by_source[, other_type_peak_signal_proportion := round(peak_signal / main_type_peak_signal, other_type_proportion_decimals)]
          
          peak_to_parse_for$best_peptide <- best_peptide$sequence
          peak_to_parse_for$best_peptide_source <- best_peptide$source
          peak_to_parse_for$best_peptide_start <- best_peptide$start
          peak_to_parse_for$best_peptide_signal <- best_peptide$mean_smoothed_signal
          peak_to_parse_for$best_peak_source <- best_peak_by_source$source
          peak_to_parse_for$best_peak_signal <- best_peak_by_source$peak_signal
          peak_to_parse_for$combined_peak_best_peptide <- ""
          peak_to_parse_for$combined_peak_best_peptide_start <- -1
          peak_to_parse_for$combined_peak_best_peptide_signal <- -1
          peak_to_parse_for$combined_peak_peak_signal <- -1
          peak_to_parse_for$peak_sequence <- peak_sequence
          
          peak_to_parse_for$other_type_best_peptide <- other_type_best_peptide$sequence
          peak_to_parse_for$other_type_best_peptide_source <- other_type_best_peptide$source
          peak_to_parse_for$other_type_best_peptide_start <- other_type_best_peptide$start
          peak_to_parse_for$other_type_best_peptide_signal <- other_type_best_peptide$mean_smoothed_signal
          peak_to_parse_for$other_type_best_ratio_original_peak_source <- other_type_best_peak_by_source$source
          peak_to_parse_for$other_type_best_ratio_original_peak_ratio <- other_type_best_peak_by_source$other_type_peak_signal_proportion
          peak_to_parse_for$other_type_combined_peak_best_peptide <- ""
          peak_to_parse_for$other_type_combined_peak_best_peptide_start <- -1
          peak_to_parse_for$other_type_combined_peak_best_peptide_signal <- -1
          peak_to_parse_for$other_type_combined_peak_peak_signal <- -1
          
          if (peak_to_parse_i == 1) {
            peaks_to_parse_output <- peak_to_parse_for
          } else {
            peaks_to_parse_output <- rbindlist(list(peaks_to_parse_output, peak_to_parse_for))
          }
        }
        
        if (output_initialized == 0) {
          output_peak_data <- peaks_to_parse_output
          output_initialized <- 1
        } else {
          output_peak_data <- rbindlist(list(output_peak_data, peaks_to_parse_output))
        }
        
        ### MULTIPLE SOURCES PEAKS
        ## Create a boolean list for each source for where it has peaks
        peak_in_position <- list()
        max_fake_last_start <- max(peaks_to_parse$fake_peak_last_start)
        unique_sources <- unique(peaks_to_parse$source)
        if (length(unique_sources) > 1) {
          for (source_for in unique_sources) {
            # source_for <- unique_sources[1]
            peak_in_position_aux <- rep(0, max_fake_last_start)
            
            sub_peaks_to_parse <- peaks_to_parse[source == source_for]
            for (peak_i in 1:sub_peaks_to_parse[,.N]) {
              # peak_i <- 1
              peak_for <- sub_peaks_to_parse[peak_i]
              peak_in_position_aux[peak_for$fake_peak_start:peak_for$fake_peak_last_start] <- 1
            }
            
            peak_in_position[[source_for]] <- peak_in_position_aux
          }
          
          #Find all posible combinations of sources of at least two sources
          if (length(unique_sources) > 1) {
            source_combinations <- list()
            for (i in 2:length(unique_sources)) {
              source_combinations <- c(source_combinations,
                                       combn(unique_sources, i, simplify = F))
            }    
            
            
          }
          
          #For each possible combination find if there is overlapping peaks
          for (i in 1:length(source_combinations)) {
            # i <- 1
            combination_for <- source_combinations[[i]]
            
            overlapping_positions <- rep(1, max_fake_last_start)
            for (source_for in combination_for) {
              overlapping_positions <- overlapping_positions * peak_in_position[[source_for]]
            }
            
            #Check if there are at least one overlapping peak
            if (sum(overlapping_positions) > 0) {
              #Find the INDEX of the starting and ending position of overlapping peaks (named fake start)
              overlapping_peak_data <- data.table(fake_peak_start = calculateIntervals(overlapping_positions))
              #Divide the starts and ends and add them to different columns
              overlapping_peak_data$is_start <- rep(c(1, 0), overlapping_peak_data[, .N] / 2)
              peak_data_aux <- overlapping_peak_data[is_start == 1]
              peak_data_aux$fake_peak_last_start <- overlapping_peak_data[is_start == 0]$fake_peak_start
              overlapping_peak_data <- peak_data_aux[, -c("is_start")]
              overlapping_peak_data[, peptide_amount := fake_peak_last_start - fake_peak_start + 1]
              
              #Add the combination data
              overlapping_peak_data$source <- paste(sort(combination_for), collapse = ", ")
              overlapping_peak_data$type <- type_for
              overlapping_peak_data$protein <- protein_for
              
              #Sort the columns to match the peak data
              setcolorder(overlapping_peak_data, colnames(overlapping_peaks_to_parse))
              
              #Filter peaks by width
              overlapping_peak_data <- overlapping_peak_data[peptide_amount >= min_num_of_peptides_in_peak]
              
              #Add the overlapping peak to the peak data
              if (overlapping_peak_data[, .N] > 0) {
                overlapping_peaks_to_parse <- rbindlist(list(overlapping_peaks_to_parse,
                                                             overlapping_peak_data))
              }
            }
          }
          
          if (overlapping_peaks_to_parse[, .N] > 0) {
            #### Parse all overlapping peaks, fixing the starts and fetching signal information
            #Get the REAL starts and ends for the peptides
            #1, 5, 9...
            overlapping_peaks_to_parse <- merge(overlapping_peaks_to_parse, real_starts, by = c("protein", "fake_peak_start"), all.x = T)
            overlapping_peaks_to_parse <- merge(overlapping_peaks_to_parse, real_last_starts, by = c("protein", "fake_peak_last_start"), all.x = T)
            overlapping_peaks_to_parse[, peak_last_end := peak_last_start + sequence_length - 1]
            
            #Sort to put similar sources together (to reduce the amount of subsets later on)
            overlapping_peaks_to_parse <- overlapping_peaks_to_parse[order(source, fake_peak_start, fake_peak_last_start)] 
            
            #For each peak, fetch the data
            last_source <- ""
            for (peak_to_parse_i in 1:overlapping_peaks_to_parse[, .N]) {
              # peak_to_parse_i <- 4
              # peak_to_parse_i <- 1
              peak_to_parse_for <- overlapping_peaks_to_parse[peak_to_parse_i]
              
              #See if you need to update the data because a new source combination
              if (peak_to_parse_for$source != last_source) {
                peak_to_parse_sources <- unlist(strsplit(peak_to_parse_for$source, ", "))
                sub_sub_all_profile_data <- sub_all_profile_data[source %in% peak_to_parse_sources]
                
                last_source <- peak_to_parse_for$source
              }
              
              #Extract the data for this start range
              peak_profile_data <- sub_sub_all_profile_data[(start >= peak_to_parse_for$peak_start) & (start <= peak_to_parse_for$peak_last_start)]
              #Combine all sources in this peak
              combined_peak_profile_data <- peak_profile_data[, .(mean_smoothed_signal = round(mean(mean_smoothed_signal), combined_mean_signal_decimals)),
                                                              by = .(type, protein, start, sequence)]
              
              #Fetch the information for this peak
              best_peptide <- peak_profile_data[type == type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
              peak_by_source_data <- peak_profile_data[type == type_for][, .(peak_signal = sum(mean_smoothed_signal)), by = source]
              best_peak_by_source <- peak_by_source_data[peak_signal == max(peak_signal)][1]
              combined_peak_best_peptide <- combined_peak_profile_data[type == type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
              peak_sequence_data <- combined_peak_profile_data[type == type_for, .(protein, start, sequence)]
              peak_sequence_data[, sequence_aux := substring(sequence, sequence_length - sequence_offset + 1, sequence_length)]
              peak_sequence <- paste(c(peak_sequence_data[1]$sequence, peak_sequence_data[-1]$sequence_aux), collapse = "")
              
              #Fetch information for the other type
              other_type_best_peptide <- peak_profile_data[type != type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
              other_type_peak_by_source_data <- peak_profile_data[type != type_for][, .(peak_signal = sum(mean_smoothed_signal)), by = source]
              sources_in_peak <- unlist(strsplit(peak_to_parse_for$source, ", "))
              min_fake_start <- peak_to_parse_for$fake_peak_start
              max_fake_start <- peak_to_parse_for$fake_peak_last_start
              original_peaks <- peaks_to_parse_output[(source %in% sources_in_peak) &
                                                        ((fake_peak_start %in% c(min_fake_start:max_fake_start)) | (fake_peak_last_start %in% c(min_fake_start:max_fake_start)) | ((fake_peak_start < min_fake_start) & (fake_peak_last_start > max_fake_start)))]
              other_type_best_original_peak_by_source <- original_peaks[other_type_best_ratio_original_peak_ratio == max(other_type_best_ratio_original_peak_ratio)][1]
              other_type_combined_peak_best_peptide <- combined_peak_profile_data[type != type_for][mean_smoothed_signal == max(mean_smoothed_signal)][1]
              other_type_best_peak_by_source[, other_type_peak_signal_proportion := round(peak_signal / main_type_peak_signal, other_type_proportion_decimals)]
              
              peak_to_parse_for$best_peptide <- best_peptide$sequence
              peak_to_parse_for$best_peptide_source <- best_peptide$source
              peak_to_parse_for$best_peptide_start <- best_peptide$start
              peak_to_parse_for$best_peptide_signal <- best_peptide$mean_smoothed_signal
              peak_to_parse_for$best_peak_source <- best_peak_by_source$source
              peak_to_parse_for$best_peak_signal <- best_peak_by_source$peak_signal
              peak_to_parse_for$combined_peak_best_peptide <- combined_peak_best_peptide$sequence
              peak_to_parse_for$combined_peak_best_peptide_start <- combined_peak_best_peptide$start
              peak_to_parse_for$combined_peak_best_peptide_signal <- combined_peak_best_peptide$mean_smoothed_signal
              peak_to_parse_for$combined_peak_peak_signal <- sum(peak_by_source_data$peak_signal)
              peak_to_parse_for$peak_sequence <- peak_sequence
              
              peak_to_parse_for$other_type_best_peptide <- other_type_best_peptide$sequence
              peak_to_parse_for$other_type_best_peptide_source <- other_type_best_peptide$source
              peak_to_parse_for$other_type_best_peptide_start <- other_type_best_peptide$start
              peak_to_parse_for$other_type_best_peptide_signal <- other_type_best_peptide$mean_smoothed_signal
              peak_to_parse_for$other_type_best_ratio_original_peak_source <- other_type_best_original_peak_by_source$source
              peak_to_parse_for$other_type_best_ratio_original_peak_ratio <- other_type_best_original_peak_by_source$other_type_best_ratio_original_peak_ratio
              peak_to_parse_for$other_type_combined_peak_best_peptide <- other_type_combined_peak_best_peptide$sequence
              peak_to_parse_for$other_type_combined_peak_best_peptide_start <- other_type_combined_peak_best_peptide$start
              peak_to_parse_for$other_type_combined_peak_best_peptide_signal <- other_type_combined_peak_best_peptide$mean_smoothed_signal
              peak_to_parse_for$other_type_combined_peak_peak_signal <- sum(other_type_peak_by_source_data$peak_signal)
              
              if (peak_to_parse_i == 1) {
                overlapping_peaks_to_parse_output <- peak_to_parse_for
              } else {
                overlapping_peaks_to_parse_output <- rbindlist(list(overlapping_peaks_to_parse_output, peak_to_parse_for))
              }
            }
            
            if (output_initialized == 0) {
              #It can't really enter here
              output_peak_data <- overlapping_peaks_to_parse_output
              output_initialized <- 1
            } else {
              output_peak_data <- rbindlist(list(output_peak_data, overlapping_peaks_to_parse_output))
            }
          }
        }
      }
    }
  }
  
  #Add the other peak proportion
  output_peak_data$other_type_combined_peak_peak_ratio <- -1
  output_peak_data[nchar(source) > 2, other_type_combined_peak_peak_ratio := round(other_type_combined_peak_peak_signal / combined_peak_peak_signal, other_type_proportion_decimals)]
  
  #### Prepare the columns for output
  output_peak_data <- output_peak_data[, -c("fake_peak_last_start", "fake_peak_start")]
  
  setcolorder(output_peak_data, c("type", "source", "protein",
                                  "peak_start", "peak_last_start", "peak_last_end",
                                  "peptide_amount"))
  
  output_peak_data <- output_peak_data[order(-type, protein, nchar(source), source)]
  output_file <- sprintf("%s/pools_peaks_cutoff%sSD_%spep.tsv", output_folder, sd_multiplier_for_cutoff, min_num_of_peptides_in_peak)

  #### Write Output
  write.table(output_peak_data, file = output_file, col.names = T, row.names = F, sep = "\t", quote = T)
  
}
