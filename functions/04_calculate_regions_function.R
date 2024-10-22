#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}

#### AUXILIAR FUNCTIONS ####
calculateRegions <- function(peaks_data,
                             fasta_data,
                             peaks_types,
                             peaks_max_other_type_signal_ratio,
                             peaks_peptide_length,
                             peaks_peptide_offset,
                             peaks_min_peptide_length,
                             peaks_border_peptide_amount,
                             peaks_tag) {
  ### Filter the peaks
  peaks_data <- peaks_data[!(grepl(",", source))]
  peaks_data <- peaks_data[type %in% peaks_types]
  
  if (peaks_max_other_type_signal_ratio != -1) {
    peaks_data <- peaks_data[other_type_best_ratio_original_peak_ratio <= peaks_max_other_type_signal_ratio]   
  }
  
  peaks_data <- peaks_data[, .(protein, peak_start, peak_last_start)]
  peaks_data <- unique(peaks_data)
  
  sub_fasta_data <- fasta_data[protein %in% unique(peaks_data$protein)]
  
  peptide_data <- splitFastaIntoPeptides(fasta_data = sub_fasta_data,
                                         peptide_length = peaks_peptide_length, 
                                         peptide_offset = peaks_peptide_offset, 
                                         min_peptide_length = peaks_min_peptide_length)
  
  peaks_border_start_amount <- peaks_border_peptide_amount * peaks_peptide_offset
  
  peptide_data <- merge(peptide_data,
                        peaks_data,
                        by = "protein",
                        all.x = T,
                        allow.cartesian = T)
  
  peptide_data <- peptide_data[(start >= (peak_start - peaks_border_start_amount)) &
                                 (start <= (peak_last_start + peaks_border_start_amount))]
  
  peptide_data$is_border <- 1
  peptide_data[(start >= peak_start) & (start <= peak_last_start), is_border := 0]
  
  peptide_data <- peptide_data[, .(is_border = as.numeric(all(as.logical(is_border)))), by = .(protein, start, peptide)]
  peptide_data <- peptide_data[, .(protein, start, peptide, is_border)]
  
  peptide_data$origin <- peaks_tag
  
  peptide_data
}

splitFastaIntoPeptides <- function(fasta_data, peptide_length, peptide_offset, min_peptide_length) {
  fasta_data[, seq_len := nchar(sequence)]
  fasta_data[, pep_amount := ceiling((seq_len - peptide_length + 1) / peptide_offset)]
  
  shorter_proteins <- fasta_data[(seq_len >= min_peptide_length) & (seq_len < peptide_length)]
  fasta_data <- fasta_data[seq_len >= peptide_length]
  
  if (fasta_data[, .N] > 0) {
    peptide_data <- data.table(protein = rep(fasta_data$protein, fasta_data$pep_amount),
                               sequence = rep(fasta_data$sequence, fasta_data$pep_amount))
    
    peptide_data[, start := .SD[,.I], by = protein]
    peptide_data[, start := ((start - 1) * peptide_offset) + 1]
    peptide_data[, end := start + peptide_length - 1]
    peptide_data[, peptide := substring(sequence, start, end)]
    
    peptide_data <- peptide_data[, -c("sequence")]
  } else {
    peptide_data <- data.table(protein = character(), 
                               start = numeric(), 
                               end = numeric(), 
                               peptide = character())
  }
  
  if (shorter_proteins[, .N] > 0) {
    shorter_proteins$start <- 1
    shorter_proteins[, end := start + nchar(sequence) - 1]
    
    shorter_proteins <- shorter_proteins[, .(protein, start, end, peptide = sequence)]
    
    peptide_data <- rbindlist(list(peptide_data, shorter_proteins))
  }
  
  peptide_data
}

#### MAIN FUNCTION ####
processAntigenicRegions <- function(project_folder, extra_files_folder, peaks_tag,
                                    chronic_peaks_types, chronic_peaks_max_other_type_signal_ratio, 
                                    chronic_peaks_peptide_length, chronic_peaks_peptide_offset, 
                                    chronic_peaks_min_peptide_length, chronic_peaks_border_peptide_amount, 
                                    detailed_region_output_file, summary_region_output_file) {
  
  chronic_peaks_data_file <- sprintf("%s/outputs/03_pools_antigenic_peaks/pools_peaks_cutoff4SD_2pep.tsv", project_folder)
  sylvio_fasta_data_file <- sprintf("%s/sylvio-x10-tritrypdb-5.Dec.2016-no-pseudo_SIMPLIFIED_TABBED.tsv", extra_files_folder)
  brener_fasta_data_file <- sprintf("%s/tcruzi-clbrener-tritrypdb-5.Dec.2016-no-pseudo_SIMPLIFIED_TABBED.tsv", extra_files_folder)
  
  chronic_peaks_data <- fread(chronic_peaks_data_file, header = TRUE, sep = "\t", na.strings = NULL)
  sylvio_fasta_data <- fread(sylvio_fasta_data_file, header = TRUE, sep = "\t", na.strings = NULL)
  brener_fasta_data <- fread(brener_fasta_data_file, header = TRUE, sep = "\t", na.strings = NULL)
  
  chronic_all_fasta <- rbindlist(list(sylvio_fasta_data, brener_fasta_data))
  
  region_data <- calculateRegions(peaks_data = chronic_peaks_data,
                                  fasta_data = chronic_all_fasta,
                                  peaks_types = chronic_peaks_types,
                                  peaks_max_other_type_signal_ratio = chronic_peaks_max_other_type_signal_ratio,
                                  peaks_peptide_length = chronic_peaks_peptide_length,
                                  peaks_peptide_offset = chronic_peaks_peptide_offset,
                                  peaks_min_peptide_length = chronic_peaks_min_peptide_length,
                                  peaks_border_peptide_amount = chronic_peaks_border_peptide_amount,
                                  peaks_tag = chronic_peaks_tag)
  
  ## Calculate the region numbers
  region_data <- region_data[order(protein, start)]
  region_data[, last_protein := c("", protein[1:.N-1])]
  region_data[, last_start := c(-1, start[1:.N-1])]
  region_data$same_region <- 1
  region_data[(protein != last_protein) | (start > last_start + chronic_peaks_peptide_offset), same_region := 0]
  region_starts <- which(region_data$same_region == 0)
  region_aux <- c()
  for (i in 1:length(region_starts)) {
    if (i == 1) {
      start_aux <- region_starts[i]
    } else {
      end_aux <- region_starts[i] - 1
      region_aux <- c(region_aux, rep(i - 1, times = end_aux - start_aux + 1))
      
      start_aux <- region_starts[i]
    }
  }
  end_aux <- region_data[, .N]
  region_aux <- c(region_aux, rep(i, times = end_aux - start_aux + 1))
  region_data$region <- region_aux
  
  region_data <- region_data[, .(protein, start, region, is_border, peptide, origin)]
  
  ## Calculate the region sequence
  region_data[, sequence_aux := substring(peptide, 16, 16)]
  region_data[, relative_pos_aux := .SD[, .I], by = region]
  region_data[relative_pos_aux == 1, sequence_aux := peptide]
  
  region_summary <- region_data[, .(protein = .SD[1]$protein,
                                    start = min(start),
                                    last_start = max(start),
                                    sequence = paste(sequence_aux, collapse = "")), by = region]
  
  region_summary[, aa_amount := last_start - start + chronic_peaks_peptide_length]
  
  ## Calculate the region sequence without OUTSIDE borders
  #Calculate the min and max start for the borders in each protein
  border_info_aux <- region_data[is_border == 0, .(no_border_start = min(start),
                                                   no_border_last_start = max(start)), by = region]
  
  #Add the data
  no_border_region_data <- merge(region_data,
                                 border_info_aux,
                                 by = "region")
  
  #Keep only data without the outside borders
  no_border_region_data <- no_border_region_data[(start >= no_border_start) & (start <= no_border_last_start)]
  
  #Calculate the rest
  no_border_region_data[, sequence_aux := substring(peptide, 16, 16)]
  no_border_region_data[, relative_pos_aux := .SD[, .I], by = region]
  no_border_region_data[relative_pos_aux == 1, sequence_aux := peptide]
  
  no_border_region_summary <- no_border_region_data[, .(protein = .SD[1]$protein,
                                                        no_border_start = min(start),
                                                        no_border_last_start = max(start),
                                                        no_border_sequence = paste(sequence_aux, collapse = "")), by = region]
  
  no_border_region_summary[, no_border_aa_amount := no_border_last_start - no_border_start + chronic_peaks_peptide_length]
  
  #Add to global summary
  region_summary <- merge(region_summary,
                          no_border_region_summary,
                          by = c("region", "protein"))
  
  #Remove extra cols
  region_data <- region_data[, .(protein, start, region, is_border, peptide)]
  
  #Save data
  write.table(region_data, file = detailed_region_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
  write.table(region_summary, file = summary_region_output_file, col.names = T, row.names = F, sep = "\t", quote = T)
}