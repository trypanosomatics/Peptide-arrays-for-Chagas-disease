############################-
#### Required libraries ####
############################-
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scico)
library(patchwork)

############################-
#### Auxiliar functions ####
############################-

process_files <- function(array_config, project_folder, smoothing) {
  
  # Raw data are inputs (ArrayExpress file names), normalized and smoothed data are pipeline outputs
  folders <- list(
    v1 = c(raw = "inputs/02_pools_raw_data", normalized = "outputs/01_pools_normalized_data", smoothed = "outputs/02_pools_smoothed_data"),
    v2 = c(raw = "inputs/04_individual_serums_raw_data", normalized = "outputs/05_individual_serums_normalized_data", smoothed = "outputs/06_individual_serums_smoothed_data")
  )[[array_config$id]]
  folders <- setNames(sprintf("%s/%s", project_folder, folders), names(folders))
  
  # Individual sera are grouped by country (AR_P1 -> AR)
  group_source <- function(source) {
    if (array_config$id == "v2") sub("_.*", "", source) else source
  }
  
  #### RAW ####
  all_data_raw <- data.frame()
  files <- list.files(path = folders[["raw"]], pattern = "_raw\\.tsv$", full.names = TRUE)
  for (file in files) {
    # file <- files[1]
    data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
    
    data_long <- data %>%
      pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")
    
    # File names are <source>_<type>_raw.tsv, and individual sera sources contain "_" (AR_P1)
    file_fields <- unlist(strsplit(sub("_raw\\.tsv$", "", basename(file)), "_"))
    data_long$source <- group_source(paste(head(file_fields, -1), collapse = "_"))
    data_long$type <- tail(file_fields, 1)
    
    all_data_raw <- bind_rows(all_data_raw, data_long)
  }
  
  #### NORMALIZED ####
  all_data_normalized <- data.frame()
  files <- find_output_files(folders[["normalized"]], array_config$id, "normalized")
  for (file in files) {
    data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
    check_array(data, array_config$id, file)
    file_fields <- parse_output_file_name(file)$fields
    
    data_long <- data %>%
      select(-array) %>%
      pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")
    data_long$source <- group_source(file_fields$source)
    data_long$type <- file_fields$type
    
    all_data_normalized <- bind_rows(all_data_normalized, data_long)
  }
  
  #### SMOOTHED ####
  all_data_smoothed <- data.frame()
  files <- find_output_files(folders[["smoothed"]], array_config$id, "smoothed", smoothing = smoothing)
  if (length(files) == 0) {
    stop(sprintf("No %s smoothed files found in %s for smoothing '%s'", array_config$id, folders[["smoothed"]], smoothing))
  }
  for (file in files) {
    # file <- files[1]
    data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
    check_array(data, array_config$id, file)
    
    all_data_smoothed <- bind_rows(all_data_smoothed, select(data, -array))
  }
  all_data_smoothed$source <- group_source(all_data_smoothed$source)
  
  for (step_data in c("all_data_raw", "all_data_normalized", "all_data_smoothed")) {
    data <- get(step_data)
    data$type <- as.factor(data$type)
    data$source <- as.factor(data$source)
    assign(step_data, data)
  }
  
  return(list(all_data_raw = all_data_raw, all_data_normalized = all_data_normalized, all_data_smoothed = all_data_smoothed))
}

trans <- function(x, y_factors) {
  ifelse(x <= y_factors[1], x * y_factors[2], y_factors[1]*y_factors[2] + y_factors[3] * (x - y_factors[1]))
}

DensityJitterplot <- function(plot_data, x_column_name, y_column_name, gradient_palette = gradient_palette, y_factors = y_factors,
                              bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                              geom_point_size = 1, geom_point_shape = 16, plot_title = "", x_label = "", y_label = "",
                              axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                              fixed_scale = 0, fixed_scale_min = 0, fixed_scale_max = 0, y_lim_min = NULL, y_lim_max = NULL, legend_position = "right") {
  
  plot_data_aux <- as.data.table(plot_data)
  setnames(plot_data_aux, x_column_name, "x_column")
  setnames(plot_data_aux, y_column_name, "y_column")
  plot_data_aux$y_column <- trans(plot_data_aux$y_column, y_factors)
  
  # Density range
  min_y_value <- floor(min(plot_data_aux$y_column))
  max_y_value <- ceiling(max(plot_data_aux$y_column))
  grid_size_y <- (max_y_value - min_y_value + 1) / bin_amount_per_axis
  plot_data_aux[, y_cell := ((y_column - min_y_value) %/% grid_size_y) + 1]
  
  # Density by cell
  point_density <- plot_data_aux[, .N, by = .(y_cell)]
  plot_data_aux <- merge(plot_data_aux, point_density[, .(y_cell, N)], by = "y_cell")
  
  # Plot
  p <- ggplot(plot_data_aux[order(plot_data_aux$y_column),]) +
    geom_jitter(aes(x = x_column, y = y_column, colour = log10(N + 1)), 
                width = 0.2, size = geom_point_size, shape = geom_point_shape) +
    scale_color_scico(palette = gradient_palette, guide = "colorbar") +
    scale_y_continuous(breaks = trans(c(0, 2500, 5000, 7500, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000), y_factors),
                       labels = c("0", "2500", "5000", "7500", "10000", "20000", "30000", "40000", "50000",  "60000", "70000", "80000", "90000")) +
    theme_bw() +
    theme(axis.title = element_text(size = axis_title_size, colour = "black"),
          axis.text = element_text(size = axis_text_size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(axis_ticks_length, "cm"),
          panel.background = element_rect(linewidth = panel_background_size, colour = "black"),
          legend.position = legend_position)
  
  # Plot personalization
  if (plot_title != "") {
    p <- p + ggtitle(label = plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (x_label != "") {
    p <- p + xlab(x_label)
  }
  if (y_label != "") {
    p <- p + ylab(y_label)
  }
  
  # Lims
  if (!is.null(y_lim_min) & !is.null(y_lim_max)) {
    p <- p + coord_cartesian(ylim = c(y_lim_min, y_lim_max))
  }
  
  return(p)
}

#############################-
####### Main function #######
#############################

plot_density <- function(array_config, project_folder, output_folder, smoothing, y_factors){
  
  #### Datasets ####
  plot_data_result <- process_files(array_config = array_config, project_folder = project_folder, smoothing = smoothing)
  has_negatives <- "NE" %in% array_config$types
  
  if (has_negatives){
    NE_raw <- plot_data_result[[1]][(plot_data_result[[1]]$type=="NE" & plot_data_result[[1]]$source!="LE"),]
    NE_norm <- plot_data_result[[2]][(plot_data_result[[2]]$type=="NE" & plot_data_result[[2]]$source!="LE"),]
    NE_smooth <- plot_data_result[[3]][(plot_data_result[[3]]$type=="NE" & plot_data_result[[3]]$source!="LE"),]
    
    max_NE <- round(max(NE_raw$Signal), 0)
  }
  
  PO_raw <- plot_data_result[[1]][(plot_data_result[[1]]$type=="PO" & plot_data_result[[1]]$source!="LE"),]
  PO_norm <- plot_data_result[[2]][(plot_data_result[[2]]$type=="PO" & plot_data_result[[2]]$source!="LE"),]
  PO_smooth <- plot_data_result[[3]][(plot_data_result[[3]]$type=="PO" & plot_data_result[[3]]$source!="LE"),]
  rm(plot_data_result)
  
  max_PO <- round(max(PO_raw$Signal), 0)
  
  #### Plots
  woaxis <- theme(
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
  )
  
  if(has_negatives){
  NE_raw_plot <- DensityJitterplot(NE_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                   geom_point_shape = 16, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                   y_lim_min = 0, y_lim_max = max_NE, legend_position = "none", y_factors = y_factors)
  
  NE_norm_plot <- DensityJitterplot(NE_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                    geom_point_shape = 16, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                    y_lim_min = 0, y_lim_max = max_NE, legend_position = "none", y_factors = y_factors)
  
  NE_smooth_plot <- DensityJitterplot(NE_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 2, 
                                      geom_point_shape = 16, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                      y_lim_min = 0, y_lim_max = max_NE,y_factors = y_factors)
  
  
  Negatives <- NE_raw_plot + (NE_norm_plot + woaxis) + (NE_smooth_plot + woaxis) + plot_annotation(title = "Negatives")
  
  #Save plot
  ggsave(filename = output_file_name(output_folder, array_config$id, "signal_distribution", smoothing = smoothing, sign = "negatives", ext = "pdf"), plot = Negatives, device = "pdf", width = 10, height = 8)
  
  }
  
  PO_raw_plot <- DensityJitterplot(PO_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                   geom_point_shape = 16, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                   y_lim_min = 0, y_lim_max = max_PO, legend_position = "none", y_factors = y_factors)
  
  PO_norm_plot <- DensityJitterplot(PO_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                    geom_point_shape = 16, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                    y_lim_min = 0, y_lim_max = max_PO, legend_position = "none", y_factors = y_factors)
  
  PO_smooth_plot <- DensityJitterplot(PO_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 2, 
                                      geom_point_shape = 16, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                      y_lim_min = 0, y_lim_max = max_PO,y_factors = y_factors)
  
  
  Positives <- PO_raw_plot + (PO_norm_plot + woaxis) + (PO_smooth_plot + woaxis) + plot_annotation(title = "Positives")
  
  #Save plot
  ggsave(filename = output_file_name(output_folder, array_config$id, "signal_distribution", smoothing = smoothing, sign = "positives", ext = "pdf"), plot = Positives, device = "pdf", width = 10, height = 8)
  
}


