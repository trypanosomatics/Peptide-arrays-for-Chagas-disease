############################-
#### Required libraries ####
############################-
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}
if (!require(readr, quietly = TRUE)) {
  writeLines("Installing library 'readr' for R")
  install.packages("readr", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(readr)
}
if (!require(dplyr, quietly = TRUE)) {
  writeLines("Installing library 'dplyr' for R")
  install.packages("dplyr", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(dplyr)
}
if (!require(tidyr, quietly = TRUE)) {
  writeLines("Installing library 'tidyr' for R")
  install.packages("tidyr", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(tidyr)
}
if (!require(ggplot2, quietly = TRUE)) {
  writeLines("Installing library 'ggplot2' for R")
  install.packages("ggplot2", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(ggplot2)
}
if (!require(scico, quietly = TRUE)) {
  writeLines("Installing library 'scico' for R")
  install.packages("scico", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(scico)
}
if (!require(patchwork, quietly = TRUE)) {
  writeLines("Installing library 'patchwork' for R")
  install.packages("patchwork", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(patchwork)
}

############################-
#### Auxiliar functions ####
############################-

process_files <- function(project_folder, type_data, raw_data_suffix, normalized_data_suffix, smoothed_data_suffix) {
  
  path_to_files <- ""
  
  directory_numbers <- list(
    pools = c(raw = "02", normalized = "01", smoothed = "02"),
    individual_serums = c(raw = "04", normalized = "05", smoothed = "06")
  )
  
  steps <- c("raw", "normalized", "smoothed")
  all_data_raw <- data.frame()
  all_data_normalized <- data.frame()
  all_data_smoothed <- data.frame()
  
  for(step in steps){
    dir_num <- directory_numbers[[type_data]][[step]]
    if (step == "raw") {
      path_to_files <- sprintf("%s/inputs/%s_%s_raw_data", project_folder, dir_num, type_data)
      files <- list.files(path = path_to_files, pattern = raw_data_suffix, full.names = TRUE)
      
      for (file in files) {
        # file <- files[1]
        data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
        
        data_long <- data %>%
          pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")
        
        if(type_data == "pools"){
          data_long$source <- sapply(strsplit(data_long[[2]], "_"), `[`, 1)
          data_long$type <- sapply(strsplit(data_long[[2]], "_"), `[`, 2)
          data_long$type <- as.factor(data_long$type)
          data_long$source <- as.factor(data_long$source)}
        else{
          data_long$source <- unlist(strsplit(tail(unlist(strsplit(file, split = "/")), n=1),split = "_"))[1]
          data_long$type <- unlist(strsplit(tail(unlist(strsplit(file, split = "/")), n=1),split = "_"))[3]
        }
        
        all_data_raw <- bind_rows(all_data_raw, data_long)
      }
      
      rm(data, data_long)
      
    } else if (step == "normalized") {
      path_to_files <- sprintf("%s/outputs/%s_%s_normalized_data", project_folder, dir_num, type_data)
      
      files <- list.files(path = path_to_files, pattern = normalized_data_suffix, full.names = TRUE)

      for (file in files) {
        data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
        
        data_long <- data %>%
          pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")
        
        if(type_data == "pools"){
          data_long$source <- sapply(strsplit(data_long[[2]], "_"), `[`, 1)
          data_long$type <- sapply(strsplit(data_long[[2]], "_"), `[`, 2)
          data_long$type <- as.factor(data_long$type)
          data_long$source <- as.factor(data_long$source)}
        else{
          data_long$source <- unlist(strsplit(tail(unlist(strsplit(file, split = "/")), n=1),split = "_"))[1]
          data_long$type <- unlist(strsplit(tail(unlist(strsplit(file, split = "/")), n=1),split = "_"))[3]
        }
        
        all_data_normalized <- bind_rows(all_data_normalized, data_long)
      }
      
      rm(data, data_long)
      
    } else if (step == "smoothed") {
      path_to_files <- sprintf("%s/outputs/%s_%s_smoothed_data", project_folder, dir_num, type_data)
      
      files <- list.files(path = path_to_files, pattern = smoothed_data_suffix, full.names = TRUE)

      for (file in files) {
        # file <- files[1]
        data <- read_delim(file, delim = "\t", trim_ws = TRUE, show_col_types = FALSE)
        
        all_data_smoothed <- bind_rows(all_data_smoothed, data)
      }
      
      if (type_data == "individual_serums"){
      all_data_smoothed$source <- sapply(strsplit(all_data_smoothed[,"source"], "_"), `[`, 1)
      }
      
      all_data_smoothed$type <- as.factor(all_data_smoothed$type)
      all_data_smoothed$source <- as.factor(all_data_smoothed$source)
      
      rm(data, file, files, path_to_files)
    }
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

plot_density <- function(project_folder, type_data, raw_data_suffix, normalized_data_suffix, smoothed_data_suffix,
                         y_factors, plot_data){
  
  #### Datasets ####
  plot_data_result <- process_files(project_folder = project_folder, type_data = type_data,
                                    raw_data_suffix = raw_data_suffix, normalized_data_suffix = normalized_data_suffix,
                                    smoothed_data_suffix = smoothed_data_suffix)
  
  if (type_data == "pools"){
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
  
  if(type_data == "pools"){
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
  ggsave(filename = sprintf("%s/outputs/08_plots/Signal_negatives_%s.pdf", project_folder, type_data), plot = Negatives, device = "pdf", width = 10, height = 8)
  
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
  ggsave(filename = sprintf("%s/outputs/08_plots/Signal_positives_%s.pdf", project_folder, type_data), plot = Positives, device = "pdf", width = 10, height = 8)
  
}


