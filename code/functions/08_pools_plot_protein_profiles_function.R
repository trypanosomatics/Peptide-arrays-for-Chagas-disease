############################-
#### Required libraries ####
############################-
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}
if (!require(ggplot2, quietly = TRUE)) {
  writeLines("Installing library 'ggplot2' for R")
  install.packages("ggplot2", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(ggplot2)
}
if (!require(grid, quietly = TRUE)) {
  writeLines("Installing library 'grid' for R")
  install.packages("grid", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(grid)
}
if (!require(gridExtra, quietly = TRUE)) {
  writeLines("Installing library 'gridExtra' for R")
  install.packages("gridExtra", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(gridExtra)
}
if (!require(patchwork, quietly = TRUE)) {
  writeLines("Installing library 'patchwork' for R")
  install.packages("patchwork", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(patchwork)
}

plot_parameters <- function(geom_point_size_parameter, geom_line_size_parameter,
                            title_size_parameter, axis_title_size_parameter, axis_text_size_parameter, 
                            axis_ticks_length_parameter, panel_background_size_parameter, 
                            plot_margin_l_parameter, plot_margin_r_parameter,
                            geom_errorbar_width_parameter, geom_errorbar_size_parameter,
                            legend_x_parameter, legend_y_parameter, legend_text_size_parameter, 
                            legend_element_width_parameter, legend_element_height_parameter)
{
  geom_point_size <<- geom_point_size_parameter
  geom_line_size <<- geom_line_size_parameter
  title_size <<- title_size_parameter
  axis_title_size <<- axis_title_size_parameter
  axis_text_size <<- axis_text_size_parameter
  axis_ticks_length <<- axis_ticks_length_parameter
  panel_background_size <<- panel_background_size_parameter
  plot_margin_l <<- plot_margin_l_parameter
  plot_margin_r <<- plot_margin_r_parameter
  geom_errorbar_width <<- geom_errorbar_width_parameter
  geom_errorbar_size <<- geom_errorbar_size_parameter
  legend_x <<- legend_x_parameter
  legend_y <<- legend_y_parameter
  legend_text_size <<- legend_text_size_parameter
  legend_element_width <<- legend_element_width_parameter
  legend_element_height <<- legend_element_height_parameter
}

#############################-
#### Auxiliary Functions ####
#############################-

prepare_plot_data <- function(
    project_folder, 
    input_folder, 
    profile_data_suffix, 
    sources, 
    protein, 
    output_folder, 
    sd_multiplier_for_cutoff, 
    only_proteins_above, 
    only_proteins_below,
    output_suffix
) {
  
  # LOAD THE DATA NECESSARY
  input_files <- list.files(input_folder, pattern = profile_data_suffix, full.names = TRUE)
  input_files <- input_files[grepl(paste(sources, collapse = "|"), input_files)]
  
  if (length(input_files) == 0) {
    stop("No files found matching the given profile_data_suffix and sources.")
  }
  
  # Get data
  all_plot_data <- do.call(rbind, lapply(input_files, read.delim, stringsAsFactors = FALSE))
  
  # Filter by protein
  if (!is.null(protein)) {
    all_plot_data <- all_plot_data[all_plot_data$protein == protein, ]
    if (nrow(all_plot_data) == 0) {
      stop(paste("No data found for the specified protein:", protein))
    }
  }
  
  # Order
  setDT(all_plot_data)
  plot_order <- all_plot_data[type == "PO", .(signal = sum(mean_smoothed_signal)), by = .(protein, start)]
  plot_order <- plot_order[, .(signal = max(signal)), by = .(protein)]
  plot_order <- plot_order[order(-signal)]
  
  proteins_to_plot <- all_plot_data[type == "PO", .(max_signal = max(mean_smoothed_signal)), by = .(protein)]
  if (only_proteins_above > 0) {
    proteins_to_plot <- proteins_to_plot[max_signal >= only_proteins_above]
  }
  if (only_proteins_below > 0) {
    proteins_to_plot <- proteins_to_plot[max_signal <= only_proteins_below]
  }
  proteins_to_plot <- plot_order[protein %in% proteins_to_plot$protein]$protein
  
  # Save the index
  index_data <- data.table(protein = proteins_to_plot)
  index_data[, page := .I]
  
  output_file_prefix <- if (is.null(protein)) {
    sprintf("protein_profiles_cruzi_%sSD", sd_multiplier_for_cutoff)
  } else {
    sprintf("protein_profiles_cruzi_%sSD_selected_proteins", sd_multiplier_for_cutoff)
  }
  
  output_file_name <- if (is.null(output_suffix)) {
    sprintf("%s/%s_index.tsv", output_folder, output_file_prefix)
  } else {
    sprintf("%s/%s_%s_index.tsv", output_folder, output_file_prefix, output_suffix)
  }
  
  write.table(index_data, 
              file = output_file_name, 
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              quote = TRUE)
  
  
  return(list(all_plot_data = all_plot_data, plot_order = plot_order, proteins_to_plot = proteins_to_plot))
}

#############################-
####    Main Function    ####
#############################-

plot_proteins = function(project_folder,
                         input_folder,
                         design_file,
                         sources,
                         profile_data_suffix,
                         output_folder,
                         protein,
                         only_proteins_above,
                         only_proteins_below,
                         sd_multiplier_for_cutoff,
                         output_suffix,
                         fixed_scale) {
  
  plot_data_result <- prepare_plot_data(
    project_folder = project_folder,
    input_folder = input_folder,
    protein = protein,
    profile_data_suffix = profile_data_suffix,
    sources = sources,
    output_folder = output_folder,
    sd_multiplier_for_cutoff = sd_multiplier_for_cutoff,
    only_proteins_above = only_proteins_above,
    only_proteins_below = only_proteins_below,
    output_suffix = output_suffix
  )
  
  all_plot_data <- plot_data_result$all_plot_data
  plot_order <- plot_data_result$plot_order
  proteins_to_plot <- plot_data_result$proteins_to_plot
  
  ############################-
  #### CONFIGURE THE PLOT ####
  ############################-
  
  #Add the errors
  all_plot_data[, error_max := mean_smoothed_signal + sd_smoothed_signal]
  all_plot_data[, error_min := mean_smoothed_signal - sd_smoothed_signal]
  
  #Calculate the max value for the fixed scale
  max_plot <- max(as.numeric(all_plot_data$error_max))
  
  # CALCULATE CUTOFF
  normalization_global_statistics_file <- sprintf("%s/outputs/01_pools_normalized_data/global_statistics.tsv", project_folder)
  normalization_global_statistics <- fread(normalization_global_statistics_file, header = TRUE, sep = "\t", na.strings = NULL)
  
  mode_aux <- normalization_global_statistics$mode
  sd_aux <- normalization_global_statistics$sd
  
#rm  # PLOT DATA
#  plot_data_folder <- input_folder
#  plot_data_file_suffix <- profile_data_suffix
  
  has_negative_data <- 0
  positive_serum_label <- "Positive serum"
  negative_serum_label <- "Negative serum"
  threshold_label <- "Threshold"
  
  color_palette <- c("#133859ff", "#e67960ff", "black")
  types <- c(negative_serum_label, positive_serum_label, threshold_label)
  types_colors <- setNames(color_palette, types)
  
  #PDF file
  output_file_prefix <- if (is.null(protein)) {
    sprintf("protein_profiles_cruzi_%sSD", sd_multiplier_for_cutoff)
  } else {
    sprintf("protein_profiles_cruzi_%sSD_selected_proteins", sd_multiplier_for_cutoff)
  }
  
  output_file <- if (is.null(output_suffix)) {
    sprintf("%s/%s.pdf", output_folder, output_file_prefix)
  } else {
    sprintf("%s/%s_%s.pdf", output_folder, output_file_prefix, output_suffix)
  }
  
  graphics.off()
  pdf(output_file, width = 22, height = 11)
  
  for (protein_for in proteins_to_plot) {
    writeLines(sprintf("Plotting %s...", protein_for))
    
    sub_plot_data <- all_plot_data[protein == protein_for]
    
    plots <- list()
    for (source_for in sources) {
      
      sub_sub_plot_data <- sub_plot_data[source == source_for]
      if (nrow(sub_sub_plot_data) == 0) {
        writeLines(sprintf("No data for source %s. Skipping...", source_for))
        next
      }
      
      sub_sub_plot_data[type == "PO", type := positive_serum_label]
      sub_sub_plot_data[type == "NE", type := negative_serum_label]
      
      cutoffs_to_plot <- mode_aux + sd_multiplier_for_cutoff * sd_aux
      
      #Plot
      p <- ggplot(sub_sub_plot_data, aes(x = start, y = mean_smoothed_signal, color = type)) +
        geom_point() +  
        geom_line() +
        geom_errorbar(aes(ymin = error_min, ymax = error_max), width = 0.2) + 
        scale_color_manual(values = types_colors) +
        geom_segment(data = data.frame(x = min(sub_sub_plot_data$start), 
                                       xend = max(sub_sub_plot_data$start),
                                       y = cutoffs_to_plot, 
                                       yend = cutoffs_to_plot, 
                                       type = threshold_label),
                     aes(x = x, xend = xend, y = y, yend = yend, color = type), 
                     linetype = "dashed", size = 0.5)+
        labs(title = source_for, x = "Peptide Position", y = "Signal", color = protein_for) +
        theme_bw() +
        theme(panel.grid.major = element_line(linewidth = 0.1, colour = "gray92"),
              panel.grid.minor = element_line(linewidth = 0.1, colour = "gray92"),
              plot.title = element_text(size = 20),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 
        
      if (source_for == "LE") {
        p <- p + theme(panel.background = element_rect(fill = "#fefbec"),
                       legend.position = "none")
      }
      
      if (fixed_scale) {
        p <- p + ylim(0, max_plot)
      }
      
      plots[[source_for]] <- p
    }
    
    final_plot <- wrap_plots(plots) + plot_layout(guides = "collect")
    
    print(final_plot)
  }
  
  dev.off()
  graphics.off()
}  