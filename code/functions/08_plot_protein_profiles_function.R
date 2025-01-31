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

##################################-
#### PLOT AUXILIARY FUNCTIONS ####
##################################-

plotProteinProfile_v3 <- function(plot_data, x_column_name, y_column_name, group_column_name, foreground_group, background_group,
                                  cutoffs_y = c(),
                                  plot_title = "", x_label = "X", y_label = "Y", legend_labels = c(), line_colors = c(),
                                  manual_cutoff_color = c(), manual_cutoff_type = c(), manual_cutoff_size = c(),
                                  show_errors = 0, error_min_column_name = "", error_max_column_name = "",
                                  show_legend = 0,
                                  fixed_scale = 0, fixed_scale_max = 0) {
  
  plotProteinProfileDetailed_v3(plot_data = plot_data, x_column_name = x_column_name, y_column_name = y_column_name, group_column_name = group_column_name, foreground_group = foreground_group, background_group = background_group,
                                cutoffs_y = cutoffs_y,
                                plot_title = plot_title, x_label = x_label, y_label = y_label,
                                geom_point_size = geom_point_size,
                                geom_line_size = geom_line_size, 
                                title_size = title_size, axis_title_size = axis_title_size, axis_text_size = axis_text_size, axis_ticks_length = axis_ticks_length, panel_background_size = panel_background_size, plot_margin_l = plot_margin_l, plot_margin_r = plot_margin_r,
                                legend_labels = legend_labels, line_colors = line_colors,
                                manual_cutoff_color = manual_cutoff_color, manual_cutoff_type = manual_cutoff_type, manual_cutoff_size = manual_cutoff_size,
                                show_errors = show_errors, error_min_column_name = error_min_column_name, error_max_column_name = error_max_column_name, geom_errorbar_width = geom_errorbar_width, geom_errorbar_size = geom_errorbar_size,
                                show_legend = show_legend, legend_x = legend_x, legend_y = legend_y, legend_text_size = legend_text_size, legend_element_width = legend_element_width, legend_element_height = legend_element_height,
                                fixed_scale = fixed_scale, fixed_scale_max = fixed_scale_max)  
}

plotProteinProfileDetailed_v3 <- function(plot_data, x_column_name, y_column_name, group_column_name, foreground_group, background_group,
                                          cutoffs_y = c(),
                                          plot_title = "", x_label = "X", y_label = "Y",
                                          geom_point_size = geom_point_size,
                                          geom_line_size = geom_line_size,
                                          title_size = title_size, axis_title_size = axis_title_size, axis_text_size = axis_text_size, axis_ticks_length = axis_ticks_length, panel_background_size = panel_background_size, plot_margin_l = plot_margin_l, plot_margin_r = plot_margin_r,
                                          legend_labels = c(), line_colors = c(),
                                          manual_cutoff_color = c(), manual_cutoff_type = c(), manual_cutoff_size = c(), cutoff_above_plot = 1,
                                          show_errors = 0, error_min_column_name = "", error_max_column_name = "", geom_errorbar_width = geom_errorbar_width, geom_errorbar_size = geom_errorbar_size,
                                          show_legend = 0, legend_x = legend_x, legend_y = legend_y, legend_text_size = legend_text_size, legend_element_width = legend_element_width, legend_element_height = legend_element_height,
                                          fixed_scale = 0, fixed_scale_max = 0) {
  plot_data_aux <- plot_data
  plot_data_aux[1] <- plot_data_aux[1]
  
  setnames(plot_data_aux, x_column_name, "x_column")
  setnames(plot_data_aux, y_column_name, "y_column")
  setnames(plot_data_aux, group_column_name, "group_column")   
  
  plot_data_aux$group_column <- factor(plot_data_aux$group_column, levels = c(background_group, foreground_group))
  
  if (show_errors) {
    setnames(plot_data_aux, error_min_column_name, "error_min_column")
    setnames(plot_data_aux, error_max_column_name, "error_max_column")
  }
  
  p <- ggplot()
  
  if ((length(cutoffs_y) > 0) & (cutoff_above_plot == 0)) {
    for (cutoff_i in 1:length(cutoffs_y)) {
      line_color_for <- "blue"
      line_type_for <- "dashed"
      line_size_for <- 0.5
      if (length(manual_cutoff_color) > 0) {
        line_color_for <- manual_cutoff_color[cutoff_i]
      }
      if (length(manual_cutoff_type) > 0) {
        line_type_for <- manual_cutoff_type[cutoff_i]
      }
      if (length(manual_cutoff_size) > 0) {
        line_size_for <- manual_cutoff_size[cutoff_i]
      }
      
      p <- p + geom_segment(data = plot_data_aux, x = min(plot_data_aux$x_column), y = cutoffs_y[cutoff_i], 
                            xend = max(plot_data_aux$x_column), yend = cutoffs_y[cutoff_i], 
                            colour = line_color_for, linetype = line_type_for, size = line_size_for)
    }
  }
  
  for (group_for in c(background_group, foreground_group)) {
    p <- p + geom_line(data = plot_data_aux[group_column == group_for], aes_string(x = "x_column", y = "y_column", color = "group_column"), size = geom_line_size) +
      geom_point(data = plot_data_aux[group_column == group_for], aes_string(x = "x_column", y = "y_column", color = "group_column"), size = geom_point_size)
    
    if (show_errors) {
      p <- p + geom_errorbar(data = plot_data_aux[group_column == group_for], aes_string(x = "x_column", y = "y_column", color = "group_column", ymin = "error_min_column", ymax = "error_max_column"), size = geom_errorbar_size, width = geom_errorbar_width)
    }
  }
  
  p <- setTheme_BW(p = p)
  p <- formatAxis(p = p, axis_text_size = axis_text_size, axis_ticks_length = axis_ticks_length)
  p <- formatBackground(p = p,
                        panel_grid_major_size = 0, panel_grid_minor_size = 0, 
                        panel_background_size = panel_background_size)
  p <- setPlotMargin(p = p, plot_margin_right = plot_margin_r, plot_margin_left = plot_margin_l)
  
  p <- setPlotTitle(p = p, plot_title = plot_title, title_size = title_size)
  p <- setXLabel(p = p, x_label = x_label, x_axis_title_size = axis_title_size)
  p <- setYLabel(p = p, y_label = y_label, y_axis_title_size = axis_title_size)
  
  if (show_legend) {
    p <- addLegend(p = p, legend_text_size = legend_text_size, 
                   legend_x = legend_x, legend_y = legend_y,
                   legend_element_width = legend_element_width, legend_element_height = legend_element_height)
  } else {
    p <- addLegend(p = p, legend_text_size = 0)
  }
  
  if (fixed_scale) {
    p <- addYScale(p = p, y_scale_max = fixed_scale_max)
  } else {
    if (show_errors) {
      p <- addYScale(p = p, y_scale_max = max(plot_data_aux$error_max_column))
    } else {
      p <- addYScale(p = p, y_scale_max = max(plot_data_aux$y_column))
    }
  }
  
  if ((length(cutoffs_y) > 0) & (cutoff_above_plot == 1)) {
    for (cutoff_i in 1:length(cutoffs_y)) {
      line_color_for <- "blue"
      line_type_for <- "dashed"
      line_size_for <- 0.5
      if (length(manual_cutoff_color) > 0) {
        line_color_for <- manual_cutoff_color[cutoff_i]
      }
      if (length(manual_cutoff_type) > 0) {
        line_type_for <- manual_cutoff_type[cutoff_i]
      }
      if (length(manual_cutoff_size) > 0) {
        line_size_for <- manual_cutoff_size[cutoff_i]
      }
      
      p <- p + geom_segment(data = plot_data_aux, x = min(plot_data_aux$x_column), y = cutoffs_y[cutoff_i], 
                            xend = max(plot_data_aux$x_column), yend = cutoffs_y[cutoff_i], 
                            colour = line_color_for, linetype = line_type_for, size = line_size_for)
    }
  }
  
  p
}

addFill <- function(p, fill_alpha) {
  if (fill_alpha > 0) {
    p <- p + geom_area(alpha = fill_alpha, position = "identity")
  }
  
  p
}
addLine <- function(p, geom_line_size, geom_line_linetype = "") {
  if (geom_line_size > 0) {
    if (geom_line_linetype != "") {
      p <- p + geom_line(linewidth = geom_line_size, linetype = geom_line_linetype)
    } else {
      p <- p + geom_line(linewidth = geom_line_size)
    }
  }
  
  p
}
addLineAndPoint <- function(p, geom_line_size, geom_line_linetype = "", geom_point_size) {
  if (geom_line_size > 0 & geom_point_size > 0) {
    if (geom_line_linetype != "") {
      p <- p + geom_line(linewidth = geom_line_size, linetype = geom_line_linetype) + geom_point(size = geom_point_size)
    } else {
      p <- p + geom_line(linewidth = geom_line_size) + geom_point(size = geom_point_size)
    }
  }
  
  p
}
addPath <- function(p, geom_line_size) {
  if (geom_line_size > 0) {
    p <- p + geom_path(aes(group = 1), size = geom_line_size, lineend = "round", linejoin = "mitre")
  }
  
  p
}
addPoint <- function(p, geom_point_size) {
  if (geom_point_size > 0) {
    p <- p + geom_point(size = geom_point_size)
  }
  
  p
}
addFacets <- function(p, row_variable, column_variable = "", strip_position = "left") {
  if (column_variable == "") {
    column_variable <- "."
  }
  if (strip_position == "left") {
    p <- p + facet_grid(as.formula(paste(row_variable,"~", column_variable)), scales = "free_y", switch = "y")
  } else {
    p <- p + facet_grid(as.formula(paste(row_variable,"~", column_variable)), scales = "free_y")
  }
  
  p
}
setTheme_BW <- function(p) {
  p <- p + theme_bw()
  
  p
}
formatAxis <- function(p,
                       axis_text_size, axis_text_color = "black",
                       axis_ticks_length, axis_ticks_color = "black",
                       axis_line_size = 0, axis_line_color = "black") {
  if (axis_text_size > 0) {
    p <- p + theme(axis.text = element_text(size = axis_text_size, colour = axis_text_color))    
  } else {
    p <- p + theme(axis.text = element_blank())
  }
  if (axis_ticks_length > 0) {
    p <- p + theme(axis.ticks = element_line(colour = axis_ticks_color),
                   axis.ticks.length = unit(axis_ticks_length, "cm"))    
  } else {
    p <- p + theme(axis.ticks = element_blank())
  }
  if (axis_line_size > 0) {
    p <- p + theme(axis.line = element_line(linewidth = axis_line_size, colour = axis_line_color))
  } else {
    p <- p + theme(axis.line = element_blank())
  } 
  
  p
}
formatBackground <- function(p,
                             panel_grid_major_size, panel_grid_major_color = "black",
                             panel_grid_minor_size, panel_grid_minor_color = "black",
                             panel_background_size, panel_background_color = "black",
                             aspect_ratio = 1/1) {
  if (panel_grid_major_size > 0) {
    p <- p + theme(panel.grid.major = element_line(linewidth = panel_grid_major_size, colour = panel_grid_major_color))
  } else {
    p <- p + theme(panel.grid.major = element_blank())
  }
  if (panel_grid_minor_size > 0) {
    p <- p + theme(panel.grid.minor = element_line(linewidth = panel_grid_minor_size, colour = panel_grid_minor_color))
  } else {
    p <- p + theme(panel.grid.minor = element_blank())
  }
  if (panel_background_size > 0) {
    p <- p + theme(panel.background = element_rect(linewidth = panel_background_size, colour = panel_background_color))
  } else {
    p <- p + theme(panel.border = element_blank(),
                   panel.background = element_blank())
  }
  if (aspect_ratio != 0) {
    p <- p + theme(aspect.ratio = aspect_ratio)    
  }
  
  p
}
formatStrips <- function(p, strip_text_size, strip_text_color = "black", strip_text_angle = 0, label_position = "left") {
  if (strip_text_size > 0) {
    if (label_position == "left") {
      p <- p + theme(strip.text = element_text(size = strip_text_size, colour = strip_text_color),
                     strip.text.y.left = element_text(angle = strip_text_angle),
                     strip.background = element_blank())    
    } else {
      p <- p + theme(strip.text = element_text(size = strip_text_size, colour = strip_text_color),
                     strip.text.y.right = element_text(angle = strip_text_angle),
                     strip.background = element_blank())
    }
  }
  
  p
}
setPlotMargin <- function(p, plot_margin_top = 0, plot_margin_right = 0, plot_margin_bottom = 0, plot_margin_left = 0) {
  #top, right, bottom, and left margins
  p <- p + theme(plot.margin = unit(c(plot_margin_top, plot_margin_right, plot_margin_bottom, plot_margin_left), "cm"))
  
  p
}
expandPlotXBorder <- function(p, expand_size) {
  p <- p + scale_x_continuous(expand = c(expand_size, 0))
  
  p
}
setPlotTitle <- function(p, plot_title, title_size = 0, title_color = "black", title_is_centered = 1) {
  if (plot_title != "") {
    p <- p + ggtitle(label = plot_title) + theme(plot.title = element_text(size = title_size, colour = title_color))
    
    if (title_is_centered) {
      p <- p + theme(plot.title = element_text(hjust = 0.5))
    }
  } else {
    p <- p + ggtitle(label = "") + theme(plot.title =  element_blank())
  }
  
  p
}
setXLabel <- function(p, x_label, x_axis_title_size = 0, x_axis_title_color = "black") {
  if (x_label != "") {
    p <- p + xlab(x_label) + theme(axis.title.x = element_text(size = x_axis_title_size, colour = x_axis_title_color))
  } else {
    p <- p + xlab("") + theme(axis.title.x = element_blank())
  }
  
  p
}
setYLabel <- function(p, y_label, y_axis_title_size = 0, y_axis_title_color = "black") {
  if (y_label != "") {
    p <- p + ylab(y_label) + theme(axis.title.y = element_text(size = y_axis_title_size, colour = y_axis_title_color))
  } else {
    p <- p + ylab("") + theme(axis.title.y = element_blank())
  }
  
  p
}

addManualColorScale <- function(p, line_colors, is_legend_visible = 1, legend_labels = c()) {
  if (length(line_colors) > 0) {
    if (is_legend_visible & (length(legend_labels) > 0)) {
      p <- p + scale_color_manual(values = line_colors,
                                  labels = legend_labels)
    } else {
      p <- p + scale_color_manual(values = line_colors)
    }    
  }
  
  p
}
addManualLinetypeScale <- function(p, line_linetypes, is_legend_visible = 1, legend_labels = c()) {
  if (length(line_linetypes) > 0) {
    if (is_legend_visible & (length(legend_labels) > 0)) {
      p <- p + scale_linetype_manual(values = line_linetypes,
                                     labels = legend_labels)
    } else {
      p <- p + scale_linetype_manual(values = line_linetypes)
    }    
  }
  
  p
}
addManualFillScale <- function(p, fill_colors, is_legend_visible = 1, legend_labels = c()) {
  if (length(fill_colors) > 0) {
    #If you don't have to show the legend you don't care about the labels
    if (is_legend_visible & (length(legend_labels) > 0)) {
      p <- p + scale_fill_manual(values = fill_colors,
                                 labels = legend_labels)
    } else {
      p <- p + scale_fill_manual(values = fill_colors)
    }    
  }
  
  p
}
addErrorBars <- function(p, error_min_column, error_max_column = "", geom_errorbar_size = 0, geom_errorbar_width = 0) {
  if ((error_min_column != "") & (error_max_column != "") & (geom_errorbar_size > 0)) {
    p <- p + geom_errorbar(aes_string(ymin = error_min_column, ymax = error_max_column), 
                           size = geom_errorbar_size, width = geom_errorbar_width)
  }
  
  p
}
addLegend <- function(p, legend_text_size, legend_text_color = "black", 
                      legend_x = 0, legend_y = 0,
                      legend_element_width = 0, legend_element_height = 0) {
  if (legend_text_size > 0) {
    p <- p + theme(legend.position = c(legend_x, legend_y),
                   legend.title = element_blank(),
                   legend.text = element_text(size = legend_text_size, colour = legend_text_color),
                   legend.background = element_rect(fill = alpha('white', 0)))
    p <- p + guides(fill = guide_legend(keywidth = legend_element_width, keyheight = legend_element_height),
                    linetype = guide_legend(keywidth = legend_element_width, keyheight = legend_element_height),
                    colour = guide_legend(keywidth = legend_element_width, keyheight = legend_element_height))    
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  p
}
addYScale <- function(p, y_scale_min = 0, y_scale_max) {
  if (y_scale_max > y_scale_min) {
    p <- p + coord_cartesian(ylim = c(y_scale_min, y_scale_max))
  }
  
  p
}
addHorizontalLines <- function(p, 
                               min_x, max_x, heights,
                               default_lines_color = "blue", default_lines_linetype = "dashed", default_lines_size = 0.5,
                               manual_lines_color = c(), manual_lines_type = c(), manual_lines_size = c()) {
  if (length(heights) > 0) {
    for (height_i in 1:length(heights)) {
      height_for <- heights[height_i]
      
      line_color_for <- default_lines_color
      line_type_for <- default_lines_linetype
      line_size_for <- default_lines_size
      if (length(manual_lines_color) > 0) {
        line_color_for <- manual_lines_color[height_i]
      }
      if (length(manual_lines_type) > 0) {
        line_type_for <- manual_lines_type[height_i]
      }
      if (length(manual_lines_size) > 0) {
        line_size_for <- manual_lines_size[height_i]
      }
      
      p <- p + geom_segment(x = min_x, y = height_for, xend = max_x, yend = height_for, colour = line_color_for, linetype = line_type_for, size = line_size_for)
    }
  }
  
  p
}


#######################-
#### MAIN FUNCTION ####-
#######################
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
                         output_suffix) {
  
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
  
  neg_line_color = "#f8766d"
  pos_line_color = "#00bfc4"
  cutoff_color = "black"
  
  #Add the errors
  all_plot_data[, error_max := mean_smoothed_signal + sd_smoothed_signal]
  all_plot_data[, error_min := mean_smoothed_signal - sd_smoothed_signal]
  
  #Calculate the max value for the fixed scale
  max_plot <- max(as.numeric(all_plot_data$error_max))
  
  #Aux Function for Plots
  extractLegend <- function(gg) {
    grobs <- ggplot_gtable(ggplot_build(gg))
    foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
    grobs$grobs[[foo]]
  }  
  
  # CALCULATE CUTOFF
  normalization_global_statistics_file <- sprintf("%s/outputs/01_pools_normalized_data/global_statistics.tsv", project_folder)
  normalization_global_statistics <- fread(normalization_global_statistics_file, header = TRUE, sep = "\t", na.strings = NULL)
  
  mode_aux <- normalization_global_statistics$mode
  sd_aux <- normalization_global_statistics$sd
  
  # PLOT DATA
  plot_data_folder <- input_folder
  plot_data_file_suffix <- profile_data_suffix
  
  output_file_prefix <- if (is.null(protein)) {
    sprintf("protein_profiles_cruzi_%sSD", sd_multiplier_for_cutoff)
  } else {
    sprintf("protein_profiles_cruzi_%sSD_selected_proteins", sd_multiplier_for_cutoff)
  }
  
  output_file_name <- if (is.null(output_suffix)) {
    sprintf("%s/%s.pdf", output_folder, output_file_prefix)
  } else {
    sprintf("%s/%s_%s.pdf", output_folder, output_file_prefix, output_suffix)
  }
  
  has_negative_data <- 0
  positive_serum_label <- "Positive serum"
  negative_serum_label <- "Negative serum"
  threshold_label <- "Threshold"
  pdf_height <- 21
  output_file <- output_file_name
  
  graphics.off()
  for (protein_for in proteins_to_plot) {
    writeLines(sprintf("Plotting %s...", protein_for))
    
    sub_plot_data <- all_plot_data[protein == protein_for]
    
    plots <- list()
    for (source_i in 1:length(sources)) {
      source_for <- sources[source_i]
      
      sub_sub_plot_data <- sub_plot_data[source == source_for]
      if (nrow(sub_sub_plot_data) == 0) {
        writeLines(sprintf("No data for source %s. Skipping...", source_for))
        next
      }
      
      sub_sub_plot_data[type == "PO", type := positive_serum_label]
      sub_sub_plot_data[type == "NE", type := negative_serum_label]
      
      plot_title <- source_for
      
      if (has_negative_data) {
        labels_aux <- c(negative_serum_label, positive_serum_label)
        colors_aux = c(neg_line_color, pos_line_color)
      } else {
        labels_aux <- c(positive_serum_label)
        colors_aux = c(pos_line_color)
      }
      
      #Extract legend of the first plot
      if (source_i == 1) {
        fake_data <- data.table(x = c(1:6),
                                y = rnorm(6),
                                cat = c(positive_serum_label, positive_serum_label,
                                        negative_serum_label, negative_serum_label,
                                        threshold_label, threshold_label))
        fake_data$cat <- factor(fake_data$cat, levels = c(positive_serum_label, negative_serum_label, threshold_label))
        p <- ggplot(data = fake_data, aes(x = x, y = y, color = cat, linetype = cat))
        p <- p + geom_point(size = 4) + geom_line(linewidth = 2)
        
        p <- p + theme_bw()
        p <- p + scale_color_manual(values = c(pos_line_color, neg_line_color, cutoff_color)) +
          scale_linetype_manual(values = c("solid", "solid", "dashed"))
        
        legend_element_width <- 4
        legend_element_height <- 3
        p <- p + theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20))
        p <- p + theme(legend.text = element_text(size = 20))
        p <- p + theme(legend.position = "right")
        
        p <- p + guides(colour = guide_legend(title = protein_for, keywidth = legend_element_width, keyheight = legend_element_height,
                                              override.aes = list(shape = c(16, 16, NA), linetype = c("solid", "solid", "dashed"),
                                                                  linewidth = c(2, 2, 1))),
                        linetype = guide_legend(title = protein_for, keywidth = legend_element_width, keyheight = legend_element_height,
                                                override.aes = list(shape = c(16, 16, NA), linetype = c("solid", "solid", "dashed"),
                                                                    linewidth = c(2, 2, 1))))
        
        if (source_i == 1 && nrow(fake_data) > 0) {
          legend_aux <- extractLegend(p)
        } else {
          legend_aux <- NULL
        }
        
      }
      
      cutoffs_to_plot <- mode_aux + sd_multiplier_for_cutoff * sd_aux
      
      #Plot
      p <- plotProteinProfile_v3(plot_data = sub_sub_plot_data, x_column_name ="start", y_column_name = "mean_smoothed_signal", group_column_name = "type", foreground_group = positive_serum_label, background_group = negative_serum_label,
                                 cutoffs_y = cutoffs_to_plot,
                                 plot_title = plot_title, x_label = "Peptide position", y_label = "Signal",
                                 legend_labels = labels_aux, line_colors = colors_aux,
                                 manual_cutoff_color = c(cutoff_color), manual_cutoff_type = c(), manual_cutoff_size = c(),
                                 show_errors = 1, error_min_column_name = "error_min", error_max_column_name = "error_max",
                                 show_legend = 0,
                                 fixed_scale = 1, fixed_scale_max = max_plot)
      p <- p + theme(panel.grid.major = element_line(linewidth = 0.1, colour = "gray92"))
      p <- p + theme(panel.grid.minor = element_line(linewidth = 0.1, colour = "gray92"))
      p <- p + theme(plot.title = element_text(size = 20))
      p <- p + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
      if (source_for == "LE") {
        p <- p + theme(panel.background = element_rect(fill = "#fefbec"))
      }
      plots[[source_i]] <- p
    }
    
    if (!is.null(legend_aux)) {
      plots[[length(plots) + 1]] <- legend_aux
    }
    
    #Initialize the PDF the first time
    if (protein_for == proteins_to_plot[1]) {
      pdf(file = output_file, width = 28, height = 14)
    }
    layout_aux <- rbind(c(1,2,5),
                        c(3,4,6))
    grid.arrange(grobs = plots, ncol = 4, layout_matrix = layout_aux)
  }
  dev.off()
  graphics.off()
  
}
