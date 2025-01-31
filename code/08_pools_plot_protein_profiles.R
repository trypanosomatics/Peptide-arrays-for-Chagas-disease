#!/usr/bin/env Rscript

# CONFIG for running in terminal
args <- commandArgs(trailingOnly = TRUE)

# Default values
main_folder <- "."
testing <- FALSE
sources <- c("AR", "BO", "BR", "CO", "MX", "US", "LE")
profile_data_suffix <- "smoothed.tsv"
protein <- NULL  
sd_multiplier_for_cutoff <- 4
only_proteins_above <- 0
only_proteins_below <- 0
output_suffix <- NULL

# User values
if (length(args) == 0) {
  cat("No arguments provided. Using default values.\n")
} else {
  # User values
  for (i in seq(1, length(args), by = 2)) {
    if (args[i] == "--main_folder") {
      main_folder <- args[i + 1]
    }
    if (args[i] == "--testing") {
      testing <- as.logical(args[i + 1])
    }
    if (args[i] == "--sources") {
      sources <- unlist(strsplit(args[i + 1], ","))
    }
    if (args[i] == "--profile_data_suffix") {
      profile_data_suffix <- args[i + 1]
      output_suffix <- args[i + 1]
    }
    if (args[i] == "--protein") {
      protein <- args[i + 1]
    }
    if (args[i] == "--sd_multiplier_for_cutoff") {
      sd_multiplier_for_cutoff <- as.numeric(args[i + 1])
    }
    if (args[i] == "--only_proteins_above") {
      only_proteins_above <- as.numeric(args[i + 1])
    }
    if (args[i] == "--only_proteins_below") {
      only_proteins_below <- as.numeric(args[i + 1])
    }
    if (args[i] == "--output_suffix") {
      output_suffix <- args[i + 1]
    }
  }
}

# Other plot parameters
geom_point_size <- 1.5
geom_line_size <- 0.5

title_size <- 16
axis_title_size <- 16
axis_text_size <- 16
axis_ticks_length <- 0.2
panel_background_size <- 1
plot_margin_l <- 0
plot_margin_r <- 0.5

geom_errorbar_width <- 1
geom_errorbar_size <- 0.5

legend_x <- 0.8
legend_y <- 0.9
legend_text_size <- 12
legend_element_width <- 2
legend_element_height <- 1.5


# Print values
cat("Main folder:", main_folder, "\n")
cat("Testing:", testing, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n")
cat("Profile data suffix:", profile_data_suffix, "\n")
cat("Protein (optional):", protein, "\n")
cat("SD multiplier for cutoff>", sd_multiplier_for_cutoff, "\n")

#### PATH CONFIG ####
if (testing == TRUE) {
  # For testing
  project_folder <- sprintf("%s/data/test_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
} else {
  # For running the actual data
  project_folder <- sprintf("%s/data/chagastope_data", main_folder)
  extra_files_folder <- sprintf("%s/data/extra_files", main_folder)
}

input_folder <- sprintf("%s/outputs/02_pools_smoothed_data", project_folder)
design_file <- sprintf("%s/inputs/01_pools_array_design/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv", project_folder)
output_folder <- sprintf("%s/outputs/08_plots", project_folder)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

functions_folder <- sprintf("%s/code/functions", main_folder)

#### CALL MAIN FUNCTION ####
source(sprintf("%s/08_plot_protein_profiles_function.R", functions_folder))

plot_proteins(project_folder = project_folder,
              input_folder = input_folder,
              design_file = design_file,
              sources = sources,
              profile_data_suffix = profile_data_suffix,
              output_folder = output_folder,
              protein = protein,
              only_proteins_above = only_proteins_above,
              only_proteins_below = only_proteins_below,
              sd_multiplier_for_cutoff = sd_multiplier_for_cutoff,
              output_suffix = output_suffix)

plot_parameters(geom_point_size_parameter =  geom_point_size, geom_line_size_parameter = geom_line_size,
                title_size_parameter = title_size, axis_title_size_parameter = axis_title_size, axis_text_size_parameter = axis_text_size, 
                axis_ticks_length_parameter = axis_ticks_length, panel_background_size_parameter = panel_background_size, 
                plot_margin_l_parameter = plot_margin_l, plot_margin_r_parameter = plot_margin_r,
                geom_errorbar_width_parameter = geom_errorbar_width, geom_errorbar_size_parameter = geom_errorbar_size,
                legend_x_parameter = legend_x, legend_y_parameter = legend_y, legend_text_size_parameter = legend_text_size, 
                legend_element_width_parameter = legend_element_width, legend_element_height_parameter = legend_element_height)

