############################-
#### Required libraries ####
############################-
renv::restore()

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggridges)
library(scico)
library(dplyr)
library(ggh4x)

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
    sprintf("individual_serums_protein_profiles_cruzi_%sSD", sd_multiplier_for_cutoff)
  } else {
    sprintf("individual_serums_protein_profiles_cruzi_%sSD_selected_proteins", sd_multiplier_for_cutoff)
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

#######################-
#### MAIN FUNCTION ####-
#######################
plot_proteins_ridge = function(project_folder, input_folder, design_file, sources, 
                               profile_data_suffix, output_folder, protein, 
                               sd_multiplier_for_cutoff, output_suffix, only_proteins_above, only_proteins_below, pdf_width, pdf_height) { 
  
  
  plot_data_result <- prepare_plot_data( project_folder = project_folder, input_folder = input_folder, 
                                         protein = protein, profile_data_suffix = profile_data_suffix, 
                                         sources = sources, output_folder = output_folder, 
                                         sd_multiplier_for_cutoff = sd_multiplier_for_cutoff, 
                                         only_proteins_above = only_proteins_above, 
                                         only_proteins_below = only_proteins_below, 
                                         output_suffix = output_suffix) 
  
  all_plot_data <- plot_data_result$all_plot_data
  all_plot_data$source <- factor(all_plot_data$source)
  all_plot_data$country <- as.factor(sub("_.*", "", all_plot_data$source))

  plot_order <- plot_data_result$plot_order

  proteins_to_plot <- plot_data_result$proteins_to_plot
  
  color_palette <- scico(6, palette = 'batlow')
  countries <- c("AR", "BO", "BR", "CO", "MX", "US")
  country_colors <- setNames(color_palette, countries)
  
  all_plot_data$country_color <- country_colors[all_plot_data$country]
  
  #Crear archivo PDF 
  output_file_prefix <- if (is.null(protein)) { 
    sprintf("individual_serums_protein_ridgeplots_%sSD", sd_multiplier_for_cutoff) 
    } else { 
      sprintf("individual_serums_protein_ridgeplots_%sSD_selected_proteins", sd_multiplier_for_cutoff) 
      } 
  output_file_name <- if (is.null(output_suffix)) { 
    sprintf("%s/%s.pdf", output_folder, output_file_prefix) 
    } else { sprintf("%s/%s_%s.pdf", output_folder, output_file_prefix, output_suffix) 
      } 
  
  pdf(file = output_file_name, width = pdf_width, height = pdf_height) 
  
  for (protein_for in proteins_to_plot) { 
    #protein_for <- proteins_to_plot[1]
    cat(sprintf("Plotting %s...\n", protein_for)) 
    sub_plot_data <- all_plot_data[protein == protein_for]
    
    unique_country_colors <- country_colors[names(country_colors) %in% unique(sub_plot_data$country)]
    strip <- strip_themed(background_y = elem_list_rect(fill = unique_country_colors))
    
    p <- ggplot(sub_plot_data, aes(x = start, y = source, height = mean_smoothed_signal, fill = mean_smoothed_signal)) + 
      geom_density_ridges_gradient(stat = "identity", scale = 3) + 
      scale_fill_scico(palette = "lipari", aesthetics = "fill") + #AGREGAR , alpha = 0.4
      theme_minimal() + 
      labs(title = sprintf("Signal plot for %s protein", protein_for), x = "Peptide position", 
           y = "Source", fill = "Signal") +
      scale_y_discrete(position = "right") + 
      facet_wrap2(~ country, ncol = 1, scales = "free_y", strip = strip, strip.position = "left")+
      theme(
        strip.text = element_text(color = "white"),  # Color del texto de la etiqueta
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        axis.ticks.y = element_line(color = "black"),
        axis.text.y = element_text(vjust = 0)
      )
    
    print(p) } 
  
  dev.off()
  graphics.off()
}
