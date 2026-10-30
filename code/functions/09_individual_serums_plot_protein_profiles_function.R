############################-
#### Required libraries ####
############################-
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
    array_config,
    input_folder, 
    smoothing, 
    sources, 
    protein, 
    index_output_file,
    only_proteins_above, 
    only_proteins_below
) {
  
  # LOAD THE DATA NECESSARY
  input_files <- find_output_files(input_folder, array_config$id, "smoothed", source = sources, smoothing = smoothing)
  
  if (length(input_files) == 0) {
    stop(sprintf("No %s smoothed files found in %s for smoothing '%s' and the given sources.", array_config$id, input_folder, smoothing))
  }
  
  # Get data
  all_plot_data <- rbindlist(lapply(input_files, function(input_file) {
    input_data <- fread(input_file, header = TRUE, sep = "\t", na.strings = NULL)
    check_array(input_data, array_config$id, input_file)
    input_data
  }))
  all_plot_data[, array := NULL]
  
  # Filter by protein
  if (!is.null(protein)) {
    selected_protein <- protein # data.table would read "protein" as the column
    all_plot_data <- all_plot_data[all_plot_data$protein == selected_protein]
    if (nrow(all_plot_data) == 0) {
      stop(paste("No data found for the specified protein:", protein))
    }
  }
  
  # Order
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
  index_data <- add_array_column(index_data, array_config$id)
  
  write.table(index_data, 
              file = index_output_file, 
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              quote = TRUE)
  
  
  return(list(all_plot_data = all_plot_data, plot_order = plot_order, proteins_to_plot = proteins_to_plot))
}  

#######################-
#### MAIN FUNCTION ####-
#######################
plot_proteins_ridge = function(array_config, input_folder, sources, 
                               smoothing, output_folder, protein, 
                               only_proteins_above, only_proteins_below, pdf_width, pdf_height) { 
  
  scope <- if (is.null(protein)) "all" else protein
  
  plot_data_result <- prepare_plot_data(array_config = array_config, input_folder = input_folder, 
                                        smoothing = smoothing, sources = sources, protein = protein, 
                                        index_output_file = output_file_name(output_folder, array_config$id, "ridgeplots_index", smoothing = smoothing, scope = scope), 
                                        only_proteins_above = only_proteins_above, 
                                        only_proteins_below = only_proteins_below) 
  
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
  pdf_output_file <- output_file_name(output_folder, array_config$id, "ridgeplots", smoothing = smoothing, scope = scope, ext = "pdf")
  
  pdf(file = pdf_output_file, width = pdf_width, height = pdf_height) 
  
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
