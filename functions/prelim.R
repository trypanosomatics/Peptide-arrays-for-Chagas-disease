plot_data_result <- prepare_plot_data( project_folder = project_folder, input_folder = input_folder, 
                                       protein = protein, profile_data_suffix = profile_data_suffix, 
                                       sources = sources, output_folder = output_folder, 
                                       sd_multiplier_for_cutoff = sd_multiplier_for_cutoff, 
                                       only_proteins_above = only_proteins_above, 
                                       only_proteins_below = only_proteins_below, 
                                       output_suffix = output_suffix ) 

all_plot_data <- plot_data_result$all_plot_data
all_plot_data$country <- as.factor(sub("_.*", "", all_plot_data$source))

plot_order <- plot_data_result$plot_order

proteins_to_plot <- plot_data_result$proteins_to_plot

color_palette <- scico(6, palette = 'batlow')

#prueba todos los paises
ggplot(all_plot_data[protein=="TcCLB.511671.60",], aes(x = start, y = source, height = mean_smoothed_signal, fill = country)) + 
  geom_density_ridges(stat = "identity", scale = 2, alpha = 0.8) + 
  scale_fill_manual(values = color_palette) +
  scale_y_discrete(limits = rev(unique(sort(all_plot_data$source)))) +
  theme_minimal() + labs(title = "Ridge Plot of Signal Changes", x = "Start Position", 
                         y = "Source", fill = "Source")

protein="TcCLB.511671.60"

ggplot(all_plot_data[protein=="TcCLB.511671.60",], aes(x = start, y = source, height = mean_smoothed_signal, fill = country)) + 
  geom_density_ridges(stat = "identity", scale = 2, alpha = 0.8) + 
  scale_fill_manual(values = color_palette) + 
  theme_minimal() + 
  labs(title = sprintf("Signal plot for %s protein", protein), x = "Start Position", y = "Serum", fill = "Country") 

