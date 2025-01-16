setwd("~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/")
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tidyr)
library(patchwork)
library(scico)
library(data.table)

#### Raw data ####

path_to_files <- "C://Users/romer/Documentos/Peptide-arrays-for-Chagas-disease/Peptide-arrays-for-Chagas-disease/test_data/inputs/02_pools_raw_data/"
files <- list.files(path = path_to_files, pattern = "raw\\.tsv$", full.names = TRUE)
all_data_raw <- data.frame()

for (file in files) {
  # file <- files[1]
  data <- read_delim(file, delim = "\t", trim_ws = TRUE)
  
  data_long <- data %>%
    pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")
  all_data_raw <- bind_rows(all_data_raw, data_long)
}

all_data_raw$source <- sapply(strsplit(all_data_raw[,2], "_"), `[`, 1)
all_data_raw$type <- sapply(strsplit(all_data_raw[,2], "_"), `[`, 2)
all_data_raw$type <- as.factor(all_data_raw$type)
all_data_raw$source <- as.factor(all_data_raw$source)

rm(data, data_long)

#### Normalized data ####

#path_to_files <- "~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/outputs/01_pools_normalized_data/"
path_to_files <- "C://Users/romer/Documentos/Peptide-arrays-for-Chagas-disease/Peptide-arrays-for-Chagas-disease/test_data/outputs/01_pools_normalized_data/"
files <- list.files(path = path_to_files, pattern = "processed\\.tsv$", full.names = TRUE)
all_data_normalized <- data.frame()

for (file in files) {
  # file <- files[1]
  data <- read_delim(file, delim = "\t", trim_ws = TRUE)
  
  data_long <- data %>%
    pivot_longer(cols = -`Reporter Name`, names_to = "Sample", values_to = "Signal")

  all_data_normalized <- bind_rows(all_data_normalized, data_long)
}

all_data_normalized$source <- sapply(strsplit(all_data_normalized[,2], "_"), `[`, 1)
all_data_normalized$type <- sapply(strsplit(all_data_normalized[,2], "_"), `[`, 2)
all_data_normalized$type <- as.factor(all_data_normalized$type)
all_data_normalized$source <- as.factor(all_data_normalized$source)

rm(data, data_long)

#### Smoothed data ####

#path_to_files <- "~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/outputs/02_pools_smoothed_data//"
path_to_files <- "C://Users/romer/Documentos/Peptide-arrays-for-Chagas-disease/Peptide-arrays-for-Chagas-disease/test_data/outputs/02_pools_smoothed_data//"
files <- list.files(path = path_to_files, pattern = "smoothed\\.tsv$", full.names = TRUE)
all_data_smoothed <- data.frame()

for (file in files) {
  # file <- files[1]
  data <- read_delim(file, delim = "\t", trim_ws = TRUE)

  all_data_smoothed <- bind_rows(all_data_smoothed, data)
}

all_data_smoothed$type <- as.factor(all_data_smoothed$type)
all_data_smoothed$source <- as.factor(all_data_smoothed$source)
rm(data, file, files, path_to_files)

#### Plot general, toda las sources juntas ####
palette_neutral <- c("#588b8b", "#ef4043", "#ffd5c2", "#f28f3b", "#c8553d", "#2d3047", "#93b7be")

raw_plot <- ggplot(all_data_raw, aes(x = type, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.6)+
  geom_violin(alpha = 0.6)+
  scale_color_manual(values = palette_neutral)+
  #geom_boxplot(alpha = 0.6)+
  labs(title = "Distribution of Signal Intensities - Raw data",
       x = "Signal Intensity",
       y = "Frequency") +
  ylim(c(0,70000))+
  theme_minimal()+
  theme(legend.position = "none")
raw_plot

normalized_plot <- ggplot(all_data_normalized, aes(x = type, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.6)+
  scale_color_manual(values = palette_neutral)+
  geom_violin(alpha = 0.6)+
  #geom_boxplot(alpha = 0.6)+
  labs(title = "Distribution of Signal Intensities - Normalized data",
       x = "Signal Intensity",
       y = "Frequency") +
  ylim(c(0,70000))+
  theme_minimal()+
  theme(legend.position = "none")
normalized_plot

smoothed_plot <- ggplot(all_data_smoothed, aes(x = type, y = mean_smoothed_signal)) +
  geom_jitter(aes(colour = source), alpha = 0.5)+
  scale_color_manual(values = palette_neutral)+
  geom_violin(alpha = 0.6)+
  #geom_boxplot(alpha = 0.6)+
  labs(title = "Distribution of Signal Intensities - Smoothed data",
       x = "Signal Intensity",
       y = "Frequency") +
  ylim(c(0,70000))+
  theme_minimal()
smoothed_plot

raw_plot + normalized_plot + smoothed_plot + plot_layout(ncol = 3)

#### Plot sources por separado ####

NE_raw <- ggplot(all_data_raw[all_data_raw$type=="NE",], aes(x = source, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Raw data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 40000)) +
  theme_minimal() +
  theme(legend.position = "none")

PO_raw <- ggplot(all_data_raw[(all_data_raw$type=="PO" & all_data_raw$source != "LE"),], aes(x = source, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Raw data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 70000)) +
  theme_minimal() +
  theme(legend.position = "none")+
  geom_bin_2d(binwidth = c(0.1, 0.1))

NE_norm <- ggplot(all_data_normalized[all_data_normalized$type=="NE",], aes(x = source, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Normalized data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 40000)) +
  theme_minimal() +
  theme(legend.position = "none")

PO_norm <- ggplot(all_data_normalized[all_data_normalized$type=="PO",], aes(x = source, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Normalized data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 70000)) +
  theme_minimal() +
  theme(legend.position = "none")

NE_smooth <- ggplot(all_data_smoothed[all_data_smoothed$type=="NE",], aes(x = source, y = mean_smoothed_signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Smoothed data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 40000)) +
  theme_minimal()

PO_smooth <- ggplot(all_data_smoothed[all_data_smoothed$type=="PO",], aes(x = source, y = mean_smoothed_signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Smoothed data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 70000)) +
  theme_minimal()


NE_raw + NE_norm + NE_smooth + plot_layout(ncol = 3)
PO_raw + PO_norm + PO_smooth + plot_layout(ncol = 3)


#### Plot con densidad de puntos ----------------------------------------

##### Funcion #####

DensityScatterplot <- function(plot_data, x_column_name, y_column_name, gradient_palette = "hawaii",
                                                  bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                                                  geom_point_size = 1, geom_point_shape = 16, plot_title = "", x_label = "X", y_label = "Y",
                                                  axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                                                  fixed_scale = 0, fixed_scale_min = 0, fixed_scale_max = 0, y_lim_min = NULL, y_lim_max = NULL) {
  
  plot_data_aux <- as.data.table(plot_data)
  setnames(plot_data_aux, x_column_name, "x_column")
  setnames(plot_data_aux, y_column_name, "y_column")
  
  # Determina una malla para calcular la densidad de puntos en cada celda
  min_y_value <- floor(min(plot_data_aux$y_column))
  max_y_value <- ceiling(max(plot_data_aux$y_column))
  grid_size_y <- (max_y_value - min_y_value + 1) / bin_amount_per_axis
  plot_data_aux[, y_cell := ((y_column - min_y_value) %/% grid_size_y) + 1]
  
  # Calcula la cantidad de puntos por celda (N)
  point_density <- plot_data_aux[, .N, by = .(y_cell)]
  
  # Une los datos de densidad con el dataframe original
  plot_data_aux <- merge(plot_data_aux, point_density[, .(y_cell, N)], by = "y_cell")
  
  # Genera el gráfico con los puntos en el eje X y la densidad en Y
  p <- ggplot(plot_data_aux[order(plot_data_aux$y_column),]) + 
    geom_point(aes(x = x_column, y = y_column, col = log10(N+1)), size = geom_point_size, shape = geom_point_shape) +  # Usa N para el mapeo del color
    #scale_color_gradient(low = "#000099", high = "#FF3100", guide = "colorbar") +  # Escala de colores para densidad
    scale_color_scico(palette = gradient_palette, guide = "colorbar")+
    theme_bw() +
    theme(axis.title = element_text(size = axis_title_size, colour = "black"),
          axis.text = element_text(size = axis_text_size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(axis_ticks_length, "cm"),
          panel.background = element_rect(size = panel_background_size, colour = "black"),
          aspect.ratio = 2/1)
  
  # Personalización del gráfico
  if (plot_title != "") {
    p <- p + ggtitle(label = plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (x_label != "") {
    p <- p + xlab(x_label)
  }
  if (y_label != "") {
    p <- p + ylab(y_label)
  }
  
  # Agregar límite en Y si se especifica
  if (!is.null(y_lim_min) & !is.null(y_lim_max)) {
    p <- p + ylim(y_lim_min, y_lim_max)
  }
  
  # Agregar escala fija en caso de estar activada
  if (fixed_scale) {
    p <- p + coord_cartesian(ylim = c(fixed_scale_min, fixed_scale_max))
  }
  
  # Mostrar gráfico
  p
}




# Jitter coloreado por densidad -------------------------------------------
# Calcular densidad en celdas
library(scales)

# Crear una transformación personalizada
custom_trans <- trans_new(
  name = "custom",
  transform = function(x) ifelse(x <= 10000, x * 2, 15000 + (x - 10000) / 3),
  inverse = function(y) ifelse(y <= 10000, y / 2, 10000 + (y - 15000) * 3)
)

# Función
DensityJitterplot <- function(plot_data, x_column_name, y_column_name, gradient_palette = "hawaii",
                              bin_amount_per_axis = 1000, max_data_density_to_consider_per_cell = 100000, gradient_steps = 200,
                              geom_point_size = 1, geom_point_shape = 16, plot_title = "", x_label = "X", y_label = "Y",
                              axis_title_size = 16, axis_text_size = 16, axis_ticks_length = 0.2, panel_background_size = 1,
                              fixed_scale = 0, fixed_scale_min = 0, fixed_scale_max = 0, y_lim_min = NULL, y_lim_max = NULL, legend_position = "right") {
  
  plot_data_aux <- as.data.table(plot_data)
  setnames(plot_data_aux, x_column_name, "x_column")
  setnames(plot_data_aux, y_column_name, "y_column")
  
  # Ajustar rango de densidad
  min_y_value <- floor(min(plot_data_aux$y_column))
  max_y_value <- ceiling(max(plot_data_aux$y_column))
  grid_size_y <- (max_y_value - min_y_value + 1) / bin_amount_per_axis
  plot_data_aux[, y_cell := ((y_column - min_y_value) %/% grid_size_y) + 1]
  
  # Calcular densidad de puntos por celda
  point_density <- plot_data_aux[, .N, by = .(y_cell)]
  plot_data_aux <- merge(plot_data_aux, point_density[, .(y_cell, N)], by = "y_cell")
  
  # Genera el gráfico con los puntos en el eje X y la densidad en Y
  p <- ggplot(plot_data_aux[order(plot_data_aux$y_column),]) +
    geom_jitter(aes(x = x_column, y = y_column, colour = log10(N + 1)), 
                width = 0.2, size = geom_point_size, shape = geom_point_shape) +
    scale_color_scico(palette = gradient_palette, guide = "colorbar") +
    scale_y_continuous(trans = custom_trans,
                       breaks = c(0, 2500, 5000, 7500, 10000, 30000, 40000, 50000, 60000, 70000),
                       labels = c("0", "2500", "5000", "7500", "10000", "30000", "40000", "50000",  "60000", "70000")) +
    theme_bw() +
    theme(axis.title = element_text(size = axis_title_size, colour = "black"),
          axis.text = element_text(size = axis_text_size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(axis_ticks_length, "cm"),
          panel.background = element_rect(size = panel_background_size, colour = "black"),
          aspect.ratio = 2/1, legend.position = legend_position)
  
  # Personalización del gráfico
  if (plot_title != "") {
    p <- p + ggtitle(label = plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (x_label != "") {
    p <- p + xlab(x_label)
  }
  if (y_label != "") {
    p <- p + ylab(y_label)
  }
  
  # Agregar límite en Y si se especifica
  if (!is.null(y_lim_min) & !is.null(y_lim_max)) {
    p <- p + coord_cartesian(ylim = c(y_lim_min, y_lim_max))
  }
  
  # Mostrar gráfico
  p
}

#### Datasets ####
NE_raw <- all_data_raw[all_data_raw$type=="NE",]
PO_raw <- all_data_raw[(all_data_raw$type=="PO" & all_data_raw$source!="LE"),]
NE_norm <- all_data_normalized[all_data_normalized$type=="NE",]
PO_norm <- all_data_normalized[(all_data_normalized$type=="PO" & all_data_normalized$source!="LE"),]
NE_smooth <- all_data_smoothed[all_data_smoothed$type=="NE",]
PO_smooth <- all_data_smoothed[(all_data_smoothed$type=="PO" & all_data_smoothed$source!="LE"),]

#### Plots ####
max(all_data_raw[all_data_raw$type=="PO", "Signal"]) #y_max PO = 66000
max(all_data_raw[all_data_raw$type=="NE", "Signal"]) #y_max NE = 37000

##### Celdas #####
NE_raw_plot <- DensityScatterplot(NE_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 3, 
                   geom_point_shape = 15, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                   y_lim_min = 0, y_lim_max = 37000)

NE_raw_plot <- DensityJitterplot(NE_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 3, 
                                  geom_point_shape = 20, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                  y_lim_min = 0, y_lim_max = 37000, legend_position = "none")

NE_norm_plot <- DensityScatterplot(NE_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 3, 
                                  geom_point_shape = 15, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                  y_lim_min = 0, y_lim_max = 37000)

NE_smooth_plot <- DensityScatterplot(NE_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 3, 
                                  geom_point_shape = 15, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                  y_lim_min = 0, y_lim_max = 37000)

NE_raw_plot + NE_norm_plot + NE_smooth_plot + plot_layout(ncol = 3)

PO_raw_plot <- DensityScatterplot(PO_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 3, 
                                  geom_point_shape = 15, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                  y_lim_min = 0, y_lim_max = 66000)
PO_norm_plot <- DensityScatterplot(PO_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 3, 
                                   geom_point_shape = 15, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                   y_lim_min = 0, y_lim_max = 66000)
PO_smooth_plot <- DensityScatterplot(PO_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 3, 
                                     geom_point_shape = 15, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                     y_lim_min = 0, y_lim_max = 66000)

PO_raw_plot + PO_norm_plot + PO_smooth_plot + plot_layout(ncol = 3)


##### Puntos #####

#Eliminar ejes
woaxis <- theme(
  axis.title.y = element_blank(),  # Remove y-axis title
  axis.text.y = element_blank(),   # Remove y-axis text
  axis.ticks.y = element_blank(),  # Remove y-axis ticks
)

NE_raw_plot <- DensityJitterplot(NE_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                 geom_point_shape = 16, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                 y_lim_min = 0, y_lim_max = 37000, legend_position = "none")

NE_norm_plot <- DensityJitterplot(NE_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                   geom_point_shape = 16, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                   y_lim_min = 0, y_lim_max = 37000, legend_position = "none")

NE_smooth_plot <- DensityJitterplot(NE_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 2, 
                                     geom_point_shape = 16, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                     y_lim_min = 0, y_lim_max = 37000)


NE_raw_plot + (NE_norm_plot + woaxis) + (NE_smooth_plot + woaxis) + plot_annotation(title = "Negatives")

PO_raw_plot <- DensityJitterplot(PO_raw, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                  geom_point_shape = 16, plot_title = "Raw", x_label = "Source", y_label = "Signal",
                                  y_lim_min = 0, y_lim_max = 66000, legend_position = "none")
PO_norm_plot <- DensityJitterplot(PO_norm, "source", "Signal", gradient_palette = "lipari", geom_point_size = 2, 
                                   geom_point_shape = 16, plot_title = "Normalized", x_label = "Source", y_label = "Signal",
                                   y_lim_min = 0, y_lim_max = 66000, legend_position = "none")
PO_smooth_plot <- DensityJitterplot(PO_smooth, "source", "mean_smoothed_signal", gradient_palette = "lipari", geom_point_size = 2, 
                                     geom_point_shape = 16, plot_title = "Smoothed", x_label = "Source", y_label = "Signal",
                                     y_lim_min = 0, y_lim_max = 66000)

PO_raw_plot + (PO_norm_plot + woaxis) + (PO_smooth_plot + woaxis) + plot_annotation(title = "Positives")
