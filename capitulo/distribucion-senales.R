setwd("~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/")
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tidyr)
library(patchwork)


#### Raw data ####

path_to_files <- "~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/inputs/02_pools_raw_data"
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

path_to_files <- "~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/outputs/01_pools_normalized_data/"
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

path_to_files <- "~/Documentos/Guada/Array/Peptide-arrays-for-Chagas-disease/test_data/outputs/02_pools_smoothed_data//"
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

PO_raw <- ggplot(all_data_raw[all_data_raw$type=="PO",], aes(x = source, y = Signal)) +
  geom_jitter(aes(colour = source), alpha = 0.3) +
  #geom_violin(alpha = 0.6) +
  geom_boxplot()+
  scale_color_manual(values = palette_neutral) +
  labs(title = "Distribution of Signal Intensities - Raw data",
       x = "Source",
       y = "Signal Intensity") +
  ylim(c(0, 70000)) +
  theme_minimal() +
  theme(legend.position = "none")

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
