#### PREPARATION STEPS ####
## You can run this code as it is and reproduce Ag2 (TcCLB.511671.60) Alanine Scan results, or you can follow these next steps for the entire dataset.

## 1. In chagastope_data/inputs/11_individual_serums_array_design place the "Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv" file (links in the paper)
## 2. In chagastope_data/inputs/12_individual_serums_raw_data place each of the raw data files from CHAGASTOPE-v2 (such as AR_P1_PO_raw.tsv) downloaded from Array Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11655/)
## 3. Set the "testing" variable in the config below to FALSE, or run this code with the "-test F" argument
## 4. If you are running this code in Rstudio, set the "main_folder" variable in the config below to the folder containing this code

#### WARNINGS ####
## This code uses large amounts of RAM

#### CONFIG ####
main_folder <- "." #When running in Rstudio, set this to the absolute path of the folder containing this code
testing <- TRUE #set this to FALSE when running the actual data

#### READ ARGUMENTS AND GET PATH (DO NOT CHANGE) ####
args <- commandArgs(TRUE)

if (length(args == 2)) {
    if (args[1] == "-test") {
        testing <- as.logical(args[2])
    }
}

if (testing == TRUE) {
    #For testing
    project_folder <- sprintf("%s/test_data", main_folder) 
    design_alanine_file <- sprintf("%s/77_extra_files/alanine_scan_design_test_subset.tsv", main_folder)  
} else {
    #For running the actual data
    project_folder <- sprintf("%s/chagastope_data", main_folder)
    design_alanine_file <- sprintf("%s/77_extra_files/alanine_scan_design.tsv", main_folder)
}

#### INTERNAL CONFIG (DO NOT CHANGE) ####
library(dplyr)
library(reshape2)
library(data.table)

design_data_file <- sprintf("%s/inputs/11_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv", project_folder)

#### FUNCTIONS ###
### estimateSignalChanges ###
## This function estimates the signal change per residue and for each sample analyzed in one protein.
estimateSignalChanges = function(selected_protein) {
  design_alanine_filtered <- filter(design_alanine, protein == selected_protein)
  AlanineScanData <- merge(Chips_roche_norm,design_alanine_filtered,by.x = "Reporter.Name",by.y = "array_express_id") #This step merges data from design and CHIPS results.
  AlanineScanData$signal <- (AlanineScanData$allData.replica1 + AlanineScanData$allData.replica2) / 2 #Signal is estimated as mean of both replicas.
  AlanineScanData$allData.replica1 <- NULL
  AlanineScanData$allData.replica2 <- NULL
  orders <- group_by(filter(AlanineScanData,alanine_position == 0),sample) %>% summarise(signal = mean(signal))
  samples <- orders[order(orders$signal,decreasing = T),]$sample
  
  
  AlanineScanData$alanine_position_start <- AlanineScanData$peptide_to_scan_start + AlanineScanData$alanine_position - 1
  AlanineScanData$alanine_position_start[AlanineScanData$alanine_position == 0] <- 0
  
  temp_protein <- filter(AlanineScanData,!is.na(signal))
  if (length(temp_protein$Reporter.Name) < 2) {next()}
  temp_protein$change <- 0
  temp_originals <- select(filter(temp_protein,alanine_position_start == 0),peptide_to_scan,signal,sample) 
  temp_originals$original_signal <- temp_originals$signal
  temp_originals$signal <- NULL
  changeDT_total <- data.frame(change = 0,sample = "",stringsAsFactors = F,mutated_aa = "",protein = "")
  changeDT_total <- changeDT_total[0,]
  for (s in samples) {
    changeDT = filter(temp_protein,sample == s)
    temp_originals_filtered <- filter(temp_originals, sample == s)
    changeDT <- merge(changeDT,temp_originals_filtered,by = c("peptide_to_scan","sample"))
    changeDT$change <- changeDT$signal - changeDT$original_signal #Here we estimate signal change as the difference between signal of mutated peptides and the ones that has the original sequence of the protein.
    
    mutated_aa <- changeDT$alanine_position_start
    mutated_aa[changeDT$alanine_position == 0] <- 0
    mutated_aa <- unique(mutated_aa)
    
    changeDT$mutated_aa <- substr(changeDT$peptide_to_scan,changeDT$alanine_position,changeDT$alanine_position) #Here we save the aa that we are mutating.
    changeDT_temp = dplyr::group_by(filter(changeDT,alanine_position != 0),alanine_position_start,sample) %>% dplyr::summarise(change = mean(change),
                                                                                                                               mutated_aa = unique(mutated_aa),
                                                                                                                               protein = unique(protein))
    changeDT_total <- rbind(changeDT_total,as.data.frame(changeDT_temp))
  }
  return(changeDT_total)
} 
### generateHeatmapMatrix ###
## This function transform data frame into a matrix suitable for heatmap visualization. Samples are separated into rows and analyzed residues are in columns.
generateHeatmapMatrix = function(dt){
  samples <- unique(dt$sample)
  temp_matrix <- acast(dt, paste(alanine_position_start,mutated_aa)~sample, value.var = "change",fun.aggregate = mean,na.rm = T)
  temp_matrix <- temp_matrix[order(as.numeric(gsub("[^0-9.]", "",rownames(temp_matrix)))),,drop = F]
  temp_matrix <- t(temp_matrix)
  temp_matrix <- temp_matrix[samples,] #this is only for ordering samples by reactivity
  return(temp_matrix)
}
## This functions takes the matrix generated, clusterize samples by similarity and renders a heatmap visualization with the same color scales used in the manuscript. Only used for ploting, not needed elsewhere.
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
make_heatmap = function(temp_matrix,selected_protein,selected_serums = "") {
  colours_options <- c("red","blue","orange","green","white","pink","darkblue","grey")
  colours_cb <- data.frame(names = colours_options,cb = c(palette_OkabeIto[c(6,2,1,3)],"white",palette_OkabeIto[c(7,5,8)]))
  temp_matrix <- temp_matrix[,as.numeric(as.character(lapply(strsplit(colnames(temp_matrix)," "), `[[`, 1))) > 0]
  paletteLength <- 250
  myColor <- c(colorRampPalette(c(colours_cb$cb[colours_cb$names == "red"], "white"))(paletteLength/2),colorRampPalette(c("white",colours_cb$cb[colours_cb$names == "blue"]))(paletteLength/2))
  myBreaks <- c(seq(min(temp_matrix), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(temp_matrix)/paletteLength, max(temp_matrix), length.out = floor(paletteLength/2)))
  if (sum(abs(colMeans(temp_matrix)) > 1300) == 0) { 
    dist <- dist(temp_matrix)} else {
      dist <- dist(temp_matrix[,abs(colMeans(temp_matrix)) > 1300])
    }
  hc = hclust(dist , method = "average")
  if (selected_serums != "") {
    heatmap <- pheatmap::pheatmap(temp_matrix,cluster_cols = F,cluster_rows = hc,color = myColor, breaks = myBreaks,border_color = NA,main = selected_protein,
                                  labels_row = make_bold_names(temp_matrix, rownames, selected_serums))
    colors_heatmap = ifelse(grepl("bold",heatmap$gtable$grobs[[5]]$label),colours_cb$cb[6],"black")
    heatmap$gtable$grobs[[5]]$gp = gpar(col = colors_heatmap)
  } else {
    heatmap <- pheatmap::pheatmap(temp_matrix,cluster_cols = F,cluster_rows = hc,color = myColor, breaks = myBreaks,border_color = NA,main = selected_protein)
  }
  return(heatmap)
} 

# setwd(project_folder)
design_alanine <- read.csv(design_alanine_file, sep = "\t")

input_data_path <- sprintf("%s/inputs/12_individual_serums_raw_data/", project_folder)
files <- list.files(input_data_path)
files <- files[grepl("raw.tsv",files)]

#Get only the ids for the AlaScan
design_data <- fread(design_data_file, header = T, sep = "\t", na.strings = NULL)
design_data <- design_data[group == "alascan_best_peaks_Tcruzi"]
ids_to_keep <- unique(design_data$array_express_id)

for (i in 1:length(files)) {
  # i <- 1
  Chips_roche_norm_dat <- read.csv(paste0(input_data_path,files[i],sep = ""),sep = "\t",stringsAsFactors = F)
  Chips_roche_norm_dat$sample <- paste(strsplit(files[i],"_")[[1]][c(1:2)],collapse = "_")
  
  #Keep only the raw data from the AlaScan
  Chips_roche_norm_dat <- Chips_roche_norm_dat[Chips_roche_norm_dat$Reporter.Name %in% ids_to_keep,]

  if (i == 1) {
    Chips_roche_norm <- Chips_roche_norm_dat
    next()
  }
  Chips_roche_norm <- rbind(Chips_roche_norm,Chips_roche_norm_dat)
} #Load all micro-array results for all samples


for (protein_to_estimate in unique(design_alanine$protein)) {
  dt <- estimateSignalChanges(protein_to_estimate)
  matrix <- generateHeatmapMatrix(dt)  
  #### SAVE RESULTS ####
  write.table(dt,paste0(project_folder, "/outputs/21_alanine_scan_raw_data/Raw_data_signal_change_alanine_scan_", protein_to_estimate,".tsv"),sep = "\t")  #long format
  write.table(matrix,paste0(project_folder, "/outputs/21_alanine_scan_raw_data/Raw_data_signal_change_matrix_alanine_scan_", protein_to_estimate,".tsv"),sep = "\t") #matrix for heatmap
  
  # OPTIONAL: using "pheatmap" package you could visualize results as in the manuscript.
  # pdf(file = paste0(project_folder, "/outputs/21_alanine_scan_raw_data/Heatmap_alanine_scan_", protein_to_estimate,".pdf"),height = 12, width = 12)
  # print(make_heatmap(temp_matrix = matrix,selected_protein = protein_to_estimate))
  # dev.off()
}
