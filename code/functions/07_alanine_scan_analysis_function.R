#### LIBRARIES ####
if (!require(data.table, quietly = TRUE)) {
  writeLines("Installing library 'data.table' for R")
  install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(data.table)
}

if (!require(dplyr, quietly = TRUE)) {
  writeLines("Installing library 'dplyr' for R")
  install.packages("dplyr", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(dplyr)
}

if (!require(reshape2, quietly = TRUE)) {
  writeLines("Installing library 'reshape2' for R")
  install.packages("reshape2", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(reshape2)
}

# if (!require(gtable, quietly = TRUE)) {
#   writeLines("Installing library 'gtable' for R")
#   install.packages("gtable", repos = "http://cran.rstudio.com/", dependencies = TRUE)
#   library(gtable) #for logo, you need version 4.2.3
# }

if (!require(pheatmap, quietly = TRUE)) {
  writeLines("Installing library 'pheatmap' for R")
  install.packages("pheatmap", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(pheatmap)
}

if (!require(ggplot2, quietly = TRUE)) {
  writeLines("Installing library 'ggplot2' for R")
  install.packages("ggplot2", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(ggplot2) #for logo
}

if (!require(ggseqlogo, quietly = TRUE)) {
  writeLines("Installing library 'ggseqlogo' for R")
  install.packages("ggseqlogo", repos = "http://cran.rstudio.com/", dependencies = TRUE)
  library(ggseqlogo) #for logo, you need version 4.2.3
}

#### FUNCTION ####
alanine_scan <- function(main_folder, testing, selected_protein, heatmap_plot, sequence_logo_plot, sequence_logo_source) {
  
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
    temp_matrix <- temp_matrix[samples, , drop = FALSE] #this is only for ordering samples by reactivity
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
  
  splitData <- function(changeDT_total) {
    setDT(changeDT_total)
    group_mapping <- unique(changeDT_total[, .(alanine_position_start)])
    group_mapping[, group := cumsum(c(1, diff(alanine_position_start) != 1))]
    changeDT_total <- merge(changeDT_total, group_mapping, by = "alanine_position_start", all.x = TRUE)
    splitted_data <- split(changeDT_total, by = "group")
    return(splitted_data)
  }

  plotLogo <- function(sub_signal_change_per_position,
                       protein_id = selected_protein, sources = sequence_logo_source,
                       relativize_logo_columns_to_most_changing_aa = T,
                       logo_default_color = "#000000", logo_default_label = "Not Fundamental",
                       logo_fundamental_color = "#D55E00", logo_fundamental_residues_label = "Fundamental", logo_fundamental_residues_min_signal = 1300,
                       logo_suboptimal_color = "#56B4E9", logo_suboptimal_residues_label = "SubOptimal", logo_suboptimal_residues_max_signal = -1300,
                       plot_title = "",
                       show_legend = F) {
    #sub_signal_change_per_position has to be already filtered.
    #This means having a limited amount of AA, usually a single peak in a single protein
    #This data can have one or more serum samples
    
    #These filters are optional if you want to filter the data inside this function
    if (protein_id != "") {
      sub_signal_change_per_position <- sub_signal_change_per_position[sub_signal_change_per_position$protein == protein_id]
    }
    
    if (sources != "average") {
      sub_signal_change_per_position <- sub_signal_change_per_position[sample %in% sources]
    }
    
    #Calculate the average change per mutation in each position
    #I put -mean because the letter shown in the plot is the one that was lost and so the signal was reduced, so it makes more sense to see it like this
    
    mean_change_per_replacement <- sub_signal_change_per_position[, .(
      mean_signal_change = -mean(change)
    ), by = .(alanine_position_start, mutated_aa, protein)]
    
    mean_change_per_replacement <- dcast(mean_change_per_replacement, mutated_aa~alanine_position_start, value.var = "mean_signal_change")
    
    mean_change_per_replacement <- as.data.table(mean_change_per_replacement)
    
    if (relativize_logo_columns_to_most_changing_aa) {
      #This normalization makes it so the height of all letters in the plot reaches the height of the largest change (instead of stacking one above the other)
      max_value <- max(mean_change_per_replacement[, -1], na.rm = TRUE)
      sum_value <- sum(mean_change_per_replacement[, -1], na.rm = TRUE)
      mean_change_per_replacement[, (names(mean_change_per_replacement)[!names(mean_change_per_replacement) %in% 'mutated_aa']) := 
                                    lapply(.SD, function(x) x * max_value / sum_value), .SDcols = !'mutated_aa']
      }
    
    #Plot logo
    #Define the colors of each amino acid based on their signal
    position_total_change <- colSums(mean_change_per_replacement[, -c("mutated_aa")], na.rm = T)
    
    is_residue_fundamental <- position_total_change > logo_fundamental_residues_min_signal
    is_residue_suboptimal <- position_total_change < logo_suboptimal_residues_max_signal
    
    dt_colors_aux <- data.table(position = 1:length(is_residue_fundamental),
                                color = logo_default_color,
                                group = logo_default_label)
    dt_colors_aux[is_residue_fundamental, color := logo_fundamental_color]
    dt_colors_aux[is_residue_fundamental, group := logo_fundamental_residues_label]
    dt_colors_aux[is_residue_suboptimal, color := logo_suboptimal_color]
    dt_colors_aux[is_residue_suboptimal, group := logo_suboptimal_residues_label]
    
    
    mt_aux <- as.matrix(mean_change_per_replacement[, -c("mutated_aa")])
    rownames(mt_aux) <- mean_change_per_replacement$mutated_aa
    
    p <- ggplot() + 
      geom_hline(yintercept = logo_fundamental_residues_min_signal, linetype = "dashed", color = "black") +
      geom_logo(mt_aux, method = "custom", col_scheme = "position", position_colors = dt_colors_aux) +
      ggtitle(plot_title) +
      scale_x_continuous(minor_breaks = NULL , breaks =  1:dim(mt_aux)[2], labels = colnames(mt_aux)) +
      scale_y_continuous(minor_breaks = NULL) +
      xlab("Residue Position") +
      ylab("Relative Signal Change per Position") +
      # theme_pubclean() +
      theme(text = element_text(size = 15),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    if (show_legend == F) {
      p <- p + theme(legend.position = "none")
    }
    
    p
  }
  
  
  
  #############################-
  #### ggseqlogo Overrides ####
  #############################-
  #The following functions are extracted from ggseqlogo. More exactly, they are the ones modified by Leonel
  #They exist in a custom package, but I used this method so I can still have the "standard" ggseqlogo in other scripts
  #It is very likely that I don't need to import ALL these functions, but finding out which ones I actually need takes time I don't have
  
  #### col_schemes.R ####
  #' List color schemes available in ggseqlogo
  #' 
  #' @param v If true, font names are printed to stderr. Otherwise, color scheme names are returned as a character vector
  #' @export
  list_col_schemes <- function(v=T){
    
    col_schemes = c('auto', 'chemistry', 'chemistry2','hydrophobicity', 'nucleotide', 'nucleotide2',
                    'base_pairing', 'clustalx', 'taylor')
    if(!v) return(col_schemes)
    message('Available ggseqlogo color schemes:')
    for(f in col_schemes) message('\t', f)
  }
  
  
  # Get color scheme
  # @param col_scheme name of color scheme
  # @param seq_type sequence type of color scheme
  get_col_scheme <- function(col_scheme, seq_type='auto'){
    
    if (col_scheme == "position") {
      cs <- data.frame(letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
                       group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
                       stringsAsFactors = F)
      return(cs)
    }
    
    # Check if user-defined color scheme
    if(is.data.frame(col_scheme)){
      if(!'ggseqlogo_cs' %in% class(col_scheme)) 
        stop('Colour scheme must be generated using "make_col_scheme" function')
      return(col_scheme)
    }
    
    # Get ambigious colour scheme
    col_scheme = match.arg(col_scheme, list_col_schemes(F))
    
    # Get default color scheme for sequence type
    if(col_scheme == 'auto'){
      if(seq_type == 'auto') stop('"col_scheme" and "seq_type" cannot both be "auto"')
      
      col_scheme = switch(tolower(seq_type), aa = 'chemistry', 
                          dna = 'nucleotide', rna = 'nucleotide', 
                          other='nucleotide')
      
    }
    
    
    # Pick from default color schemes
    cs = switch(col_scheme, 
                # Color scheme based on chemistry of amino acids
                chemistry2 = data.frame(
                  letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
                  group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
                  col = c(rep('#058644', 5), rep('#720091', 2), rep('#0046C5', 3), rep('#C5003E', 2), rep('#2E2E2E', 8)),
                  stringsAsFactors = F
                ), 
                
                # Color scheme based on chemistry of amino acids
                chemistry = data.frame(
                  letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
                  group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
                  col = c(rep('#109648', 5), rep('#5E239D', 2), rep('#255C99', 3), rep('#D62839', 2), rep('#221E22', 8)),
                  stringsAsFactors = F
                ), 
                
                # Hydrophobicity index (PMID: 7108955) from -4.5 to 4.5
                hydrophobicity = data.frame(
                  letter = c('I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'W', 
                             'S', 'Y', 'P', 'H', 'D', 'E', 'N', 'Q', 'K', 'R'),
                  group = c(4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.9, -0.8,
                            -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5),
                  stringsAsFactors=F
                ), 
                
                # Colour based on nucleotide
                nucleotide2 = data.frame(
                  letter = c('A', 'C', 'G', 'T', 'U'),
                  col = c('darkgreen', 'blue', 'orange', 'red', 'red'),
                  stringsAsFactors = F
                ), 
                
                #alt red BA1200
                nucleotide = data.frame(
                  letter = c('A', 'C', 'G', 'T', 'U'),
                  col = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                  stringsAsFactors = F
                ), 
                
                base_pairing = data.frame(
                  letter = c('A', 'T', 'U', 'G', 'C'),
                  group = c(rep('Weak bonds', 3), rep('Strong bonds', 2)),
                  col = c(rep('darkorange', 3), rep('blue', 2)),
                  stringsAsFactors = F
                ),
                
                # ClustalX color scheme: 
                # http://www.jalview.org/help/html/colourSchemes/clustal.html
                clustalx = data.frame(
                  letter = c('W', 'L', 'V', 'I', 'M', 'F', 'A', 'R', 'K', 'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'Y', 'C', 'G', 'P'),
                  col = c(rep('#197FE5', 7), rep('#E53319', 2), rep('#19CC19', 4), rep('#CC4CCC', 2), 
                          rep('#19B2B2', 2), '#E57F7F', '#E5994C', '#B0B000'),
                  stringsAsFactors = F
                ),
                
                # Taylor color scheme (PMID: 9342138)
                taylor = data.frame(
                  letter = c('D','S','T','G','P','C','A','V','I','L','M','F','Y','W','H','R','K','N','Q','E'),
                  col = c('#FF0000','#FF3300','#FF6600','#FF9900','#FFCC00','#FFFF00','#CCFF00','#99FF00',
                          '#66FF00','#33FF00','#00FF00','#00FF66','#00FFCC','#00CCFF','#0066FF','#0000FF',
                          '#6600FF','#CC00FF','#FF00CC','#FF0066'),
                  stringsAsFactors = F
                )
    )
    
    if(!'group' %in% names(cs)) cs$group = cs$letter
    
    # Set attributes
    attr(cs, 'cs_label') = col_scheme
    class(cs) = c('data.frame','ggseqlogo_cs')
    
    return(cs)
  }
  
  
  
  
  
  #' Create new sequence logo color scheme
  #' 
  #' @param chars Vector of one letter characters 
  #' @param groups Vector of groups for letters with same length as chars (optional if cols parameter is provided) 
  #' @param cols Vector of colors with same length as chars (optional if values parameter is provided) 
  #' @param values Vector of numerical values with same length as chars
  #' @param name Name of color scheme
  #' 
  #' @export
  #' 
  #' @importFrom grDevices col2rgb
  #' @examples 
  #' 
  #' # Discrete color scheme examples
  #' cs1 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), groups=c('g1', 'g1', 'g2', 'g2'), 
  #'                       cols=c('red', 'red', 'blue', 'blue'), name='custom1')
  #' 
  #' cs2 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), cols=c('red', 'red', 'blue', 'blue'), 
  #'                       name='custom2')
  #' 
  #' # Quantitative color scheme
  #' cs3 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), values=1:4, name='custom3')
  make_col_scheme <- function(chars=NULL, groups=NULL, cols=NULL, values=NULL, name=''){
    
    
    if(is.null(chars) | any(nchar(chars) != 1) | !is.character(chars))
      stop('"chars" must be a character vector of one letter characters')
    
    
    if(is.null(values)){
      # Discrete colour scheme
      
      # Error check lengths
      if(length(chars) != length(cols)) stop('"chars" and "cols" must have same length')
      # Error check types
      if(!is.character(cols)) stop('"cols" must be a character vector')
      
      # Check valid colours
      tmp = col2rgb(cols); rm(tmp)
      
      if(is.null(groups)) groups = chars
      
      cs = data.frame( letter=chars, group=groups, col=cols, stringsAsFactors = F )
      
    }else{
      
      # Quantitative color scheme
      if(length(chars) != length(values)) stop('"chars" and "values" must have same length')
      cs = data.frame( letter=chars, group=values, stringsAsFactors=F )
    }
    
    # Remove duplicate letters
    cs = cs[!duplicated(cs$letter),]
    
    # Set attributes
    attr(cs, 'cs_label') = name
    class(cs) = c('data.frame','ggseqlogo_cs')
    
    return(cs)
  }
  
  #### ggseqlogo.R ####
  # if(T){
  #   require(ggplot2)
  #   setwd('~/Development/ggseqlogo/')
  #   source('R/heights.r')
  #   source('R/col_schemes.r')
  #   GGSEQLOGO_FONT_BASE = '~/Development/ggseqlogo/inst/fonts/'
  # }
  
  
  # Change range of values
  newRange <- function(old_vals, new_min=0, new_max=1){
    old_min = min(old_vals)
    old_max = max(old_vals)
    
    new_vals = (((old_vals - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
    new_vals
  }
  
  
  #' List fonts available in ggseqlogo
  #' 
  #' @param v If true, font names are printed to stderr. Otherwise, font names are returned as a character vector
  #' @export
  list_fonts <- function(v=T){
    
    fonts = c('helvetica_regular','helvetica_bold', 'helvetica_light',
              'roboto_medium','roboto_bold', 'roboto_regular',
              'akrobat_bold', 'akrobat_regular', 
              'roboto_slab_bold', 'roboto_slab_regular', 'roboto_slab_light', 
              'xkcd_regular')
    if(!v) return(fonts)
    message('Available ggseqlogo fonts:')
    for(f in fonts) message('\t', f)
  }
  
  
  # Read font from file if not in global envir.
  get_font <- function(font){
    
    GGSEQLOGO_FONT_BASE = getOption('GGSEQLOGO_FONT_BASE')
    if(is.null(GGSEQLOGO_FONT_BASE)){
      GGSEQLOGO_FONT_BASE=system.file("extdata", "", package = "ggseqlogo")
      options(GGSEQLOGO_FONT_BASE=GGSEQLOGO_FONT_BASE)
    }
    
    #all_fonts = c('sf_bold', 'sf_regular', 'ms_bold', 'ms_regular', 'xkcd_regular')
    font = match.arg(tolower(font), list_fonts(F))
    font_filename = paste0(font, '.font')
    font_obj_name = sprintf('.ggseqlogo_font_%s', font)
    
    font_obj = getOption(font_obj_name)
    if(is.null(font_obj)){
      # Not loaded into global env yet - load it into options
      font_path = file.path(GGSEQLOGO_FONT_BASE, font_filename)
      font_obj_list = list( tmp=readRDS(font_path) )
      names(font_obj_list) = font_obj_name
      options(font_obj_list)
      font_obj = font_obj_list[[1]]
    }
    
    # Return font data
    font_obj
  }
  
  
  # Generate height data for logo
  logo_data <- function( seqs, method='bits', stack_width=0.95, 
                         rev_stack_order=F, font, seq_group=1, 
                         seq_type = 'auto', namespace=NULL ){
    
    # Get font 
    font_df = get_font(font)
    
    # TODO
    # hh = twosamplelogo_method(seqs, seqs_bg, pval_thresh=0.05, seq_type = seq_type, namespace = namespace)
    
    # Generate heights based on method
    if(method == 'bits'){
      hh = bits_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
    }else if(method == 'probability'){
      hh = probability_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
    }else if(method == 'custom'){
      if(seq_type == 'auto') seq_type = guessSeqType(rownames(seqs))
      hh = matrix_to_heights(seqs, seq_type, decreasing = rev_stack_order)
    }else{
      stop('Invalid method!')
    }
    
    # Merge font df and heights
    ff = merge(font_df, hh, by = 'letter')
    # Scale x and ys to new positions
    x_pad = stack_width/2
    ff$x = newRange(ff$x, ff$position - x_pad, ff$position + x_pad)
    ff$y = newRange(ff$y, ff$y0, ff$y1)
    
    # Rename columns
    ff = as.data.frame(ff)[,c('x', 'y', 'letter', 'position', 'order')]
    ff$seq_group = seq_group
    
    # Set sequence type as attribute, to be used downstream
    attr(ff, 'seq_type') = attr(hh, 'seq_type')
    
    # Return data table
    ff
  }
  
  #' ggseqlogo custom theme
  #' 
  #' @param base_size font base size
  #' @param base_family font base family
  #' 
  #' @import ggplot2
  #' @export
  theme_logo <- function(base_size=12, base_family=''){
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace% 
      theme(panel.grid = element_blank(), legend.position = 'bottom', 
            axis.text.x=element_text(colour="black"),
            axis.text.y=element_text(colour="black"))
  }
  
  #' Plots sequence logo as a layer on ggplot 
  #' 
  #' @param data Character vector of sequences or named list of sequences. All sequences must have same width.
  #' @param method Height method, can be one of "bits" or "probability" (default: "bits")
  #' @param seq_type Sequence type, can be one of "auto", "aa", "dna", "rna" or "other" 
  #' (default: "auto", sequence type is automatically guessed)
  #' @param namespace Character vector of single letters to be used for custom namespaces. Can be alphanumeric, including Greek characters.
  #' @param font Name of font. See \code{list_fonts} for available fonts.
  #' @param stack_width Width of letter stack between 0 and 1 (default: 0.95)
  #' @param rev_stack_order If \code{TRUE}, order of letter stack is reversed (default: FALSE)
  #' @param col_scheme Color scheme applied to the sequence logo. See \code{list_col_schemes} for available fonts.
  #' (default: "auto", color scheme is automatically picked based on \code{seq_type}). 
  #' One can also pass custom color scheme objects created with the \code{make_col_scheme} function
  #' @param low_col,high_col Colors for low and high ends of the gradient if a quantitative color scheme is used (default: "black" and "yellow").
  #' @param na_col Color for letters missing in color scheme (default: "grey20")
  #' @param plot If \code{FALSE}, plotting data is returned 
  #' @param ... Additional arguments passed to layer params
  #' 
  #' @export
  #' @import ggplot2
  #' 
  #' @examples
  #' # Load sample data
  #' data(ggseqlogo_sample)
  #' 
  #' # Produce single sequence logo using geom_logo
  #' p1 = ggseqlogo( seqs_dna[[1]] ) 
  #' 
  geom_logo <- function (data = NULL, method = "bits", seq_type = "auto", namespace = NULL, 
                         font = "roboto_medium", stack_width = 0.95, rev_stack_order = F, 
                         col_scheme = "auto", low_col = "black", high_col = "yellow", 
                         na_col = "grey20", plot = T, position_colors = NULL, ...) 
  {
    if (stack_width > 1 | stack_width <= 0) 
      stop("\"stack_width\" must be between 0 and 1")
    if (is.null(data)) 
      stop("Missing \"data\" parameter!")
    if (!is.null(namespace)) 
      seq_type = "other"
    if (is.null(position_colors) & col_scheme != "position") 
      stop("Missing \"Position Colors\" parameter!")
    all_methods = c("bits", "probability", "custom")
    pind = pmatch(method, all_methods)
    method = all_methods[pind]
    if (is.na(method)) 
      stop("method must be one of 'bits' or 'probability', or 'custom'")
    if (is.character(data) | is.matrix(data)) 
      data = list(`1` = data)
    if (is.list(data)) {
      if (is.null(names(data))) 
        names(data) = seq_along(data)
      lvls = names(data)
      data_sp = lapply(names(data), function(n) {
        curr_seqs = data[[n]]
        logo_data(seqs = curr_seqs, method = method, stack_width = stack_width, 
                  rev_stack_order = rev_stack_order, seq_group = n, 
                  seq_type = seq_type, font = font, namespace = namespace)
      })
      data = do.call(rbind, data_sp)
      data$seq_group = factor(data$seq_group, levels = lvls)
    }
    if (!plot) 
      return(data)
    seq_type = attr(data, "seq_type")
    cs = get_col_scheme(col_scheme, seq_type)
    legend_title = attr(cs, "cs_label")
    if (col_scheme == "position") {
      data = merge(data, position_colors, by = "position", 
                   all.x = T)
    }
    else {
      data = merge(data, cs, by = "letter", all.x = T)
    }
    data = data[order(data$order), ]
    colscale_gradient = is.numeric(cs$group)
    colscale_opts = NULL
    if (col_scheme == "position") {
      colscale_opts = scale_fill_manual(breaks = position_colors$group, 
                                        values = position_colors$col, name = "Fundamental Residues", 
                                        na.value = na_col)
    }
    else {
      if (colscale_gradient) {
        colscale_opts = scale_fill_gradient(name = legend_title, 
                                            low = low_col, high = high_col, na.value = na_col)
      }
      else {
        tmp = cs[!duplicated(cs$group) & !is.na(cs$group), 
        ]
        col_map = unlist(split(tmp$col, tmp$group))
        colscale_opts = scale_fill_manual(values = col_map, 
                                          name = legend_title, na.value = na_col)
      }
    }
    guides_opts = NULL
    if (identical(cs$letter, cs$group)) 
      guides_opts = guides(fill = F)
    y_lim = NULL
    extra_opts = NULL
    if (method == "tsl") {
      y_lab = "Depleted    Enriched"
      tmp = max(abs(data$y))
      row_a = row_b = data[1, ]
      row_a$y = -tmp
      row_b$y = tmp
      data = rbind(data, row_a, row_b)
      data$facet = factor(data$y > 0, c(T, F), c("Enriched", 
                                                 "Depleted"))
      extra_opts = NULL
    }
    else if (method == "custom") {
      y_lab = ""
    }
    else {
      y_lab = method
      substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
    }
    data$group_by = with(data, interaction(seq_group, letter, 
                                           position))
    data$x = data$x
    logo_layer = layer(stat = "identity", data = data, mapping = aes_string(x = "x", 
                                                                            y = "y", fill = "group", group = "group_by"), geom = "polygon", 
                       position = "identity", show.legend = NA, inherit.aes = F, 
                       params = list(na.rm = T, ...))
    breaks_fun = function(lim) {
      1:floor(lim[2]/1.05)
    }
    list(logo_layer, scale_x_continuous(breaks = breaks_fun, 
                                        labels = identity), ylab(y_lab), xlab(""), colscale_opts, 
         guides_opts, coord_cartesian(ylim = y_lim), extra_opts)
  }
  
  
  
  #' Quick sequence logo plot
  #' 
  #' @description \code{ggseqlogo} is a shortcut for generating sequence logos. 
  #' It adds the ggseqlogo theme \code{\link{theme_logo}} by default, and facets when multiple input data are provided. 
  #' It serves as a convenient wrapper, so to customise logos beyond the defaults here, please use \code{\link{geom_logo}}.
  #' 
  #' @param data Character vector of sequences or named list of sequences. All sequences must have same width
  #' @param facet Facet type, can be 'wrap' or 'grid'
  #' @param scales Facet scales, see \code{\link{facet_wrap}}
  #' @param ncol Number of columns, works only when \code{facet='wrap'}, see \code{\link{facet_wrap}}
  #' @param nrow Number of rows, same as \code{ncol}
  #' @param ... Additional arguments passed to \code{\link{geom_logo}}
  #' 
  #' @export
  #' @examples
  #' # Load sample data
  #' data(ggseqlogo_sample)
  #' 
  #' # Plot a single DNA sequence logo
  #' p1 = ggseqlogo( seqs_dna[[1]] )
  #' print(p1)
  #' 
  #' # Plot multiple sequence logos at once
  #' p2 = ggseqlogo( seqs_dna )
  #' print(p2)
  ggseqlogo <- function(data, facet='wrap', scales='free_x', ncol=NULL, nrow=NULL, ...){
    
    # Generate the plot with default theme
    p = ggplot() + geom_logo(data = data, ...) + theme_logo() 
    
    # If it's an inidivdual sequence logo, return plot
    if(!'list' %in% class(data)) return(p)
    
    # If we have more than one plot, facet
    facet_opts = c('grid', 'wrap')
    pind = pmatch(facet, facet_opts)
    facet = facet_opts[pind]
    if(is.na(facet)) stop("facet option must be set to 'wrap' or 'grid'")
    
    if(facet == 'grid'){
      p = p + facet_grid(~seq_group, scales = scales)
    }else if(facet == 'wrap'){
      p = p + facet_wrap(~seq_group, scales = scales, nrow = nrow, ncol = ncol)
    }
    
    # Return plot
    return(p)
  }
  
  
  #' List of aligned transcription factor binding sequences 
  #'
  #' @name seqs_dna
  #' @docType data
  #' @keywords data
  NULL
  
  #' List of aligned kinase-substrate binding sequences 
  #'
  #' @name seqs_aa
  #' @docType data
  #' @keywords data
  NULL
  
  #' List of position frequency matrices for transcription factors
  #'
  #' @name pfms_dna
  #' @docType data
  #' @keywords data
  NULL
  
  
  # message('-- running example')
  # load('data/ggseqlogo_sample.rda')
  # p = ggseqlogo(sample_data$seqs_dna, nrow=3)
  # d = p$layers[[1]]$data
  # print(p)
  
  #### heights.R ####
  # Namespaces
  .AA_NAMESPACE = function() c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  .DNA_NAMESPACE = function() c('A', 'T', 'G', 'C')
  .RNA_NAMESPACE = function() c('A', 'U', 'G', 'C')
  
  # Generate letter matrix from vector of sequences
  # 
  # @param input vector of sequences
  letterMatrix <- function(input){
    # Ensure kmers are the same length characters 
    seq.len = sapply(input, nchar)
    num_pos = seq.len[1]
    if(! all(seq.len == num_pos)) stop('Sequences in alignment must have identical lengths')
    
    # Construct matrix of letters
    split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )
    
    t( matrix(split, seq.len, length(split)/num_pos) )
  }
  
  # Guess sequence type based on letter matrix
  # 
  # @param sp letters
  guessSeqType <- function(sp){
    # Ensure we have something
    if(length( intersect(sp, c(.AA_NAMESPACE(), .DNA_NAMESPACE(),.RNA_NAMESPACE())) ) == 0)
      stop('Could not get guess seq_type. Please explicitly define sequence type or use "other" with custom namespaces.')
    
    dat = setdiff(intersect(sp, .AA_NAMESPACE()), c(.DNA_NAMESPACE(),.RNA_NAMESPACE()))
    if(length(dat) > 0){
      return('AA')
    }else if('U' %in% sp){
      return('RNA')
    }
    return('DNA')
  }
  
  
  # Find namespace
  # 
  # @param letter_mat Matrix of latters
  # @param seq_type Sequence type
  # @param namespace Alphabet
  findNamespace <- function(letter_mat, seq_type, namespace){
    
    # Get all letters in our alignment
    sp = as.character(letter_mat)
    
    # Other namespace
    if(seq_type == "other"){
      if(is.null(namespace)) 
        stop('seq_type of "other" must have a defined namespace')
      
      namespace = as.character(namespace)
      # Get unique
      namespace = unique( unlist(strsplit(namespace, '')) )
      
      
      # Validate
      non_alphanumeric = grepl('[^a-zA-Z0-9\u03bb\u03b1\u03b2\u0393\u03b3\u0394\u03b4\u03b5\u03b6\u03b7\u03b8\u0398\u03b9\u03ba\u039b\u039b\u03bc\u039e\u03be\u03a0\u03c0\u03c1\u03c3\u03c4\u03c5\u03a6\u03c6\u03c7\u03c8\u03a8\u03a9\u03c9]', namespace)
      if( any( non_alphanumeric ) )
        stop('All letters in the namespace must be alphanumeric')
      
      # Ensure there is something in each column
      # apply(letter_mat, 2, function(column_letters){
      #   int = intersect(namespace, column_letters)
      #   if(length(int) == 0)
      #     stop('The alignment has no letters in namespace match aligned sequences in at least one column')
      # })
      
    }else{
      if(!is.null(namespace)) 
        stop('For custom namespaces please set seq_type to "other"')
      
      # Guess sequence type
      if(seq_type == "auto")
        seq_type = guessSeqType(sp)
      
      # Get predefined namespace
      namespace = get( sprintf('.%s_NAMESPACE', toupper(seq_type)) )()
    }
    
    return(list(seq_type = toupper(seq_type), 
                namespace = namespace))
  }
  
  # Calcualte bits
  #
  # @param pwm Position weight matrix
  # @param N Number of letters in namespace
  # @param Nseqs Number of sequences in PWM
  computeBits <- function(pwm, N=4, Nseqs=NULL){
    Nseqs = attr(pwm, 'nongapped')
    H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
    e_n = 0
    if(!is.null(Nseqs)) e_n = (1/logb(2)) * (N-1)/(2*Nseqs) 
    
    R_i = log2(N) - (H_i  + e_n)
    # Set any negatives to 0
    R_i = pmax(R_i, 0)
    return(R_i)
  }
  
  # Construct relative frequency matrix
  # @param seqs aligned sequences as vector
  # @param seq_type sequence type
  # @param namespace letters used for matrix construction
  # @param keep_letter_mat Keep letter matrix for some height methods
  makePFM <- function(seqs, seq_type='auto', namespace=NULL, keep_letter_mat=F){
    
    letter_mat = NA
    if(is.matrix(seqs)){
      # Process matrix
      if(is.null(rownames(seqs))) stop('Matrix must have letters for row names')
      
      num_pos = ncol(seqs)
      
      # Get namespace
      ns = findNamespace(rownames(seqs), seq_type, namespace)
      namespace = ns$namespace
      seq_type = ns$seq_type
      
      nseqs = NULL
      
      bg_prob = NA
      pfm_mat = seqs
      pfm_mat = apply(pfm_mat, 2, function(x) x / sum(x, na.rm=T))
      
      missing_rows = setdiff(namespace, rownames(pfm_mat))
      
      if(length(missing_rows) > 0){
        miss = matrix(rep(0, length(missing_rows) * ncol(pfm_mat)), nrow=length(missing_rows), dimnames = list(missing_rows))
        pfm_mat = rbind(pfm_mat, miss)
      }
      
      pfm_mat = pfm_mat[namespace,]
      
    }else{
      # Process sequences
      
      # Number of positions in alignment
      num_pos = nchar(seqs[1])
      # Number of sequences
      nseqs = length(seqs)
      # Letter matrix
      letter_mat = letterMatrix(seqs)
      
      
      # Get namespace
      ns = findNamespace(letter_mat, seq_type, namespace=namespace)
      namespace = ns$namespace
      seq_type = ns$seq_type
      
      # Construct PWM
      pfm_mat = apply(letter_mat, 2, function(pos.data){
        # Get frequencies 
        t = table(pos.data)
        # Match to aa
        ind = match(namespace, names(t))
        # Create column
        col = t[ind]
        col[is.na(col)] = 0
        names(col) = namespace
        # Do relative frequencies
        col = col / sum(col)
        col
      })
      
      mat = matrix((letter_mat %in% namespace), nrow=nrow(letter_mat))
      attr(pfm_mat, 'nongapped') = apply(mat, 2, sum)
      attr(pfm_mat, 'nseqs') = nseqs
    }
    
    # Number of letters in ns
    N = length(namespace)
    
    # Assign seq type and namespace as attributes
    attr(pfm_mat, 'seq_type') = seq_type
    attr(pfm_mat, 'namespace') = namespace
    
    # Non-gapped columns
    if(seq_type == 'aa') namespace = c(namespace, 'X', 'B', 'Z')
    
    # Information content
    attr(pfm_mat, 'bits') = computeBits(pfm_mat, N, nseqs)
    
    # Assign AA names to rows/pos col
    rownames(pfm_mat) = namespace
    colnames(pfm_mat) = 1:num_pos
    
    if(keep_letter_mat) return(list(letter_mat = letter_mat, pfm=pfm_mat))
    
    return(pfm_mat)
  }
  
  
  
  ######################
  # Matrix to heights
  ####
  
  # General function to convert matrix of heights to polygon data frame 
  # @param mat matrix of heghts
  # @param seq_type sequence type
  # @decreasing Sets order of letters, high to low or low to high
  matrix_to_heights <- function(mat, seq_type, decreasing=T){
    
    mat[is.infinite(mat)] = 0 
    
    if(any(duplicated(rownames(mat)))) stop('Matrix input must have unique row names')
    
    dat = lapply(1:ncol(mat), function(i){
      vals = mat[,i]
      
      pos = sort( vals[vals >= 0], decreasing = decreasing)
      neg = sort(vals[vals < 0], decreasing = !decreasing)
      #vals = sort(vals, decreasing = T)
      cs_pos = cumsum( pos )
      cs_neg = cumsum( neg )
      
      df_pos = df_neg = NULL
      
      if(length(pos) > 0)
        df_pos = data.frame(letter=names(pos), position=i,  y0=c(0, cs_pos[-length(cs_pos)]), 
                            y1=cs_pos, stringsAsFactors = F)
      
      if(length(neg) > 0)
        df_neg = data.frame(letter=names(neg), position=i, y0=cs_neg, y1=c(0, cs_neg[-length(cs_neg)]), 
                            stringsAsFactors = F)
      
      rbind(df_pos, df_neg)
    })
    
    dat = do.call(rbind, dat)
    
    # Adjust y spacing 
    space_factor = 0.004
    y_pad = max(dat$y1) * space_factor
    dat$y0 = dat$y0 + y_pad
    dat = subset(dat, dat$y1 > dat$y0)
    
    # Dummy points to make sure full plot is drawn
    # Make sure position 1 and n have a dummy empty letter missing
    dummy = data.frame(letter=dat$letter[1], position=NA, y0=0, y1=0)
    
    # Missing first position
    if(dat$position[1] != 1){
      dummy$position = 1
      dat = rbind( dummy, dat )
    }
    
    # Missing last position
    if(dat$position[nrow(dat)] != ncol(mat)){
      dummy$position = ncol(mat)
      dat = rbind( dat, dummy )
    }
    
    rownames(dat) = NULL
    
    attr(dat, 'seq_type') = seq_type
    
    dat
  }
  
  
  
  # Shannon entropy method
  bits_method <- function(seqs, decreasing, ...){
    # Make PFM
    pfm = makePFM(seqs, ...)
    
    # Get ic
    ic = attr(pfm, 'bits')
    if(all(ic == 0)){
      warning('All positions have zero information content perhaps due to too few input sequences. Setting all information content to 2.')
      ic = (ic * 0)+2
    }
    heights = t(t(pfm) * ic)
    
    seq_type = attr(pfm, 'seq_type')
    matrix_to_heights(heights, seq_type, decreasing)
  } 
  
  # Probability method
  probability_method <- function(seqs, decreasing, ...){
    # Make PFM
    pfm = makePFM(seqs, ...)
    seq_type = attr(pfm, 'seq_type')
    matrix_to_heights(pfm, seq_type, decreasing)
  }
  
  
  
  
  
  
  
  
  
  #####
  # MAIN
  # setwd(project_folder)
  design_alanine <- read.csv(design_alanine_file, sep = "\t")
  
  input_data_path <- sprintf("%s/inputs/04_individual_serums_raw_data/", project_folder)
  files <- list.files(input_data_path)
  files <- files[grepl("raw.tsv",files)]
  
  #Get only the ids for the AlaScan
  design_data <- fread(file = design_data_file, header = T, sep = "\t", na.strings = NULL)
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
  
  proteins_to_process <- if (!is.null(selected_protein)) {
    unique(design_alanine$protein[design_alanine$protein == selected_protein])
  } else {
    unique(design_alanine$protein)
  }
  
  for (protein_to_estimate in proteins_to_process) {
    dt <- estimateSignalChanges(protein_to_estimate)
    split_dt <- splitData(dt)
    # sequence_logos <- lapply(split_dt, plotLogo)
    
    matrix <- lapply(split_dt, generateHeatmapMatrix)
    #### SAVE RESULTS ####
    write.table(dt,paste0(project_folder, "/outputs/07_alanine_scan_raw_data/Raw_data_signal_change_alanine_scan_", protein_to_estimate,".tsv"),sep = "\t")  #long format
    
    for (i in seq_along(matrix)) {
      write.table(
        matrix[[i]],
        paste0(project_folder, "/outputs/07_alanine_scan_raw_data/Raw_data_signal_change_matrix_alanine_scan_", protein_to_estimate, "_group_", i, ".tsv"), 
        sep = "\t", row.names = TRUE) #matrix for heatmap
    }
    
    # using "pheatmap" package you could visualize results as in the manuscript.
    if (heatmap_plot) {
      for (i in seq_along(matrix)) {
        pdf(
          file = paste0(project_folder, "/outputs/07_alanine_scan_raw_data/Heatmap_alanine_scan_", protein_to_estimate, "_group_", i, ".pdf"),
          height = 12, width = 12
        )
        print(make_heatmap(temp_matrix = matrix[[i]], selected_protein = paste(protein_to_estimate, "group", i)))
        dev.off()
      }
    }
    
    if (sequence_logo_plot) {
      for (i in seq_along(split_dt)){
      # for (i in seq_along(sequence_logos)){
        pdf(
          file = paste0(project_folder, "/outputs/07_alanine_scan_raw_data/sequence_logo_alanine_scan_", protein_to_estimate, "_group_", i, "_", sequence_logo_source ,".pdf"),
          height = 12, width = 12
        )
        print(plotLogo(sub_signal_change_per_position = split_dt[[i]],
        # print(plotLogo(sub_signal_change_per_position = sequence_logos[[i]], 
                            protein_id = protein_to_estimate, 
                            sources = sequence_logo_source,
                            plot_title = paste(protein_to_estimate, "group", i, sequence_logo_source, "source")))
        dev.off()
      }
    }
    
  }
  
}