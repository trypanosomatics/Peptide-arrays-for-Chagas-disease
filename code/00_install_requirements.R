#!/usr/bin/env Rscript
# Shebang, interpreter for executable

# Install BiocManager if it is not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  writeLines("Installing 'BiocManager'...")
  install.packages("BiocManager", repos = "http://cran.rstudio.com/")
}

# Install preprocessCore from Bioconductor
writeLines("Installing 'preprocessCore' from Bioconductor...")
BiocManager::install("preprocessCore", ask = FALSE)

# List of CRAN packages to install
packages <- c("data.table", "zoo", "dplyr", "reshape2", "pheatmap", "ggplot2", 
              "ggseqlogo", "gridExtra", "patchwork", "ggridges", "scico", 
              "ggh4x", "readr", "tidyr")

# Function to install packages if they are not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    writeLines(paste("Installing package:", pkg, "..."))
    install.packages(pkg)
  } else {
    writeLines(paste("Package already installed:", pkg))
  }
}

# Install all packages
sapply(packages, install_if_missing)
