#### INSTALLING LIBRARIES ####

writeLines("Installing library 'data.table' for R")
install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'preprocessCore' from Bioconductor")
install.packages("BiocManager")
BiocManager::install("preprocessCore")

writeLines("Installing library 'zoo' for R")
install.packages("zoo", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'dplyr' for R")
install.packages("dplyr", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'reshape2' for R")
install.packages("reshape2", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'pheatmap' for R")
install.packages("pheatmap", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'ggplot2' for R")
install.packages("ggplot2", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'ggseqlogo' for R")
install.packages("ggseqlogo", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'grid' for R")
install.packages("grid", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'gridExtra' for R")
install.packages("gridExtra", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'patchwork' for R")
install.packages("patchwork", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'ggridges' for R")
install.packages("ggridges", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'scico' for R")
install.packages("scico", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'ggh4x' for R")
install.packages("ggh4x", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'readr' for R")
install.packages("readr", repos = "http://cran.rstudio.com/", dependencies = TRUE)

writeLines("Installing library 'tidyr' for R")
install.packages("tidyr", repos = "http://cran.rstudio.com/", dependencies = TRUE)