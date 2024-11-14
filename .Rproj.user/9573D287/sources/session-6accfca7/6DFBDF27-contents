#' Install Required Packages 
#' This function checks if the necessary packages are installed and installs them if not.
#' It uses BiocManager to install both Bioconductor and CRAN packages.
#' @export
install_required_packages <- function() {
  required_packages <- c("DESeq2", "tidyverse", "Biobase", "stringi", "plyr", "PCAtools", "scales", "ggplot2")
  
  # Check for missing packages
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  
  # Install missing packages
  if(length(missing_packages)) {
    # Install Bioconductor and CRAN packages using BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing_packages)
  }
  
  # Load all required packages
  invisible(lapply(required_packages, library, character.only = TRUE))
}
