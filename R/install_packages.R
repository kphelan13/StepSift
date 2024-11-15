#' Install Required Packages 
#' This function checks if the necessary packages are installed and installs them if not.
#' It uses BiocManager to install both Bioconductor and CRAN packages.
#' @export
install_required_packages <- function() {
  # Define the required packages
  required_packages <- c("DESeq2", "tidyverse", "Biobase", "stringi", "dplyr", "PCAtools", "scales", "ggplot2")
  
  # Install BiocManager if not already installed
  if (!"BiocManager" %in% rownames(installed.packages())) {
    install.packages("BiocManager")
  }
  
  # Identify packages not yet installed
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
  
  # Use BiocManager to install both CRAN and Bioconductor packages
  if (length(missing_packages)) {
    BiocManager::install(missing_packages, update = FALSE)
  }
  
  # Load all required packages
  invisible(lapply(required_packages, library, character.only = TRUE))
}