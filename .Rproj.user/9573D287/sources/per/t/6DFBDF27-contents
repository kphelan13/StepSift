#' Install Required Packages
#' This function checks if the necessary packages are installed and installs them if not.
#' @export
install_required_packages <- function() {
  required_packages <- c("DESeq2", "tidyverse", "Biobase", "stringi", "plyr", "PCAtools", "scales", "ggplot2")
  
  # Check for missing packages
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  
  # Install missing packages
  if(length(missing_packages)) {
    install.packages(missing_packages)
  }
  
  # Load all required packages
  invisible(lapply(required_packages, library, character.only = TRUE))
}
