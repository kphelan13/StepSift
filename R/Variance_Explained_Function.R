#' Variance Explained Function
#' Calculate the percent variation associated with any covariate of interest
#'
#' @param data Matrix of normalized RNA-seq data with genes in rows and samples in columns.
#' @param metadata Dataframe of metadata for RNA-seq samples. Samples should be in rows and covariates in columns.
#' @param covariates List of covariates in `metadata` to be assessed.
#' @param numeric_covariates Character vector of covariates that are numeric; others will be treated as factors.
#' @param cutoff Optional threshold for a horizontal line in the variance barplot. Default is NULL.
#' @param output_dir Path to save output files.
#' @param rank Integer specifying the number of PCs to retain (from PCAtools). Default is 10. Use NULL for all PCs.
#' 
#' @return Dataframe with `variable` and `variance_explained` columns; generates a screeplot and a variance barplot.
#' @export
#' @examples
#' # Example usage:
#' PVE(data = my_counts_matrix, metadata = my_metadata_df,
#'     covariates = c("Age", "Sex"), numeric_covariates = "Age",
#'     cutoff = 5, output_dir = "path/to/output", rank = 10)

PVE <- function(data, metadata, covariates, numeric_covariates, cutoff = NULL, output_dir, rank = 10) {
  library(DESeq2)
  library(tidyverse)
  library(Biobase)
  library(stringi)
  library(plyr)
  library(PCAtools)
  library(scales)
  library(ggplot2)
  
  # Convert specified covariates to numeric
  metadata[, numeric_covariates] <- lapply(metadata[, numeric_covariates], as.numeric)
  
  # Ensure metadata columns are treated as factors
  metadata <- metadata %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.integer, as.factor)
  
  # Perform PCA
  pca_data <- pca(data, metadata = metadata, removeVar = 0.1, scale = TRUE, rank = rank)
  
  # Create scree plot
  screeplot(pca_data, axisLabSize = 9)
  ggsave(path = output_dir, filename =  "scree_plot.png", width = 7, height = 5)
  
  # Initialize data frames and variables
  df_iter <- list()
  data1 <- merge(pca_data$rotated, pca_data$metadata, by = 'row.names')
  data0 <- data1
  current_variables <- character(0)
  
  # Initialize a data frame to store results
  df <- data.frame(variable = character(0), variance_explained = numeric(0))
  
  for (v in covariates) {
    if (!(v %in% df$variable)) {
      var_explained <- 0.0
      
      for (p in names(pca_data$variance)) {
        f <- formula(paste(p, paste(v, sep = '+'), sep = '~'))
        m <- lm(f, data = data1)
        s <- summary(m)
        var_explained <- var_explained + (s$r.squared * pca_data$variance[p])
      }
      
      # Append results to the data frame
      df <- rbind(df, data.frame(variable = v, variance_explained = var_explained))
      df <- df[order(df$variance_explained, decreasing = TRUE), ]
      rownames(df) <- NULL
    }
  }
  
  # Save the result as a separate CSV file
  write.csv(df, file.path(output_dir, paste0('Variance_Explained', '.csv')), row.names = FALSE)
  
  # Create and save a bar graph
  gg <- ggplot(df, aes(x = variable, y = variance_explained)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
    xlab("Covariate") +
    ylab("% Variance Explained") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Variance Explained by Covariate")
  
  ggsave(path = output_dir, filename = "variance_explained_bar_graph.png", plot = gg, width = 7, height = 5)
  
  return(df)
}
