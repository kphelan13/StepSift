#' Stepwise Variance Regression
#' Iterative regression to assess covariate variance colinearilty
#'
#' @param data Matrix of normalized RNA-seq data with genes in rows and samples in columns.
#' @param metadata Dataframe of metadata for RNA-seq samples. Samples should be in rows and covariates in columns.
#' @param covariates List of covariates in `metadata` to be assessed.
#' @param numeric_covariates Character vector of covariates that are numeric; others will be treated as factors.
#' @param output_dir Path to save output files.
#' @param rank Integer specifying the number of PCs to retain (from PCAtools). Default is 10. Use NULL for all PCs.
#' 
#' @return Dataframe of covariates ordered by greatest association with variance after each iterative regression.
#' @export
#' @examples
#' # Example usage:
#' SVR(data = my_counts_matrix, metadata = my_metadata_df,
#'     covariates = c("Age", "Sex"), numeric_covariates = "Age",
#'     output_dir = "path/to/output", rank = 10)

SVR <- function(data, metadata, covariates, numeric_covariates, output_dir, rank = 10) {
  
  #Format Metadata and specify numeric covariates
  metadata[,numeric_covariates] <- lapply(metadata[,numeric_covariates], as.numeric)
  metadata <- metadata %>% mutate_if(is.character, as.factor) %>% 
    mutate_if(is.integer, as.factor)
  
  # Perform PCA
  pca_data <- pca(data, metadata = metadata, scale = TRUE, rank = rank)
  
  # Initialize data structures
  df_iter <- list()
  data1 <- merge(pca_data$rotated, pca_data$metadata, by = 'row.names')
  data0 <- data1
  current_variables <- c()
  
  df <- NULL  # Initialize an empty dataframe to store results
  
  # Initialize a list to store dataframes for each iteration
  df_iter <- list()
  
  # Iterate 10 times
  for (iteration in 1:length(covariates)) {
    # Reset df to an empty dataframe at the beginning of each iteration
    df <- NULL
    
    # Iterate through a list of variables
    for (v in covariates) {
      
      if (!(v %in% current_variables)) {
        var_explained = 0.0  # Initialize the variance explained
        
        # Iterate through principal components
        for (p in names(pca_data$variance)) {
          f <- formula(paste(p, paste(v, sep = '+'), sep = '~'))  # Create a formula
          m <- lm(f, data = data1)  # Fit a linear regression model
          s <- summary(m)  # Get model summary
          var_explained = var_explained + s$r.squared * pca_data$variance[p]  # Calculate variance explained
        }
        
        ddf <- data.frame(variable = v, variance_explained = var_explained, row.names = v)  # Create a dataframe with results
        df <- rbind(df, ddf)  # Add results to the main dataframe
      } else {
        df <- rbind(df, data.frame(variable = v, variance_explained = NA, row.names = v))  # If the variable is already in use, add NA
      }
    }
    
    row.names(df) <- df$variable  # Set row names to the variable names
    df <- subset(df[order(df$variance_explained, decreasing = TRUE), ], select = -variable)  # Order dataframe by variance explained
    head(df, 1)  # Display the first row of the dataframe
    
    current_variables <- c(current_variables, row.names(df)[1])  # Update the list of current variables
    colnames(df) <- paste0(colnames(df), length(current_variables))  # Update column names
    df_iter[[length(current_variables)]] <- df  # Store the dataframe in a list
    
    # Iterate through principal components and update data1 with residuals
    for (p in names(pca_data$variance)) {
      f <- formula(paste(p, paste(current_variables, collapse = '+'), sep = '~'))  # Create a formula
      m <- lm(f, data = data0, na.action = na.exclude)  # Fit a linear regression model
      data1[, p] <- residuals(m)  # Update data1 with residuals
    }
    
    df_iter[[iteration]] <- df
  }
  
  m <- merge(df_iter[[1]], df_iter[[2]], by = 'row.names')
  row.names(m) <- m$Row.names
  m <- subset(m, select = -Row.names)
  
  for (i in 3:length(covariates)) {
    m <- merge(m, df_iter[[i]], by = 'row.names')
    row.names(m) <- m$Row.names
    m <- subset(m, select = -Row.names)
  }
  
  # Create a vector of new column labels
  new_col_names <- paste("Step", 1:ncol(m))
  
  # Assign the new column labels to 'm'
  colnames(m) <- new_col_names
  
  # Find the variable with the highest value in each column
  highest_vars <- apply(m, 2, function(col) rownames(m)[which.max(col)])
  
  # Reorder the rows in 'm' based on the identified variables
  m_reordered <- m[match(highest_vars, rownames(m)), ]
  
  
  # Save the results as a CSV file
  write.csv(m_reordered, file.path(output_dir, paste0('Stepwise_Variance_Regression', '.csv')))
  
  # Return the results or any other relevant information
  return(list(pca_data = pca_data, variance_explained_data = m))
}
