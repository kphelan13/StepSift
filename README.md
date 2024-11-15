# StepSift
## About
Adjusting for unwanted technical and biologic effects in large-scale genomics sequcing data is common practice. Quantification of the variance attributable to any given covariate and identification of covariates which are explaining similar variance provides a data-driven approach for covariate adjustment in a dataset-dependent manner. There are two major functions in this package:

1. PVE(Percent variance explained): This function uses linear dimensionality reduction to quantify the amount of variance due to any given covariate.

   
2. SVR (Stepwise variance regression): This function iteravely performs PVE, but at each step regresses the covariate explaining the most variance from the sequencing dataset and moves forward with the residual matrix. In effect, this function provides an ordered list of covariates accounting for those which may be explaining the same or similar variance.

## Methodology: Calculating Variance Explained by Covariates

This package calculates the variance explained by covariates in sequencing data using the following steps:

1. **Principal Components Analysis (PCA)**: PCA is performed on the sequencing data matrix `X`, where samples are rows and features are columns. Each principal component `PC_k` is associated with a variance explained, denoted as `VarExp_k`.

2. **Regression on Covariates**: For each sample `i`, let `s_ik` be the score of sample `i` on principal component `k`, and let `C_i` represent the covariate value for that sample. A regression is conducted between `s_ik` and `C_i` to obtain a regression coefficient `β_kC` for each covariate `C` on each principal component `k`.

3. **Weighting by Variance Explained**: Each regression coefficient `β_kC` is weighted by the variance explained by the corresponding principal component, giving a weighted coefficient: `β_kC * VarExp_k`.

4. **Summing Weighted Coefficients**: Finally, the weighted coefficients are summed across all principal components to calculate the total variance explained by each covariate `C`.

The formula for the final variance explained by each covariate `C` is:

Variance Explained_C = Σ ( β_kC * VarExp_k ) for k = 1 to K

where:

- `β_kC` is the regression coefficient between the score of principal component `k` and the covariate `C`,
- `VarExp_k` is the variance explained by principal component `k`, and
- `K` is the total number of principal components considered.

This final value represents the percentage of variance in the sequencing data that can be attributed to each covariate.

## Notes for use
The counts matrix should be any normalized dataset (ie Log2 normalized counts per million) with observations (ie genes) in rows and the samples in columns. The metadata should have matched sample IDs as rows with covarite IDs in the columns.

These functions have some dependencies, please make sure the following are installed: DESeq2, tidyverse, Biobase, stringi, dplyr, PCAtools, scales
ggplot2. The function 'install_required_packages' can do this for you.

Any comments, questions, or suggestions? Please email me at phelankj@mail.uc.edu


## Installation

```r
# Install the devtools package if needed
install.packages("devtools")

# Install package from GitHub
devtools::install_github("kphelan13/StepSift")

# Load Package and run 'install_packages' to install dependencies
library(StepSift)
install_packages()
```

## Example
```r
# Load the package
library(StepSift)

# Run the PVE function
PVE(counts = my_counts_matrix, metadata = my_metadata_df,
    covariates = c("Age", "Sex"), numeric_covariates = "Age",
    cutoff = 5, output_dir = "path/to/output", rank = 10)

