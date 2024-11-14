# StepSift

## Installation

```r
# Install the devtools package if needed
install.packages("devtools")

# Install your package from GitHub
devtools::install_github("kphelan13/StepSift")
```

## Example
```r
# Load the package
library(StepSift)

# Run the PVE function
PVE(counts = my_counts_matrix, metadata = my_metadata_df,
    covariates = c("Age", "Sex"), numeric_covariates = "Age",
    cutoff = 5, output_dir = "path/to/output", rank = 10)

