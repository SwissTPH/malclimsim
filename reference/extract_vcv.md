# Extract Variance-Covariance Matrix and Restart Values

This function extracts the variance-covariance matrix and restart values
from the results of an MCMC run. It can optionally save the computed
proposal matrix and start values to specified file paths.

## Usage

``` r
extract_vcv(
  results,
  S_prev,
  save = TRUE,
  param_names,
  file_proposal = "",
  file_start = ""
)
```

## Arguments

- results:

  A list containing MCMC output with elements `mcmc_run$pars` (a matrix
  of parameter values) and `mcmc_run$chain` (a vector indicating chain
  membership).

- S_prev:

  Number of previous samples to include from each chain.

- save:

  Logical. If `TRUE`, saves the proposal matrix and start values to
  files.

- param_names:

  A character vector of parameter names to label the proposal matrix.

- file_proposal:

  File path to save the proposal covariance matrix (only if
  `save = TRUE`).

- file_start:

  File path to save the start values (only if `save = TRUE`).

## Value

A list with two elements:

- start_values:

  A matrix of start values (medians) for each parameter.

- proposal_matrix:

  A variance-covariance matrix used as a proposal distribution.

## Examples

``` r
# Simulate dummy MCMC output
set.seed(42)
results <- list(
  mcmc_run = list(
    pars = matrix(rnorm(1000), ncol = 5),
    chain = rep(1:5, each = 40)
  )
)
param_names <- paste0("param_", 1:5)
S_prev <- 20

# Run without saving files
output <- extract_vcv(
  results = results,
  S_prev = S_prev,
  save = FALSE,
  param_names = param_names,
  file_proposal = tempfile(fileext = ".rds"),
  file_start = tempfile(fileext = ".rds")
)
str(output)
#> List of 2
#>  $ : num [1:5, 1] -0.0881 0.0424 -0.0109 -0.1169 0.0987
#>  $ : num [1:5, 1:5] 0.1 0 0 0 0 0 0.1 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:5] "param_1" "param_2" "param_3" "param_4" ...
#>   .. ..$ : chr [1:5] "param_1" "param_2" "param_3" "param_4" ...
```
