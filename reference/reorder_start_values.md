# Reorder MCMC Start Values to Match Prior Specification

Ensures that the rows of a matrix of MCMC start values are ordered
consistently with the names in the `param_priors` list. This is useful
when passing initial values to an MCMC routine that expects a specific
parameter order.

## Usage

``` r
reorder_start_values(start_values, param_priors)
```

## Arguments

- start_values:

  A matrix of initial parameter values for MCMC, where each row
  corresponds to a parameter and each column to a chain or
  initialization. The row names should be the parameter names.

- param_priors:

  A named list of prior specifications (e.g., as used by
  `mcstate::pmcmc_parameters$new()`), where names correspond to
  parameter names.

## Value

A reordered matrix of start values with the row order matching the order
of names in `param_priors`. Parameters missing from `start_values` are
dropped (with a warning optionally).

## Details

The function matches parameter names from `start_values` to those in
`param_priors`, and reorders them accordingly. Any parameters present in
`param_priors` but missing in `start_values` are ignored in the returned
matrix.

## Examples

``` r
start_vals <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, dimnames = list(c("beta", "gamma"), NULL))
priors <- list(gamma = list(min = 0, max = 1, prior = dunif),
               beta = list(min = 0, max = 1, prior = dunif))
reordered <- reorder_start_values(start_vals, priors)
```
