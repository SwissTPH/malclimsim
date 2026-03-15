# Extract and Print Maximum Log-Likelihood and Log-Posterior

This function extracts and optionally prints the maximum values of the
log-likelihood and log-posterior from the results of an MCMC run. It is
useful for quickly evaluating the fit of a model.

## Usage

``` r
max_ll_post(results, print_ll = TRUE)
```

## Arguments

- results:

  A list containing the MCMC results from an inference run, including
  the `coda_pars` matrix with columns for `log_likelihood` and
  `log_posterior`.

- print_ll:

  Boolean. If `TRUE`, the maximum log-likelihood and log-posterior
  values will be printed to the console (default: `TRUE`).

## Value

A list with two named elements:

- max_log_likelihood:

  The maximum value of the log-likelihood.

- max_log_posterior:

  The maximum value of the log-posterior.
