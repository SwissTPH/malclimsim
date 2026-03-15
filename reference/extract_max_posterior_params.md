# Extract Parameters with Maximum Log Posterior

This function finds and extracts the parameter set corresponding to the
maximum log posterior value from MCMC output. It assumes that the MCMC
results object includes a `coda_pars` matrix or data frame with columns
for `log_prior`, `log_likelihood`, and `log_posterior`, along with
parameter names.

## Usage

``` r
extract_max_posterior_params(results)
```

## Arguments

- results:

  A list containing MCMC output, with a `coda_pars` element (a matrix or
  data frame) that includes a column named `"log_posterior"` and columns
  for sampled parameters.

## Value

A named numeric vector of parameter values corresponding to the maximum
log posterior.

## Examples

``` r
# Create dummy MCMC output for illustration
set.seed(123)
results <- list(
  coda_pars = data.frame(
    param1 = rnorm(100),
    param2 = rnorm(100),
    log_prior = rnorm(100),
    log_likelihood = rnorm(100),
    log_posterior = rnorm(100)
  )
)

# Extract parameters corresponding to the maximum log posterior
extract_max_posterior_params(results)
#>      param1    param2
#> 27 0.837787 0.2353866
```
