# Simulate Models Using Sampled Parameters

Runs the epidemiological model multiple times using different parameter
sets sampled from MCMC results. Each simulation is optionally prewarmed
using a specified number of years prior to the start date.

## Usage

``` r
simulate_models(
  model,
  param_inputs,
  sampled_params,
  start_date,
  end_date,
  prewarm_years = 2,
  days_per_year = 360
)
```

## Arguments

- model:

  The model function to simulate from.

- param_inputs:

  A list of baseline parameters to be updated with each sampled
  parameter set.

- sampled_params:

  A list of named parameter sets, typically obtained from posterior
  samples.

- start_date:

  Character string giving the start date for the simulation (e.g.,
  "2021-01-01").

- end_date:

  Character string giving the end date for the simulation (e.g.,
  "2021-12-31").

- prewarm_years:

  Integer specifying how many years before `start_date` to simulate for
  model prewarming (default is 2).

- days_per_year:

  Integer specifying the number of days per simulation year (default is
  360).

## Value

A list of data frames. Each data frame contains the simulation results
for a different parameter set.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example setup with a dummy model and parameters
model <- your_model_function  # Replace with your actual model function
param_inputs <- return_default_inputs()  # Replace with your actual input generator

# Simulate posterior parameter samples
posterior_samples <- matrix(rnorm(30), nrow = 10, ncol = 3)
colnames(posterior_samples) <- c("alpha", "beta", "gamma")

# Convert to list of named parameter sets
sampled_params <- lapply(1:nrow(posterior_samples), function(i) {
  as.list(posterior_samples[i, ])
})

# Run simulations
simulations <- simulate_models(
  model = model,
  param_inputs = param_inputs,
  sampled_params = sampled_params,
  start_date = "2021-01-01",
  end_date = "2021-12-31"
)
} # }
```
