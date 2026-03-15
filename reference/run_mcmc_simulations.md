# Run Model Simulations for Each Parameter Set

This function runs the epidemiological model using different parameter
sets sampled from the MCMC posterior.

## Usage

``` r
run_mcmc_simulations(
  model,
  param_inputs,
  param_samples,
  start_date,
  end_date,
  prewarm_years = 2,
  days_per_year = 360
)
```

## Arguments

- model:

  A function representing the epidemiological model, typically compiled
  via `odin.dust`.

- param_inputs:

  A list of default model parameters to be updated with sampled values.

- param_samples:

  A matrix or data frame of parameter sets sampled from the MCMC
  posterior. Each row is a distinct sample.

- start_date:

  Character string or Date object indicating the simulation start date
  (e.g., "2014-01-01").

- end_date:

  Character string or Date object indicating the simulation end date
  (e.g., "2022-12-31").

- prewarm_years:

  Number of years for model prewarming to reach equilibrium before the
  start date (default: 2).

- days_per_year:

  Number of simulation time steps per year (default: 360).

## Value

A list of data frames, each containing the simulated compartment outputs
for one parameter set.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume model is loaded via `load_model()` and priors are defined
model <- load_model("model_det_1")

# Get default parameter inputs and samples from MCMC
param_inputs <- return_default_param_inputs()
param_samples <- matrix(
  c(0, 0, 2, 5, 3, 0, 1, 0.5, 0.8, 1),  # Example parameter set (one row)
  nrow = 1,
  byrow = TRUE
)
colnames(param_samples) <- c("lag_R", "lag_T", "alpha", "sigma_LT", "sigma_RT",
                              "R_opt", "k1", "eff_SMC", "s", "b")  # Adjust as needed

# Run simulation
simulations <- run_mcmc_simulations(
  model,
  param_inputs,
  param_samples,
  start_date = "2014-01-01",
  end_date = "2022-12-31",
  prewarm_years = 2,
  days_per_year = 360
)

# View first simulation output
head(simulations[[1]])
} # }
```
