# Simulate Model Using Parameters with Maximum Log Posterior

Runs a model simulation using the parameter set corresponding to the
highest posterior density (i.e., the maximum log posterior value) from
the MCMC output. The function includes a prewarm period prior to the
specified `start_date` to stabilize the system before collecting
simulation results.

## Usage

``` r
simulate_with_max_posterior_params(
  results,
  start_date,
  end_date,
  model,
  prewarm_years = 2,
  days_per_year = 360,
  mu_transform_A = NULL,
  mu_transform_C = NULL,
  covariate_matrix = NULL
)
```

## Arguments

- results:

  A list containing MCMC results. Must include `log_posterior`,
  `param_inputs`, and samples of estimated parameters.

- start_date:

  Simulation start date (as `Date` or character string) for the output
  period.

- end_date:

  Simulation end date (as `Date` or character string).

- model:

  The model function to simulate. Should be compatible with
  [`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md).

- prewarm_years:

  Number of years to simulate prior to `start_date` to allow system
  stabilization (default is 2).

- days_per_year:

  Number of days in a simulation year (default is 360).

- mu_transform_A:

  Optional transformation function for adult incidence (default is
  `NULL`).

- mu_transform_C:

  Optional transformation function for child incidence (default is
  `NULL`).

- covariate_matrix:

  Optional covariate matrix used in incidence transformation (default is
  `NULL`).

## Value

A data frame containing model simulation results from `start_date` to
`end_date`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run simulation using max-posterior parameter set
sim_df <- simulate_with_max_posterior_params(
  results = results,
  start_date = "2021-01-01",
  end_date = "2021-12-31",
  model = data_sim
)
} # }
```
