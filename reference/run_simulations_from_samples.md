# Run Simulations from Sampled Parameter Sets

This function runs a simulation model multiple times, each time using a
different row from a sampled parameter matrix (e.g., from posterior
draws). It is used to generate predictive simulations or counterfactual
scenarios based on parameter uncertainty. Each simulation returns output
from the
[`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md)
function, which integrates optional transformations and covariates.

## Usage

``` r
run_simulations_from_samples(
  model,
  param_inputs,
  param_samples,
  start_date,
  end_date,
  prewarm_years = 2,
  mu_transform_C = NULL,
  mu_transform_A = NULL,
  covariate_matrix = NULL,
  month = TRUE,
  noise = FALSE
)
```

## Arguments

- model:

  A simulation model function (e.g., a compartmental ODE model) that
  takes in a parameter list and returns state trajectories.

- param_inputs:

  A named list of baseline input parameters. These serve as the default
  values before updating with each row from `param_samples`.

- param_samples:

  A matrix or data frame of sampled parameter values, where each row
  corresponds to a unique parameter set to simulate.

- start_date, end_date:

  The simulation date range as `Date` or character objects coercible to
  `Date`. Typically, these define the time window for prediction after
  prewarming.

- prewarm_years:

  Number of years before `start_date` used for prewarming the model
  (e.g., to reach equilibrium or load past covariates). Defaults to 2.

- mu_transform_C:

  Optional function to transform the underlying `mu` incidence into a
  predicted count (e.g., logistic or log link for cluster C).

- mu_transform_A:

  Optional transformation function for cluster A (if used).

- covariate_matrix:

  Optional covariate matrix used for transforming incidence (e.g., for
  climate-driven models or spatial heterogeneity).

- month:

  Logical. If `TRUE`, the output will be aggregated or indexed monthly
  (rather than by week or day).

- noise:

  Logical. If `TRUE`, draws the simulated incidence from a negative
  binomial distribution using `size_1` as the dispersion parameter.

## Value

A list of data frames, each containing a simulation run with the same
structure as the output of
[`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md).
The list has one element per row in `param_samples`.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_list <- run_simulations_from_samples(
  model = malaria_model,
  param_inputs = base_params,
  param_samples = posterior_draws,
  start_date = "2016-01-01",
  end_date = "2022-12-31",
  mu_transform_C = logit_transform,
  covariate_matrix = climate_covariates,
  noise = TRUE
)
} # }
```
