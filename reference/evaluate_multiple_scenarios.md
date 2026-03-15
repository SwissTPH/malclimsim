# Evaluate Multiple SMC Scenarios

Runs model simulations for a list of SMC coverage patterns and
summarizes outcomes.

## Usage

``` r
evaluate_multiple_scenarios(
  patterns,
  smc_day_of_month = 1,
  model,
  param_inputs,
  param_samples,
  start_date,
  end_date,
  avg_cov,
  years,
  exclude_years = 2023,
  mu_transform_C = NULL,
  mu_transform_A = NULL,
  outcome_fn = function(y1, y0) sum(y1$inc_C_transformed),
  o1 = NULL,
  ci_level = 0.95,
  out_dir = NULL,
  month = FALSE,
  apply_decay = TRUE,
  use_SMC_as_covariate = FALSE,
  noise = FALSE
)
```

## Arguments

- patterns:

  A named list of 12-element binary vectors or lists of such vectors
  (per year), indicating months when SMC is applied.

- smc_day_of_month:

  Integer specifying the day of the month to begin SMC each round
  (default is 1).

- model:

  A compiled model object to simulate from.

- param_inputs:

  A named list of baseline model parameters.

- param_samples:

  A matrix or data frame of sampled parameter sets (one row per sample).

- start_date:

  Simulation start date (character or Date).

- end_date:

  Simulation end date (character or Date).

- avg_cov:

  Average SMC coverage to apply during active months.

- years:

  A vector of years to apply SMC coverage (e.g., 2018:2023).

- exclude_years:

  Years to exclude when computing summaries (default is 2023).

- mu_transform_C:

  Optional transformation function for compartment C incidence.

- mu_transform_A:

  Optional transformation function for compartment A incidence.

- outcome_fn:

  Function computing the scalar outcome of interest from two simulation
  outputs (default: sum of transformed incidence).

- o1:

  Optional baseline simulation for comparison.

- ci_level:

  Confidence level for summary intervals (default is 0.95).

- out_dir:

  Output directory to save plots. If NULL, plots are not saved.

- month:

  Logical. Whether the model and summaries are in monthly (TRUE) or
  weekly (FALSE) time (default is FALSE).

- apply_decay:

  Logical. Whether to apply time-decay to SMC coverage values (default
  is TRUE).

- use_SMC_as_covariate:

  Logical. Whether SMC is used as a covariate in the observation model
  instead of being built into the transmission model (default is FALSE).

- noise:

  Logical. Whether to add stochastic noise to simulated incidence
  (default is FALSE).

## Value

A list with three elements:

- outputs:

  A named list of lists with estimates, plots, and summaries for each
  scenario.

- summaries:

  A named list of time series summaries for each scenario.

- estimates:

  A named list of scalar outcome estimates for each scenario.
