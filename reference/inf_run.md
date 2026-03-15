# Run MCMC Inference Simulation

Generates synthetic incidence data if required, sets up MCMC parameters
and priors, defines model comparison functions, and runs the MCMC
process.

## Usage

``` r
inf_run(
  model,
  param_inputs,
  control_params,
  params_to_estimate,
  proposal_matrix,
  adaptive_params,
  start_values,
  noise = FALSE,
  seed = 24,
  month_unequal_days = FALSE,
  dates,
  synthetic = TRUE,
  incidence_df = NULL,
  save_trajectories = TRUE,
  rerun_n = Inf,
  rerun_random = FALSE,
  param_priors = NULL,
  override_priors = NULL,
  n_years_warmup = 3,
  obs_config
)
```

## Arguments

- model:

  The model object used for simulation and inference.

- param_inputs:

  List of input parameters for the model, including temperature, decay,
  and coverage details.

- control_params:

  Named list of control parameters for MCMC (e.g., number of steps,
  burn-in, chains).

- params_to_estimate:

  Character vector of parameter names to be estimated.

- proposal_matrix:

  Covariance matrix used for the proposal distribution in MCMC.

- adaptive_params:

  List of parameters controlling adaptation in MCMC (e.g., update
  interval, scale factor).

- start_values:

  Named vector of initial values for MCMC chains.

- noise:

  Logical; if TRUE, adds noise to the synthetic incidence data.

- seed:

  Numeric; random seed for reproducibility of simulations.

- month_unequal_days:

  Logical; if TRUE, accounts for different numbers of days in each
  month.

- dates:

  Character or Date vector of length two specifying the start and end
  date of the simulation period.

- synthetic:

  Logical; if TRUE, synthetic data is generated using the model and
  parameters.

- incidence_df:

  Data frame of observed incidence data. Required if
  `synthetic = FALSE`.

- save_trajectories:

  Logical; if TRUE, saves intermediate state trajectories during MCMC.

- rerun_n:

  Numeric; frequency (in iterations) at which to re-run the
  deterministic model (default is Inf).

- rerun_random:

  Logical; if TRUE, re-run times are randomized instead of fixed
  interval.

- param_priors:

  Optional named list of prior distributions for parameters being
  estimated.

- override_priors:

  Optional list to override specific priors in `param_priors`.

- n_years_warmup:

  Number of years to run the model prior to the start date (default is
  3).

- obs_config:

  A named list specifying configuration of the observation model. Use
  [`make_obs_config()`](https://swisstph.github.io/malclimsim/reference/make_obs_config.md)
  to construct this.

## Value

A list containing MCMC results, including posterior samples, fixed
parameters, priors, and incidence data.
