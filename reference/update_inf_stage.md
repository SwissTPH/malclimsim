# Update proposal matrix, start values, and MCMC settings for next stage

Update proposal matrix, start values, and MCMC settings for next stage

## Usage

``` r
update_inf_stage(
  results_obj,
  proposal_matrix,
  param_names,
  S_prev = 3000,
  draw_n = 4,
  shrink = 0.8,
  stage,
  n_steps = NULL
)
```

## Arguments

- results_obj:

  Result object returned by inf_run (with MCMC traces)

- proposal_matrix:

  Current proposal variance matrix

- param_names:

  Character vector of parameter names (rownames of proposal_matrix)

- S_prev:

  Integer: number of samples to use for variance-covariance extraction

- draw_n:

  Integer: number of random draws for new start values per chain

- shrink:

  Numeric: shrinkage factor for bounds (e.g. 0.8 retains 80% of
  interval)

- stage:

  Character: which stage to pull from create_mcmc_params()

- n_steps:

  Integer or NULL: if not NULL, overrides control_params\$n_steps

## Value

A list with elements:

- proposal_matrix : updated proposal variance matrix

- start_values : matrix of new start values (chains x parameters)

- adaptive_params : list of adaptive settings for this stage

- control_params : list of control settings for this stage (with n_steps
  possibly overridden)
