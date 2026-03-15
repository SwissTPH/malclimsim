# Create Adaptive Proposal and MCMC Control Parameters

This function defines the default values for adaptive proposal control
parameters and MCMC control parameters, while allowing users to specify
custom values.

## Usage

``` r
create_mcmc_params(
  stage = "stage1",
  initial_vcv_weight = 1,
  initial_scaling = 2,
  initial_scaling_weight = NULL,
  min_scaling = 0.1,
  scaling_increment = NULL,
  log_scaling_update = TRUE,
  acceptance_target = 0.234,
  forget_rate = 0.6,
  forget_end = Inf,
  adapt_end = Inf,
  pre_diminish = 40000,
  n_steps = 10000,
  n_burnin = 0,
  n_chains = 4,
  n_workers = 4,
  n_threads_total = 8
)
```

## Arguments

- stage:

  possible values are "NULL", stage1", "stage2", or "noadapt"

- initial_vcv_weight:

  Weight for the initial variance-covariance matrix (default = 1).

- initial_scaling:

  Scaling factor for the proposal (default = 2).

- initial_scaling_weight:

  Optional weight for initial scaling (default = NULL).

- min_scaling:

  Minimum scaling factor for the proposal (default = 0.1).

- scaling_increment:

  Increment for scaling factor (default = NULL).

- log_scaling_update:

  Logical, should the scaling be updated on a log scale? (default =
  TRUE).

- acceptance_target:

  The target acceptance rate for the proposal (default = 0.234).

- forget_rate:

  Rate at which the proposal forgets past history (default = 0.6).

- forget_end:

  Time step at which the forgetting stops (default = Inf).

- adapt_end:

  Time step at which the adaptation stops (default = Inf).

- pre_diminish:

  Time steps before the diminishing starts (default = 40000).

- n_steps:

  Total number of MCMC steps (default = 10000).

- n_burnin:

  Number of burn-in steps (default = 0).

- n_chains:

  Number of MCMC chains (default = 4).

- n_workers:

  Number of workers for parallel execution (default = 4).

- n_threads_total:

  Total number of threads for parallel execution (default = 8).

## Value

A list containing the adaptive proposal parameters and MCMC control
parameters.
