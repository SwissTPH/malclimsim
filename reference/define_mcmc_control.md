# Define MCMC Control Settings

Creates two sets of MCMC control objects for use with
[`mcstate::pmcmc`](https://rdrr.io/pkg/mcstate/man/pmcmc.html), one with
a short run for adaptation, and one for the main run with full settings.

## Usage

``` r
define_mcmc_control(
  control_params,
  adaptive_params,
  save_trajectories,
  rerun_n,
  rerun_random
)
```

## Arguments

- control_params:

  A list with MCMC control settings including `n_steps`, `n_burnin`,
  `n_chains`, `n_workers`, and `n_threads_total`.

- adaptive_params:

  A list or object specifying the adaptive proposal settings (e.g., from
  `mcstate::adaptive_proposal()`).

- save_trajectories:

  Logical; whether to save state trajectories from the MCMC run.

- rerun_n:

  Integer; frequency of rerunning adaptive proposal fitting.

- rerun_random:

  Logical; whether rerun intervals are randomized.

## Value

A list with two elements:

- control1:

  An initial short `pmcmc_control` object (e.g., for proposal
  adaptation).

- control2:

  The full `pmcmc_control` object used for the main MCMC run.
