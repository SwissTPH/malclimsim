# Run Model Simulation with Pre-Warm Period

Extends time-varying covariate inputs and runs the model over a pre-warm
period to allow system stabilization before the target simulation
window.

## Usage

``` r
run_simulation_with_prewarm(
  model,
  param_inputs,
  dates,
  prewarm_years = 2,
  days_per_year = 360
)
```

## Arguments

- model:

  A model object used for simulation (compatible with
  [`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md)).

- param_inputs:

  A list of model parameters, including time-varying covariates.

- dates:

  A vector of two dates: start and end of the target simulation window
  (as `Date` or string in `"YYYY-MM-DD"` format).

- prewarm_years:

  Integer; number of years to pre-warm the model before the actual
  simulation window. Defaults to 2.

- days_per_year:

  Integer; number of time steps per year (typically 360 or 365).
  Defaults to 360.

## Value

A data frame of simulated results filtered to include only dates from
the original `start_date` onward.

## Details

Time-varying parameters are extended backward in time by repeating the
first year's values for `prewarm_years`. The simulation is then run over
this extended window, and the pre-warm results are discarded to avoid
initialization artifacts.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_results <- run_simulation_with_prewarm(
  model = malaria_model,
  param_inputs = parameters,
  dates = c("2020-01-01", "2022-12-31"),
  prewarm_years = 2
)
} # }
```
