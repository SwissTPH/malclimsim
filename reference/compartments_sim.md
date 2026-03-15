# Simulate Compartments Over Time with Prewarm Period

This function simulates the compartments of an epidemiological model
over a specified time period, including an optional prewarm period to
allow the model to reach equilibrium before the main simulation period.

## Usage

``` r
compartments_sim(
  model,
  param_inputs,
  start_date,
  end_date,
  prewarm_years = 2,
  days_per_year = 360
)
```

## Arguments

- model:

  An epidemiological model object used for simulation.

- param_inputs:

  A named list or vector of model parameters to be used in the
  simulation.

- start_date:

  A character string or Date object specifying the start date of the
  main simulation (e.g., "2014-01-01").

- end_date:

  A character string or Date object specifying the end date of the
  simulation (e.g., "2022-12-31").

- prewarm_years:

  Integer specifying the number of years to prewarm the model before the
  main simulation (default: 2).

- days_per_year:

  Integer specifying the number of days in a model year (default: 360).

## Value

A data frame containing the simulated values of each compartment (SC,
EC, IC, etc.) over the simulation period with a weekly time resolution.

## Details

The function simulates the dynamics of susceptible, exposed, infected,
treated, and recovered compartments for both children and adults. It
first runs a prewarm period (if specified) to allow the model to
stabilize, before simulating the main study period. The results are
aggregated at weekly intervals.
