# Extend Time-Varying Inputs to Match Simulation Duration

Ensures that all relevant time-varying parameters (e.g., temperature,
SMC coverage) are extended to match the required number of simulation
days by prepending repeated values from the first year of data.

## Usage

``` r
extend_time_varying_inputs_to_length(
  param_inputs,
  days_needed,
  days_per_year = 360
)
```

## Arguments

- param_inputs:

  A named list of model input vectors. Some entries should be
  time-varying inputs.

- days_needed:

  Integer specifying the total number of days required for simulation
  (e.g., `n_days`).

- days_per_year:

  Integer number of days per model year (default is 360).

## Value

A list of parameter inputs where time-varying entries have been extended
(and trimmed) to exactly `days_needed` days.

## Details

The following keys are treated as time-varying: `"cov_SMC"`, `"SMC"`,
`"decay"`, `"c_R_D"`, `"temp"`. These are extended by repeatedly
prepending the first year of values until the total length is
sufficient. Non-time-varying keys are returned unchanged.

## Examples

``` r
inputs <- list(temp = rep(25, 360), constant = 0.1)
extended <- extend_time_varying_inputs_to_length(inputs, days_needed = 720)
length(extended$temp)  # should be 720
#> [1] 720
```
