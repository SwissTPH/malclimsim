# Extend Time-Varying Inputs Backwards in Time

This function extends selected time-varying inputs by repeating the
first year's data backwards for a specified number of years. This is
useful when simulating models that require initialization over an
extended period.

## Usage

``` r
extend_time_varying_inputs(
  param_inputs,
  days_per_year = 360,
  years_to_extend = 2
)
```

## Arguments

- param_inputs:

  A named list of model input vectors, where some elements may be
  time-varying (e.g., daily values).

- days_per_year:

  Integer specifying the number of time steps (e.g., days) per year.
  Defaults to 360.

- years_to_extend:

  Number of years to prepend to each time-varying input by repeating the
  first year's data. Defaults to 2.

## Value

A list with the same structure as `param_inputs`, where selected
time-varying inputs have been extended backward in time.

## Details

The following keys are treated as time-varying and extended:
`"cov_SMC"`, `"SMC"`, `"decay"`, `"c_R_D"`, `"temp"`. Other keys are
left unchanged.

## Examples

``` r
inputs <- list(
  cov_SMC = rep(0.5, 360 * 3),
  temp = rep(25, 360 * 3),
  constant = 42
)
extended <- extend_time_varying_inputs(inputs, days_per_year = 360, years_to_extend = 1)
length(extended$cov_SMC)  # should be 360 * 4
#> [1] 1440
```
