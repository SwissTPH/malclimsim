# Generate Model-Observation Comparison Function

Wraps a call to `generate_incidence_comparison()` to produce a
likelihood or loss function for comparing model output to observed data.

## Usage

``` r
generate_comparison_function(
  month,
  age_for_inf,
  incidence_observed,
  include_prev
)
```

## Arguments

- month:

  Logical; whether the model and data are aggregated by calendar month.

- age_for_inf:

  Character; specifies the age group used for inference (e.g., `"inc_C"`
  or `"inc_A"`).

- incidence_observed:

  A data frame of observed incidence data.

- include_prev:

  Logical; whether to include prevalence data in the comparison
  function.

## Value

A comparison function object that can be passed to
[`mcstate::pmcmc()`](https://rdrr.io/pkg/mcstate/man/pmcmc.html) or
other inference routines.
