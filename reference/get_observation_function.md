# Get Observation Function Based on User Configuration

This helper function allows users to specify key options such as time
resolution, age group, inclusion of prevalence, and whether to model SMC
as a covariate using `beta_1`.

## Usage

``` r
get_observation_function(
  month = TRUE,
  age_for_inf = "u5",
  incidence_df,
  include_prev = TRUE,
  use_SMC_as_covariate = FALSE
)
```

## Arguments

- month:

  Logical. TRUE for monthly data, FALSE for weekly.

- age_for_inf:

  Character. One of: "total", "sep_ages", "u5", "o5", "all_ages".

- incidence_df:

  A data frame of observed incidence and prevalence values.

- include_prev:

  Logical. If TRUE, includes prevalence likelihood.

- use_SMC_as_covariate:

  Logical. If TRUE, uses `beta_1` covariate effect on incidence.

## Value

A likelihood function that takes (state, observed, pars) as arguments.
