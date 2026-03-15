# Create Observation Function Configuration

Constructs a configuration list specifying how to compare observed and
modeled data. This configuration is passed to the observation likelihood
generator function (`generate_incidence_comparison()`) and determines
which outcomes are compared and how.

## Usage

``` r
make_obs_config(
  use_monthly = TRUE,
  age_group = "u5",
  include_prev = FALSE,
  use_SMC_as_covariate = FALSE,
  log_link = FALSE,
  include_pop_growth = FALSE
)
```

## Arguments

- use_monthly:

  Logical. If `TRUE`, observation and model outputs are matched at
  monthly resolution; otherwise, weekly data is assumed.

- age_group:

  Character. One of `"u5"`, `"o5"`, `"sep_ages"`, `"total"`, or
  `"all_ages"`. This determines which age group(s) are included in the
  likelihood.

- include_prev:

  Logical. If `TRUE`, the observation function includes a prevalence
  likelihood term.

- use_SMC_as_covariate:

  Logical. If `TRUE`, the incidence is adjusted using a regression
  effect (`beta_1`) on SMC coverage.

- log_link:

  Logical. If `TRUE`, a log-link is used to transform the predicted
  incidence when adjusting by SMC coverage. If `FALSE`, a linear form is
  used.

- include_pop_growth:

  Logical. If `TRUE`, model-predicted incidence is rescaled to account
  for population growth using the `r_C` and/or `r_A` outputs from the
  model.

## Value

A named list containing configuration options to pass to the observation
likelihood function.

## Examples

``` r
obs_config <- make_obs_config(
  use_monthly = TRUE,
  age_group = "u5",
  include_prev = FALSE,
  use_SMC_as_covariate = TRUE,
  log_link = TRUE,
  include_pop_growth = TRUE
)
```
