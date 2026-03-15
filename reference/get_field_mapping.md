# Get Model Output Field Names for Observation Function

Returns a named list of model state variables used for incidence
comparisons, based on the time scale and age group configuration.

## Usage

``` r
get_field_mapping(time, age_group)
```

## Arguments

- time:

  Character string. Either "month" or "week".

- age_group:

  Character. One of: "total", "sep_ages", "u5", "o5", "all_ages".

## Value

A named list with one or more of: mu, mu_C, mu_A, mu_C1, mu_C2, etc.
These correspond to expected incidence variables in the model output.
