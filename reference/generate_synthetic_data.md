# Generate Synthetic Incidence Data (Optional)

Optionally generates synthetic incidence data from a model, with the
ability to add Poisson noise.

## Usage

``` r
generate_synthetic_data(
  model,
  param_inputs,
  dates,
  month,
  month_unequal_days,
  noise,
  seed,
  synthetic,
  incidence_df
)
```

## Arguments

- model:

  A simulation model object used to generate synthetic data.

- param_inputs:

  A list of model parameters to use for simulation.

- dates:

  A vector of two dates (`start_date`, `end_date`) specifying the
  simulation period.

- month:

  Logical; whether to aggregate simulated data by calendar months.

- month_unequal_days:

  Logical; whether to weight months by unequal days during simulation.

- noise:

  Logical; if `TRUE`, adds Poisson noise to the simulated outputs.

- seed:

  Integer; random seed used for reproducible noise.

- synthetic:

  Logical; if `TRUE`, synthetic data will be generated. If `FALSE`,
  returns the input `incidence_df`.

- incidence_df:

  A data frame of observed incidence; returned unchanged if
  `synthetic = FALSE`.

## Value

A data frame of synthetic or observed incidence data with columns such
as `inc`, `inc_C`, `inc_A`.
