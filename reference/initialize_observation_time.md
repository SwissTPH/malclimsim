# Initialize Observation Time Index for Incidence Data

Aligns simulated and observed data by setting an initial time index
(`month_no`) for observations.

## Usage

``` r
initialize_observation_time(simulated_result, incidence_df)
```

## Arguments

- simulated_result:

  A data frame from simulation output, with a `date_ymd` column.

- incidence_df:

  A data frame of observed incidence data, with a `date_ymd` column.

## Value

The integer `0`, and updates the `month_no` column of `incidence_df`
in-place (used for alignment).
