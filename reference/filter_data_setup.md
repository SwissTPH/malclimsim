# Set Up Filtered Data for Likelihood Evaluation

Prepares the observed incidence data for comparison to simulated data.

## Usage

``` r
filter_data_setup(incidence_observed, month, initial_time_obs)
```

## Arguments

- incidence_observed:

  A data frame of observed incidence values.

- month:

  Logical; whether the data is aggregated by month.

- initial_time_obs:

  Integer; starting time index for aligning simulated and observed data.

## Value

A filtered data object suitable for use in likelihood comparison.
