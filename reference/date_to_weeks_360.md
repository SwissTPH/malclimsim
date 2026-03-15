# Generate weekly model dates (by 7-day blocks) in 360-day calendar

Generate weekly model dates (by 7-day blocks) in 360-day calendar

## Usage

``` r
date_to_weeks_360(start_year, n_days)
```

## Arguments

- start_year:

  integer, year of day

- n_days:

  integer, total days (must be divisible by 7)

## Value

A data.frame with one row per week: year, month, day (all numeric)
