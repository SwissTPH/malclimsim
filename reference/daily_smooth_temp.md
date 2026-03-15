# Daily Smoothed Temperature Function

This function takes a data frame with monthly temperature data and
produces a daily smoothed temperature series using a smoothing spline.
The output includes predicted temperatures for each day from the
earliest to the latest date in the original data.

## Usage

``` r
daily_smooth_temp(temp_df)
```

## Arguments

- temp_df:

  A data frame with two columns:

  - `Date`: A character vector in "YYYY-MM" format (e.g., "2024-11") or
    a Date object representing the first day of the month.

  - `Temperature`: A numeric vector giving the temperature for each
    month.

## Value

A data frame with two columns:

- `Date`: A Date vector representing daily dates.

- `Temperature`: A numeric vector with the smoothed daily temperature
  predictions.
