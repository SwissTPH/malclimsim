# Impute Missing Climate Data Using Climatological Daily Means

This function imputes missing dates in a climate time series using the
average value for each (month, day) combination across all prior years
in the dataset. It is useful for filling in short gaps at the end of a
time series, such as imputing values for December 2023 using earlier
years' averages.

## Usage

``` r
impute_climate_to_end_date(data, end_date)
```

## Arguments

- data:

  A data frame with at least the following columns:

  - `dates` (Date): The observation date.

  - `CumulativeRainfall`, `anom`, `temp`: Numeric climate variables to
    impute.

- end_date:

  A `Date` object indicating the final date that should be present. The
  function will impute values from the last available date up to
  `end_date`. Assumes the `dates` column ends before this value.

## Value

A data frame with the original data and appended imputed values, sorted
by `dates`.

## Examples

``` r
# Suppose your climate data ends at 2023-12-01
# met_360_completed <- impute_climate_to_end_date(met_360, as.Date("2023-12-31"))
```
