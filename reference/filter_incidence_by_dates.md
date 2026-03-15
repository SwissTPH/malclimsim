# Filter Incidence Data by Date Range

Filters a time series of incidence data to only include rows within a
specified date range.

## Usage

``` r
filter_incidence_by_dates(incidence_df, dates)
```

## Arguments

- incidence_df:

  A data frame with a `date_ymd` column (must be coercible to Date).

- dates:

  A vector of two dates (start and end) as `Date` objects or strings in
  `"YYYY-MM-DD"` format.

## Value

A filtered data frame containing only rows with `date_ymd` within the
specified range.

## Examples

``` r
df <- data.frame(date_ymd = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 100),
                 incidence = rpois(100, 5))
filtered <- filter_incidence_by_dates(df, c("2020-01-10", "2020-01-20"))
```
