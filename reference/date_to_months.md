# Convert Date Range to Monthly Format

This function takes a start date and an end date, and converts the range
into a vector of the first day of each month in the range, formatted as
"YYYY-MM".

## Usage

``` r
date_to_months(start_date, end_date)
```

## Arguments

- start_date:

  The start date as a character string or `Date` object.

- end_date:

  The end date as a character string or `Date` object.

## Value

A vector of `Date` objects representing the first day of each month
within the given date range.

## Examples

``` r
date_to_months("2022-01-01", "2022-06-01")
#> [1] "2022-01-01" "2022-02-01" "2022-03-01" "2022-04-01" "2022-05-01"
#> [6] "2022-06-01"
```
