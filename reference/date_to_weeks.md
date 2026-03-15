# Convert Date Range to Weekly Format (Sundays)

This function takes a start date and an end date, and returns a vector
of Dates representing the Sunday of each week in the range.

## Usage

``` r
date_to_weeks(start_date, end_date)
```

## Arguments

- start_date:

  The start date as a character string or `Date` object.

- end_date:

  The end date as a character string or `Date` object.

## Value

A vector of `Date` objects representing each Sunday between the dates.

## Examples

``` r
date_to_weeks("2022-01-01", "2022-03-01")
#>  [1] "2021-12-26" "2022-01-02" "2022-01-09" "2022-01-16" "2022-01-23"
#>  [6] "2022-01-30" "2022-02-06" "2022-02-13" "2022-02-20" "2022-02-27"
```
