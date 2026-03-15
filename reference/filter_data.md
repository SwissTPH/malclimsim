# Relate Time Steps to Observed Data for mcstate

This function formats observed data into a structure usable by the
`mcstate` particle filter, adjusting for time steps (weekly or monthly)
and ensuring compatibility with the filtering framework.

## Usage

``` r
filter_data(incidence_observed, month = FALSE, initial_time_obs = 0, rate = 7)
```

## Arguments

- incidence_observed:

  A data frame of observed incidence data. Must include a column named
  `week_no` (if `month = FALSE`) or `month_no` (if `month = TRUE`), and
  a column named `value` for the observed incidence.

- month:

  Logical. If `TRUE`, assumes data is monthly; otherwise, weekly.

- initial_time_obs:

  Integer. The starting time index for monthly data (default is 0).

- rate:

  Integer. The time step rate (e.g., 7 for weekly, 30 for monthly).
  Default is 7.

## Value

A data frame formatted for use with
[`mcstate::particle_filter`](https://rdrr.io/pkg/mcstate/man/particle_filter.html),
including time start and end steps.

## Examples

``` r
if (requireNamespace("mcstate", quietly = TRUE)) {
  # Weekly example
  incidence_observed <- data.frame(
    week_no = 0:10,
    value = rpois(11, lambda = 5)
  )
  formatted_data <- filter_data(incidence_observed, month = FALSE, initial_time_obs = 0, rate = 7)
  head(formatted_data)

  # Monthly example
  incidence_observed_month <- data.frame(
    month_no = 0:10,
    value = rpois(11, lambda = 10)
  )
  formatted_month_data <- filter_data(incidence_observed_month, month = TRUE,
  initial_time_obs = 0, rate = 30)
  head(formatted_month_data)
}
#>   month_no_start month_no_end time_start time_end value
#> 2              1            2          0       30    16
#> 3              2            3         30       60    15
#> 4              3            4         60       90    10
#> 5              4            5         90      120     6
#> 6              5            6        120      150     8
#> 7              6            7        150      180    11
```
