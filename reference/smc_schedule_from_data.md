# Generate SMC Coverage and Decay Schedule

Constructs a daily schedule of SMC (Seasonal Malaria Chemoprevention)
coverage and decay values over a specified time period. Missing values
are filled forward using the last observed non-zero coverage, and decay
is applied according to a specified exponential decay constant.

## Usage

``` r
smc_schedule_from_data(smc_cov, months_30_days, years, const = -0.1806)
```

## Arguments

- smc_cov:

  A data frame with at least two columns:

  - `date_start`: A `Date` or character vector indicating the start of
    coverage.

  - `coverage`: Numeric values indicating the coverage on each date.

- months_30_days:

  Logical. If `TRUE`, generates a synthetic 360-day calendar (12 months
  of 30 days). Defaults to `FALSE`.

- years:

  A numeric vector of length 2 specifying the range of years to include
  (e.g., `c(2014, 2022)`).

- const:

  A numeric scalar for the exponential decay constant (default:
  `-0.1806`).

## Value

A data frame with the following columns:

- `dates`: A daily sequence of dates across the specified time range.

- `SMC`: Binary indicator (1 if SMC applied on that day, 0 otherwise).

- `cov`: Propagated SMC coverage values after filling and masking
  non-SMC months.

- `decay`: Daily decay values calculated using the specified constant.

## Examples

``` r
if (FALSE) { # \dontrun{
smc_data <- data.frame(
  date_start = as.Date(c("2014-06-01", "2014-07-01", "2015-06-01")),
  coverage = c(0.8, 0.85, 0.9)
)
schedule <- smc_schedule_from_data(smc_cov = smc_data,
months_30_days = FALSE, years = c(2014, 2015))
head(schedule)
} # }
```
