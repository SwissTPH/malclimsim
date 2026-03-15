# Generate an SMC deployment schedule with optional month-specific coverages

Generate an SMC deployment schedule with optional month-specific
coverages

## Usage

``` r
gen_smc_schedule(
  start_date,
  end_date,
  years,
  months_active,
  months_30_days = FALSE,
  coverage = 0.9,
  smc_day_of_month = 1
)
```

## Arguments

- start_date:

  Character. First date of the schedule ("YYYY-MM-DD").

- end_date:

  Character. Last date of the schedule ("YYYY-MM-DD").

- years:

  Integer vector. Years included in the schedule.

- months_active:

  Matrix (n_years x 12) with 1 = SMC month, 0 = not.

- months_30_days:

  Logical. Use 360-day calendar (every month = 30 days)?

- coverage:

  Either **(i)** a single numeric (legacy behaviour) or **(ii)** a
  numeric vector of length 12 *or* a named vector whose names are month
  numbers (1-12) / abbreviations ("Jan", "Feb", ...). Values give
  coverage for each month.

- smc_day_of_month:

  Integer. Day of month SMC is administered.

## Value

data.frame with columns:

- dates - sequence of dates in the period

- SMC - 1 on round day, 0 otherwise

- cov - coverage applied (piecewise-constant blocks)

- decay - efficacy decay computed from round days
