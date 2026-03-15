# Compute weekly metrics from SMC schedule (epidemiological weeks, starting Sunday)

Compute weekly metrics from SMC schedule (epidemiological weeks,
starting Sunday)

## Usage

``` r
calculate_weekly_metrics(schedule, exclude_years = NULL)
```

## Arguments

- schedule:

  Data frame with SMC schedule (must include `dates`, `SMC`, `cov`,
  `decay`)

- exclude_years:

  Vector of years to exclude (e.g., c(2023))

## Value

Weekly summarized schedule with columns: week_start (Date), SMC, cov,
decay
