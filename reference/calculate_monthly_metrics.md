# Compute monthly metrics from SMC schedule

Compute monthly metrics from SMC schedule

## Usage

``` r
calculate_monthly_metrics(schedule, exclude_years = NULL)
```

## Arguments

- schedule:

  Data frame with SMC schedule (must include `dates`, `SMC`, `cov`,
  `decay`)

- exclude_years:

  Vector of years to exclude (e.g., c(2023))

## Value

Monthly summarized schedule with columns: month, SMC, cov, decay
