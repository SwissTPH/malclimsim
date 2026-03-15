# Generate SMC Coverage Matrix for Simulation

Creates a monthly covariate matrix for SMC coverage from a deployment
pattern.

## Usage

``` r
generate_smc_coverage_matrix(
  start_date,
  end_date,
  years,
  month_pattern,
  avg_cov,
  exclude_years = NULL
)
```

## Arguments

- start_date:

  Simulation start date ("YYYY-MM-DD").

- end_date:

  Simulation end date ("YYYY-MM-DD").

- years:

  Numeric vector of years to apply the SMC schedule.

- month_pattern:

  A 12-length binary vector representing monthly SMC rounds.

- avg_cov:

  Numeric average SMC coverage to apply.

- exclude_years:

  Optional vector of years to exclude from coverage.

## Value

A data frame with `date_ymd` and `cov_SMC` columns.
