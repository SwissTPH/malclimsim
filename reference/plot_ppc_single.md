# Plot Single Time Series Comparison (Uncomplicated or Severe Cases)

Displays time series with optional credible intervals, observed points,
and consistent legend formatting across scenarios (e.g., With/Without
SMC, June/July 2023, etc.)

## Usage

``` r
plot_ppc_single(
  plot_data,
  ci_data = NULL,
  obs_data = NULL,
  obs_col = NULL,
  plot_title = "Time Series Comparison",
  xlim = NULL,
  ylim = NULL,
  severity = "Uncomplicated",
  scale_severe = 1,
  scale_severe_by_year = NULL
)
```

## Arguments

- plot_data:

  Data frame with columns: date_ymd (Date), value (numeric), label

- ci_data:

  Optional data frame: date_ymd, lower, upper, label

- obs_data:

  Optional observed data (must have date_ymd)

- obs_col:

  Name of obs column in obs_data

- plot_title:

  Title

- xlim:

  x‐axis limits (Date vector)

- ylim:

  y‐axis limits (numeric vector)

- severity:

  Either "Uncomplicated" or "Severe"

- scale_severe:

  Scaling factor for severe cases (used if no per-year mapping is
  provided)

- scale_severe_by_year:

  Named vector of multipliers by year (optional)

## Value

ggplot2 plot object
