# Plot Cases Averted Across SMC Scenarios (with 95% Credible Intervals)

Generates a bar plot of estimated cases averted for each SMC
pattern/scenario, including 95% credible interval error bars. Supports
either:

- A list of numeric vectors (posterior samples), or

- A list of lists with `mean`, `lower`, and `upper`.

## Usage

``` r
plot_cases_averted_barplot(
  estimates,
  title = "Estimated Cases Averted for Different SMC Timings",
  x_lab = "SMC Timing",
  y_lab = "Mean Cases Averted (95% CI)",
  horizontal = FALSE,
  by_year = FALSE,
  n_years = 1,
  bar_colors = c(`4 rounds (July start)` = "#0072B2", `4 rounds (June start)` =
    "#E69F00", `5 rounds (July start)` = "#009E73", `5 rounds (June start)` = "#D55E00")
)
```

## Arguments

- estimates:

  A named list where each element is either:

  - A numeric vector of posterior draws, or

  - A list with elements `mean`, `lower`, and `upper`.

- title:

  Character string for the main plot title.

- x_lab:

  Character string for the x-axis label.

- y_lab:

  Character string for the y-axis label.

- horizontal:

  Logical; if TRUE, produces a horizontal bar plot with patterns on the
  y-axis.

- by_year:

  Logical; if TRUE, divides all values by `n_years` to yield per-year
  estimates.

- n_years:

  Numeric; the number of years used for scaling if `by_year = TRUE`.
  Must be positive.

- bar_colors:

  Named character vector specifying fill colors for each SMC pattern.
  The names should match the names in `estimates`.

## Value

A `ggplot2` bar plot object showing mean cases averted with 95% CI.
