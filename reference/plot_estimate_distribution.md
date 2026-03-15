# Plot Histogram of Estimated Effects with Confidence Interval and Trimmed X-Axis

Plot Histogram of Estimated Effects with Confidence Interval and Trimmed
X-Axis

## Usage

``` r
plot_estimate_distribution(
  estimates,
  x_label = "Cases Averted",
  title = "Posterior Distribution of Estimated Cases Averted",
  bins = 30,
  ci_level = 0.95
)
```

## Arguments

- estimates:

  Numeric vector of estimates (e.g., cases averted).

- x_label:

  Label for the x-axis.

- title:

  Title for the plot.

- bins:

  Number of histogram bins (default = 30).

- ci_level:

  Confidence interval level (default = 0.95).

## Value

A ggplot object showing histogram, density, and CI lines.
