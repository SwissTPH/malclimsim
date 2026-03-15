# Plot Posterior Predictive Check (PPC)

Plot Posterior Predictive Check (PPC)

## Usage

``` r
plot_ppc(
  ppc_data,
  ci_data = NULL,
  title = "Posterior Predictive Check",
  xlab = "Date",
  ylab = "Value",
  show_obs = TRUE,
  show_best = TRUE,
  show_ribbon = TRUE
)
```

## Arguments

- ppc_data:

  Data frame from
  [`prepare_ppc_data()`](https://swisstph.github.io/malclimsim/reference/prepare_ppc_data.md)
  with columns: `date_ymd`, `value`, `variable`, and `label`.

- ci_data:

  Optional data frame with columns: `date_ymd`, `lower`, `upper`,
  `variable`.

- title:

  Plot title.

- xlab:

  Label for x-axis.

- ylab:

  Label for y-axis.

- show_obs:

  Logical. Whether to show observed data.

- show_best:

  Logical. Whether to show deterministic prediction.

- show_ribbon:

  Logical. Whether to show posterior CI ribbon.

## Value

A `ggplot2` object.
