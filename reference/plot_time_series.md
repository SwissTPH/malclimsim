# Plot Time Series of Malaria Incidence (Raw and Transformed)

This function creates time series plots for malaria incidence,
optionally including transformed incidence values (e.g., for observation
models). Also supports plotting climate covariates, with faceting and
overlay options.

## Usage

``` r
plot_time_series(
  results,
  met = NULL,
  plot_title = "Malaria Incidence Time Series",
  incidence_y_label = "Malaria Incidence",
  climate_y_label = "Climate",
  climate_facet = FALSE,
  select_incidence = c(">=5", "<5", "total"),
  select_climate = c("temp", "cumrain"),
  incidence_colors = c(`>=5` = "blue", `<5` = "red", total = "green", `transformed_>=5` =
    "skyblue", `transformed_<5` = "orange", transformed_total = "darkgreen"),
  climate_colors = c(temp = "orange", cumrain = "purple"),
  climate_alpha = 0.7,
  base_size = 15
)
```

## Arguments

- results:

  A data frame from
  [`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md)
  containing columns like inc_C, inc_A, inc_C_transformed,
  inc_A_transformed, and date_ymd.

- met:

  Optional data frame of climate covariates.

- plot_title:

  Title of the plot.

- incidence_y_label:

  Y-axis label for incidence.

- climate_y_label:

  Y-axis label for climate plot.

- climate_facet:

  Whether to facet incidence and climate plots.

- select_incidence:

  Vector of incidence types to include. Options: `"<5"`, `">=5"`,
  `"total"`, `"transformed_<5"`, `"transformed_>=5"`,
  `"transformed_total"`.

- select_climate:

  Climate vars to include (e.g., `"temp"`, `"cumrain"`).

- incidence_colors:

  Named vector of colors for incidence types.

- climate_colors:

  Named vector of colors for climate vars.

- climate_alpha:

  Transparency for climate lines.

- base_size:

  Base font size.

## Value

A ggplot or grid of ggplots
