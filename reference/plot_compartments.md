# Plot Compartmental Simulation Results

This function visualizes the output of a compartmental model simulation
with options for line plots, stacked area plots, and faceted
visualizations. Users can also plot proportions relative to the total
population.

## Usage

``` r
plot_compartments(
  compart_df,
  plot_type = "line",
  compartments = NULL,
  plot_proportion = FALSE,
  log_scale = FALSE,
  title = "Epidemic Compartments Over Time",
  color_palette = "Set2",
  y_label = "Population Count"
)
```

## Arguments

- compart_df:

  Data frame containing compartments simulation output.

- plot_type:

  Character string specifying the type of plot. Options are "line",
  "stacked", or "facet".

- compartments:

  Optional character vector of compartments to plot (e.g., c("SC",
  "IC")).

- plot_proportion:

  Logical indicating whether to plot proportions relative to the total
  population (`P` column).

- log_scale:

  Logical indicating whether to apply a logarithmic scale to the y-axis.

- title:

  Character string specifying the plot title.

- color_palette:

  Character string specifying the color palette from RColorBrewer
  (default: "Set2").

- y_label:

  Character string specifying the Y-axis label (default: "Population
  Count").

## Value

A ggplot object visualizing the compartmental data.
