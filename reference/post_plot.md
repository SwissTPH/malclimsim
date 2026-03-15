# Plot Posterior Distributions of Estimated Parameters

Generates density or histogram plots for each estimated parameter from
MCMC output, with annotated quantiles and optional true values if
available. It helps users visualize the posterior distribution of model
parameters.

## Usage

``` r
post_plot(
  results_list,
  params_to_estimate,
  dim_plot,
  show_true = TRUE,
  true_value = NULL,
  show_prior = FALSE,
  prior_n = 10000,
  title = "",
  run_labels = NULL,
  plot_type = "histogram"
)
```

## Arguments

- results_list:

  A list containing MCMC results from multiple runs. Each element should
  contain posterior samples.

- params_to_estimate:

  A character vector of parameter names to include in the plot.

- dim_plot:

  A numeric vector of length 2 specifying the number of rows and columns
  for arranging plots.

- show_true:

  Logical; if `TRUE`, true parameter values will be indicated on the
  plots if provided.

- true_value:

  A named numeric vector of true values for parameters (optional). Must
  match `params_to_estimate` if used.

- show_prior:

  Logical; if `TRUE`, prior distributions will be overlaid on the plots.

- prior_n:

  Integer; the number of points to sample from the prior distribution.

- title:

  A character string specifying the title for the entire plot layout
  (optional).

- run_labels:

  A character vector specifying labels for each MCMC run (optional).
  Must match the length of `results_list`.

- plot_type:

  A character string specifying the type of plot: `"histogram"` or
  `"density"`.

## Value

A combined plot object displaying posterior distributions for each
parameter. The plot includes optional lines for true values if
`show_true = TRUE`, and prior distributions if `show_prior = TRUE`.
