# MCMC Diagnostic and Summary Plots with Save Options

Same functionality as before, now with file saving capabilities for
diagnostics.

## Usage

``` r
MCMC_diag(
  results,
  params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"),
  thin = 1,
  save = FALSE,
  output_dir = "mcmc_diagnostics",
  file_prefix = "mcmc",
  plot_width = 7,
  plot_height = 5
)
```

## Arguments

- results:

  MCMC result object.

- params:

  Diagnostics to compute.

- thin:

  Thinning interval.

- save:

  Logical; whether to save output (default: FALSE).

- output_dir:

  Directory to save outputs (default: "mcmc_diagnostics").

- file_prefix:

  Optional prefix for filenames.

- plot_width:

  Width of plots in inches (default: 7).

- plot_height:

  Height of plots in inches (default: 5).

## Value

Named list of diagnostics.
