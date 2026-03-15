# Plot Default Priors

Plots the default or user-specified prior distributions for the
parameters of interest.

## Usage

``` r
plot_priors(param_inputs, proposal_matrix, params_to_estimate, priors = NULL)
```

## Arguments

- param_inputs:

  A list of model inputs (used only if default priors are initialized).

- proposal_matrix:

  A covariance matrix used to generate default priors if `priors` is not
  supplied.

- params_to_estimate:

  A character vector specifying which parameters to include in the plot.

- priors:

  A named list of priors. Each prior should be a list containing `min`,
  `max`, and a `prior` function that returns the log-density for a given
  numeric input. If `NULL`, default priors will be computed.

## Value

A `ggplot` object showing the prior distributions of the selected
parameters.
