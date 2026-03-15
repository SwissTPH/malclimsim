# Create Starting Values for Parameters

This function generates starting values for parameters based on the
specified method. It can either draw random values from a uniform
distribution or assign values from a specified interval.

## Usage

``` r
create_start_values(
  params_to_estimate,
  control_params,
  min_max_start_values = NULL,
  random = TRUE,
  seed = 10,
  model,
  param_inputs
)
```

## Arguments

- params_to_estimate:

  A character vector of parameter names that are to be estimated.

- control_params:

  A list of control parameters (must include `n_chains` to specify the
  number of chains).

- min_max_start_values:

  A named list with minimum and maximum start values for each parameter.

- random:

  A logical value. If `TRUE`, random values are drawn from a uniform
  distribution; if `FALSE`, values are evenly spaced between the lower
  and upper bounds.

- seed:

  An integer specifying the seed for random number generation (default
  is 10).

- model:

  A model object with a `param()` method that returns all parameter
  names.

- param_inputs:

  A list of parameter values passed to initialize the model, used to
  extract parameters.

## Value

A matrix of starting values for each parameter and chain.
