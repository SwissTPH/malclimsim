# Update Default Priors

Updates the default prior specifications for a set of parameters,
allowing users to replace or modify elements such as distribution type,
bounds, and initial values.

## Usage

``` r
update_priors(
  param_inputs,
  proposal_matrix,
  params_to_estimate,
  priors = NULL,
  new_priors
)
```

## Arguments

- param_inputs:

  A named list of model parameter inputs used to initialize default
  priors.

- proposal_matrix:

  A proposal covariance matrix used to set up default priors.

- params_to_estimate:

  A character vector of parameters to estimate and update priors for.

- priors:

  Optional. A list of existing priors. If NULL (default), default priors
  are generated using `initialize_priors()`.

- new_priors:

  A named list of new prior specifications. Each element should be a
  list with updated fields (e.g., `initial`, `min`, `max`, `prior`).

## Value

A list of prior specifications with selected entries updated based on
`new_priors`.
