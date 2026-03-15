# Create a Proposal Matrix

Generates a proposal covariance matrix for use in MCMC algorithms. The
matrix is initialized with specified variances for parameters to be
estimated and optionally incorporates parameter correlations.

## Usage

``` r
create_proposal_matrix(
  params_to_estimate,
  proposal_variance = NULL,
  correlation_matrix = NULL,
  model,
  param_inputs
)
```

## Arguments

- params_to_estimate:

  A character vector of parameter names to estimate. Only scalar-valued
  model parameters will be considered.

- proposal_variance:

  A named list of proposal variances for each parameter. If `NULL`, a
  default set of variances is used.

- correlation_matrix:

  Optional correlation matrix for a subset of parameters. If provided,
  the correlation structure will be incorporated into the proposal
  matrix.

- model:

  A model object with a [`new()`](https://rdrr.io/r/methods/new.html)
  constructor and a `param()` method that returns all parameters in the
  model.

- param_inputs:

  A named list of parameter values used to initialize the model. This is
  used to identify scalar parameters that are eligible for estimation.

## Value

A square numeric matrix with dimensions equal to the number of scalar
parameters in the model. The diagonal contains the proposal variances,
and the off-diagonal entries contain covariances if a correlation matrix
is provided.
