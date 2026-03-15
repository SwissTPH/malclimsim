# Build a list of mcstate::pmcmc_parameter() objects, starting from the return_default_priors() template and then layering on any overrides.

Build a list of mcstate::pmcmc_parameter() objects, starting from the
return_default_priors() template and then layering on any overrides.

## Usage

``` r
build_priors(
  param_inputs,
  proposal_matrix,
  params_to_estimate,
  override_priors = NULL
)
```

## Arguments

- param_inputs:

  Named list of initial values (so that we know which parameters
  actually exist).

- proposal_matrix:

  Numeric matrix (rownames must match parameter names).

- params_to_estimate:

  Character vector of names we actually want to estimate.

- override_priors:

  Optional named list of lists: each element must match the structure
  returned by return_default_priors().

## Value

Named list of
[`mcstate::pmcmc_parameter`](https://rdrr.io/pkg/mcstate/man/pmcmc_parameter.html)
objects for all parameters present in both `param_inputs` and
`proposal_matrix`.
