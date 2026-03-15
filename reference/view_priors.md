# View Details of Specified Priors

Returns a data frame summarizing the details of the specified priors,
including initial values, bounds, and the functional form of the prior
distribution.

## Usage

``` r
view_priors(param_inputs, proposal_matrix, params_to_estimate, priors = NULL)
```

## Arguments

- param_inputs:

  A named list of base parameter values. Used if `priors` is not
  supplied.

- proposal_matrix:

  Covariance matrix used in proposal distribution. Required if `priors`
  is NULL.

- params_to_estimate:

  A character vector specifying which parameters' priors should be
  viewed.

- priors:

  Optional list of priors. If `NULL`, the function calls
  `initialize_priors()` to generate defaults.

## Value

A data frame with one row per parameter and the following columns:

- Name:

  The parameter name.

- Initial:

  The initial value used in inference.

- Min:

  The lower bound of the prior.

- Max:

  The upper bound of the prior.

- Description:

  A textual representation of the prior function.

## Examples

``` r
# Define dummy prior functions
dummy_prior <- function(p) dnorm(p, mean = 0, sd = 1, log = TRUE)

# Create a minimal list of priors matching the expected structure
dummy_priors <- list(
  a_R = list(name = "a_R", initial = 0.1, min = -2, max = 2, prior = dummy_prior),
  b_R = list(name = "b_R", initial = 1.5, min = 0, max = 3, prior = dummy_prior),
  qR  = list(name = "qR",  initial = 0.01, min = 0, max = 1, prior = dummy_prior),
  z   = list(name = "z",   initial = 0.2, min = 0, max = 1, prior = dummy_prior),
  eff_SMC = list(name = "eff_SMC", initial = 0.5, min = 0, max = 1, prior = dummy_prior),
  phi = list(name = "phi", initial = 0.3, min = -1, max = 1, prior = dummy_prior),
  size = list(name = "size", initial = 5, min = 0.01, max = 20, prior = dummy_prior)
)

# View priors for a subset of parameters using a supplied list
view_priors(
  param_inputs = NULL,
  proposal_matrix = NULL,
  params_to_estimate = c("a_R", "b_R", "qR", "z", "eff_SMC", "phi", "size"),
  priors = dummy_priors
)
#>            Name Initial   Min Max
#> a_R         a_R    0.10 -2.00   2
#> b_R         b_R    1.50  0.00   3
#> qR           qR    0.01  0.00   1
#> z             z    0.20  0.00   1
#> eff_SMC eff_SMC    0.50  0.00   1
#> phi         phi    0.30 -1.00   1
#> size       size    5.00  0.01  20
#>                                                  Description
#> a_R     function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> b_R     function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> qR      function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> z       function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> eff_SMC function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> phi     function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
#> size    function (p)  dnorm(p, mean = 0, sd = 1, log = TRUE)
```
