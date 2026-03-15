# Update a Parameter List with New Values

Replaces entries in a parameter list using a named vector of new values.
Typically used to inject sampled or optimized parameter values (e.g.,
from an MCMC posterior) into a base parameter list.

## Usage

``` r
update_param_list(param_inputs, param_values)
```

## Arguments

- param_inputs:

  A named list containing parameter values, where each name corresponds
  to a model parameter.

- param_values:

  A named vector or list of parameter values to use for updating. Each
  name in this object should correspond to a name in `param_inputs`.

## Value

A modified version of `param_inputs` where parameters matching names in
`param_values` are replaced by the corresponding values.

## Examples

``` r
# Original parameter list
param_inputs <- list(alpha = 1, beta = 2, gamma = 3)

# New parameter values, e.g., from posterior sample
params_at_max_posterior <- c(alpha = 0.9, gamma = 2.7)

# Update the parameter list
updated_params <- update_param_list(param_inputs, params_at_max_posterior)
print(updated_params)
#> $alpha
#> [1] 0.9
#> 
#> $beta
#> [1] 2
#> 
#> $gamma
#> [1] 2.7
#> 
```
