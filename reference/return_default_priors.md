# Return Default Priors Used for Inference Procedure

This function returns a list of default priors used in the inference
procedure. Each parameter has an associated prior distribution, an
initial value, and range constraints.

## Usage

``` r
return_default_priors()
```

## Value

A list where each element corresponds to a parameter. Each parameter
element is a list containing:

- initial:

  Initial value for the parameter.

- min:

  Minimum value for the parameter.

- max:

  Maximum value for the parameter.

- prior:

  A function defining the prior distribution of the parameter, returning
  the log probability density.

## Examples

``` r
# Retrieve the default priors
priors <- return_default_priors()

# List the names of all parameters with priors
names(priors)
#>  [1] "lag_R"        "lag_T"        "alpha"        "sigma_LT"     "sigma_RT"    
#>  [6] "R_opt"        "k1"           "size_1"       "eff_SMC"      "s"           
#> [11] "b"            "qR"           "size_2"       "kappa_C"      "kappa_A"     
#> [16] "z"            "beta_1"       "beta_2"       "T_opt"        "c_X"         
#> [21] "fT_C"         "mu_TS"        "mu_IR"        "mu_RS"        "mu_EI"       
#> [26] "delta_b"      "delta_d"      "delta_a"      "N"            "percAdult"   
#> [31] "pi_s_1"       "c_s"          "clim_SMC_lag" "r_C_0"        "r_A_0"       

# Access the prior definition for a specific parameter
if ("phi" %in% names(priors)) {
  phi_prior <- priors$phi
  print(phi_prior)

  # Evaluate the log-prior density at a specific value
  phi_value <- 0.5
  log_density <- phi_prior$prior(phi_value)
  print(log_density)
}
```
