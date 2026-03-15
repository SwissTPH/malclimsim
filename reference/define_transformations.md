# Define Transformation Function for MCMC Simulation

Wraps the construction of a transformation function for model covariates
(e.g., temperature, SMC coverage).

## Usage

``` r
define_transformations(temp, c_R_D, SMC, decay, cov_SMC)
```

## Arguments

- temp:

  Numeric vector or matrix of temperature values.

- c_R_D:

  Numeric vector representing case reduction from treatment or
  diagnostics.

- SMC:

  Numeric vector of SMC coverage.

- decay:

  Numeric vector for SMC efficacy decay over time.

- cov_SMC:

  Numeric vector of SMC target coverage.

## Value

A transformation function compatible with `mcstate` model input.
