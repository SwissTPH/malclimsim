# Calculate Normalized SMC Decay Array Over Time

Generates a vector of normalized decayed SMC efficacies for a time
series where SMC is applied intermittently. Each decay curve is scaled
so that its maximum value is 1.

## Usage

``` r
calc_decay_arr(SMC, decay_func = decay_SMC, const = -0.1806)
```

## Arguments

- SMC:

  Numeric vector. Time series indicating the presence (positive value)
  or absence (zero) of SMC at each time point.

- decay_func:

  Function. Function defining the decay curve. Default is `decay_SMC`.

- const:

  Numeric. Decay constant to pass to `decay_func`. Default is -0.1806.

## Value

Numeric vector of the same length as `SMC`, with decay values normalized
to a max of 1 per SMC round.

## Examples

``` r
smc_series <- c(0, 1, 0, 0, 0, 1, 0, 0)
calc_decay_arr(smc_series)
#> [1] 0.0000000 0.8976710 0.8972115 0.8966616 0.8960037 0.8976710 0.8972115
#> [8] 0.8966616
```
