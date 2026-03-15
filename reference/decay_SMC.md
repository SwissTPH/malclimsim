# Compute SMC Efficacy Decay Over Time

Models the decay of Seasonal Malaria Chemoprevention (SMC) efficacy as a
sigmoid function of time since last administration.

## Usage

``` r
decay_SMC(x, const = -0.1806)
```

## Arguments

- x:

  Numeric. Number of days since the last SMC treatment.

- const:

  Numeric. The decay constant that controls the steepness of the sigmoid
  curve. Default is -0.1806.

## Value

Numeric. Remaining efficacy of SMC treatment at day `x`.

## Examples

``` r
decay_SMC(0)      # Should return near maximum efficacy
#> [1] 0.897671
decay_SMC(60)     # Should return reduced efficacy
#> [1] 0.006774873
```
