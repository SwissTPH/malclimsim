# Calculate Outcome Across Matched Simulation Lists

Calculate Outcome Across Matched Simulation Lists

## Usage

``` r
calculate_estimate(o1, o2, outcome_fn)
```

## Arguments

- o1:

  List of data frames (e.g., with SMC).

- o2:

  List of data frames (e.g., no SMC).

- outcome_fn:

  Function that takes two data frames and returns a scalar outcome.

## Value

Numeric vector of outcomes (e.g., effectiveness or cases averted) for
each parameter set.
