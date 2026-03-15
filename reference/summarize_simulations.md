# Summarize Simulation Outputs with Medians and Confidence Intervals

Computes the pointwise median and confidence intervals for selected
variables across multiple simulation runs. Useful for summarizing
posterior predictive distributions of model outputs.

## Usage

``` r
summarize_simulations(simulation_results, ci_level = 0.95, variables = NULL)
```

## Arguments

- simulation_results:

  A list of data frames, where each data frame contains the model output
  for a different parameter set. Each data frame must include a
  `date_ymd` column.

- ci_level:

  Numeric value between 0 and 1 specifying the confidence level for the
  intervals (e.g., `0.95` for 95% confidence intervals).

- variables:

  A character vector indicating which compartment names to summarize
  (e.g., `c("SC", "EC", "prev_total_with_R")`). If `NULL`, all variables
  except `date_ymd` and `simulation_id` will be included.

## Value

A data frame with one row per date and summary statistics (median,
lower, and upper bounds) for each selected variable. Columns are named
using the pattern `<variable>_median`, `<variable>_lower`, and
`<variable>_upper`.
