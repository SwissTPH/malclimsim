# Simulate Data for Inference

This function generates a simulated dataset for inference based on a
specified model and parameter values. It supports options for monthly
aggregation, adjusting for unequal days in months, and adding stochastic
noise to the incidence data.

## Usage

``` r
data_sim_for_inference(
  model,
  param_inputs,
  dates,
  noise = FALSE,
  month = FALSE,
  month_unequal_days = FALSE
)
```

## Arguments

- model:

  A compiled model object compatible with
  [`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md).
  Typically created using
  [`odin.dust::odin_dust()`](https://rdrr.io/pkg/odin.dust/man/odin_dust.html).

- param_inputs:

  A named list of parameter values required by the model.

- dates:

  A character vector of length two specifying the start and end dates
  (e.g., `c("2022-01-01", "2022-12-31")`).

- noise:

  Logical. If `TRUE`, random noise is added using a negative binomial
  distribution. Defaults to `FALSE`.

- month:

  Logical. If `TRUE`, the output is aggregated by month. Defaults to
  `FALSE`.

- month_unequal_days:

  Logical. If `TRUE`, accounts for unequal month lengths when
  aggregating. Defaults to `FALSE`.

## Value

A data frame of simulated incidence data, in the same format as the
output of
[`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md).
If `noise = TRUE`, the `inc` column will include added random variation.

## Details

- When `month = TRUE`, output is aggregated by month. If
  `month_unequal_days = TRUE`, aggregation is adjusted to reflect actual
  month lengths.

- If `noise = TRUE`, the function adds variability to the incidence
  estimates using a negative binomial distribution with size = 100.

## See also

[`data_sim`](https://swisstph.github.io/malclimsim/reference/data_sim.md)
for the underlying simulation function.

## Examples

``` r
# Example (requires a model created using odin.dust, not run here):
if (FALSE) { # \dontrun{
model <- odin.dust::odin_dust("path/to/model.R")
param_inputs <- list(beta = 0.3, gamma = 0.1)
dates <- c("2022-01-01", "2022-12-31")

# Simulate data without noise
sim_data <- data_sim_for_inference(model, param_inputs, dates, noise = FALSE)

# Simulate data with noise and monthly aggregation
sim_data_monthly <- data_sim_for_inference(model, param_inputs, dates,
                                           noise = TRUE, month = TRUE)
} # }
```
