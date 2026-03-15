# Create Index Mapping for Model Variables

This function creates a list of indices for different model variables,
organizing them into `run` and `state` categories based on the provided
info.

## Usage

``` r
index(info)
```

## Arguments

- info:

  A list containing index information, typically provided from a model
  setup. It must include `info$index$...` entries for each variable
  named below.

## Value

A list with two named elements:

- `run`: Named integer vector of indices for running the model (monthly
  and weekly incidence, prevalence, and rates)

- `state`: Integer vector of state indices, in a predefined order
  matching model state variables.

## Examples

``` r
# Create a dummy info object with index entries
info <- list(
  index = as.list(setNames(seq_along(c(
    "month_inc_C", "month_inc_A", "month_inc_total",
    "wk_inc_C",   "wk_inc_A",    "wk_inc_total",
    "prev_C_2",   "prev_A_2",    "prev_C_1",
    "prev_A_1",   "r_C",         "r_A"
  )),
  paste0("index$", c(
    "month_inc_C", "month_inc_A", "month_inc_total",
    "wk_inc_C",   "wk_inc_A",    "wk_inc_total",
    "prev_C_2",   "prev_A_2",    "prev_C_1",
    "prev_A_1",   "r_C",         "r_A"
  ))))
)

# Extract mapping
model_indices <- index(info)
str(model_indices)
#> List of 2
#>  $ run  : NULL
#>  $ state: NULL
```
