# Run the model once at the maximum-posterior parameter set (deterministic)

Run the model once at the maximum-posterior parameter set
(deterministic)

## Usage

``` r
simulate_best(model, param_inputs, best_params, start_date, end_date, ...)
```

## Arguments

- model:

  A dust or odin model object.

- param_inputs:

  Named list of baseline inputs.

- best_params:

  Named vector or list of best-fit parameters.

- start_date:

  Date or string start of sim window.

- end_date:

  Date or string end of sim window.

- ...:

  Additional args passed to
  [`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md).

## Value

A data.frame in long form with `date_ymd`, `variable`, `value`.
