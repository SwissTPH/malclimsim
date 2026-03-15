# Combine Posterior Simulations, Deterministic Prediction, and Observed Data for PPC Plotting

Combine Posterior Simulations, Deterministic Prediction, and Observed
Data for PPC Plotting

## Usage

``` r
prepare_ppc_data(ci_data, best_df, obs_df, sim_var, obs_var)
```

## Arguments

- ci_data:

  Data frame with posterior median and credible intervals. Must include
  `date_ymd`, `variable`, `median`, `lower`, and `upper`.

- best_df:

  Long-format data frame from deterministic model run, including columns
  `date_ymd`, `variable`, and `value`.

- obs_df:

  Wide-format observed data with columns `date_ymd` and `obs_var`.

- sim_var:

  Character. Name of the variable used for simulated predictions (e.g.,
  "inc_C_transformed").

- obs_var:

  Character. Name of the observed variable in `obs_df` (e.g., "inc_C").

## Value

A tibble combining observed, deterministic, and simulated prediction
data, ready for
[`plot_ppc()`](https://swisstph.github.io/malclimsim/reference/plot_ppc.md).
