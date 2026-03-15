# Simulate Incidence Data Over Time (Monthly or Weekly)

Runs a deterministic dust model simulation with an optional pre‑warm
period, then aggregates and returns incidence (and optional transformed
incidence) at either monthly or weekly resolution. Weekly output may use
a 360‑day or 365‑day year for the prewarm and simulation period.

## Usage

``` r
data_sim(
  model,
  param_inputs,
  start_date,
  end_date,
  prewarm_years = 2,
  month = FALSE,
  round = TRUE,
  save = TRUE,
  file = "",
  month_unequal_days = FALSE,
  return_EIR = FALSE,
  return_compartments = FALSE,
  mu_transform_A = NULL,
  mu_transform_C = NULL,
  covariate_matrix = NULL,
  noise = FALSE,
  size = NULL
)
```

## Arguments

- model:

  A dust model object (as returned by `load_model`) to simulate.

- param_inputs:

  Named list of model parameters, including any time‑varying vectors
  (e.g. `temp`, `c_R_D`, `cov_SMC`).

- start_date:

  Character or `Date` giving the first day of the analysis.

- end_date:

  Character or `Date` giving the last day of the analysis.

- prewarm_years:

  Integer number of years (360‑ or 365‑day) to run before `start_date`
  as warm‑up.

- month:

  Logical; if `TRUE`, returns monthly aggregates (30‑day or unequal
  month lengths if `month_unequal_days = TRUE`); if `FALSE`, returns
  weekly aggregates.

- round:

  Logical; if `TRUE`, round all incidence counts to integers.

- save:

  Logical; if `TRUE`, save the returned data frame to disk as an RDS
  file.

- file:

  Character; path/filename to use when saving (if `save = TRUE`).

- month_unequal_days:

  Logical; when `month = TRUE`, if `TRUE` use actual month boundaries
  from `param_inputs\$day_count`, otherwise use fixed 30‑day months.

- return_EIR:

  Logical; if `TRUE`, include monthly EIR values.

- return_compartments:

  Logical; currently unused.

- mu_transform_A:

  Optional function to transform adult incidence after simulation.
  Should accept `(inc_df, param_inputs)` and return a vector.

- mu_transform_C:

  Optional function to transform child incidence after simulation.
  Should accept `(inc_df, param_inputs)` and return a vector.

- covariate_matrix:

  Optional data frame of external covariates (must have a date column in
  the same name/location as `date_ymd` in the output); will be joined by
  `date_ymd`.

- noise:

  Logical; if `TRUE`, add negative‑binomial noise to each weekly or
  monthly incidence draw, using `size` as dispersion.

- size:

  Optional numeric; dispersion parameter for negative‑binomial noise.

## Value

A data frame with columns:

- date_ymd:

  Date of each week (Sunday) or month (start day).

- week_no, month_no:

  Index of each period, starting at 0.

- inc_A, inc_C, inc:

  Adult, child, and total incidence.

- EIR_monthly:

  Monthly EIR, if `return_EIR = TRUE`.

- inc_A_transformed, inc_C_transformed:

  Transformed incidence, if `mu_transform_A` or `mu_transform_C`
  supplied.

plus any merged covariates from `covariate_matrix`.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_df <- data_sim(
  model             = model_1,
  param_inputs      = param_inputs_wk,
  start_date        = "2018-01-01",
  end_date          = "2023-12-31",
  prewarm_years     = 5,
  month             = FALSE,
  round             = FALSE,
  save              = FALSE,
  mu_transform_C    = mu_transform_C_wk,
  covariate_matrix  = covariate_matrix_wk,
  noise             = TRUE,
  size              = 10,
  year_360_days     = FALSE
)
} # }
```
