# Observed Malaria Incidence Used for Model Calibration

Weekly simulated malaria incidence data used to represent observed cases
for calibration and validation of the transmission model. Generated
deterministically using maximum a posteriori (MAP) parameter values.

## Usage

``` r
obs_cases
```

## Format

A data frame with the following columns:

- date_ymd:

  Date corresponding to each week (class: Date)

- inc_A:

  Simulated incidence in adults

- inc_C:

  Simulated incidence in children

- week_no:

  Week index relative to simulation start

## Source

Generated using
[`data_sim()`](https://swisstph.github.io/malclimsim/reference/data_sim.md)
with MAP-calibrated model parameters.
