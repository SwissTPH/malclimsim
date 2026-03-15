# Lag and Trim SMC and Climate Data.frames

After validation, trims & reorders SMC and Climate so that each row *i*
corresponds to the same calendar date, with Climate lagged by
`lag_days`.

## Usage

``` r
lag_and_trim_smc_climate(smc_df, clim_df, lag_days)
```

## Arguments

- smc_df:

  Data.frame with `dates` (Date) and SMC columns.

- clim_df:

  Data.frame with `dates` (Date) and climate columns.

- lag_days:

  Integer \>= 0: days Climate lags behind SMC.

## Value

List with:

- `smc` - ordered SMC data.frame,

- `climate` - climate aligned so row_i = SMC_date_i - lag_days.
