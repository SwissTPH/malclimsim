# Validate SMC and Climate for Lagged Alignment

Ensures:

1.  SMC `dates` are continuous daily with no gaps.

2.  Climate `dates` start no earlier than (SMC_start - lag_days) and
    cover through SMC_end with no gaps.

3.  For the intended `lag_days`, every SMC_date - lag_days exists in
    climate.

## Usage

``` r
validate_smc_climate_alignment(smc_df, clim_df, lag_days)
```

## Arguments

- smc_df:

  Data.frame with Date column `dates`.

- clim_df:

  Data.frame with Date column `dates`.

- lag_days:

  Integer \>= 0: days Climate lags behind SMC.

## Value

Invisibly TRUE if all checks pass; stops with an error otherwise.
