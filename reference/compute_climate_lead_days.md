# Compute how many days Climate must lead SMC

Given SMC and Climate data.frames (both with a Date column `dates`),
returns the integer number of days by which the Climate series must
start earlier than the SMC series start.

## Usage

``` r
compute_climate_lead_days(smc_df, clim_df)
```

## Arguments

- smc_df:

  Data.frame with a Date column `dates`.

- clim_df:

  Data.frame with a Date column `dates`.

## Value

Integer \>= 0: days by which Climate must precede SMC.
