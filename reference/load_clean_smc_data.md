# Load and clean raw SMC coverage data

Cleans and formats SMC coverage data. Allows for flexibility in input:
the user can supply either a data frame directly or a path to a file
(Excel or RDS). If both `data` and paths are `NULL`, the function
returns an error.

## Usage

``` r
load_clean_smc_data(data = NULL, path_to_excel = NULL, path_to_rds = NULL)
```

## Arguments

- data:

  A data frame containing raw SMC coverage data (default = NULL).

- path_to_excel:

  Optional. Path to an Excel file containing SMC data (default = NULL).

- path_to_rds:

  Optional. Path to an RDS file containing SMC data (default = NULL).

## Value

A cleaned data frame with columns: `date_start` and `coverage`. Only one
row per Year-Month is retained (the first entry chronologically).
