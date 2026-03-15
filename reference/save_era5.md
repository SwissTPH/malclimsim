# Save ERA5 Reanalysis Data

Downloads ERA5 reanalysis data for specified years and saves it to the
specified path.

## Usage

``` r
save_era5(years, path, file_name = NULL)
```

## Arguments

- years:

  A vector of years for which the data is to be downloaded (e.g.,
  `c(2000, 2001)`).

- path:

  A string specifying the path where the downloaded data file should be
  stored.

- file_name:

  Optional string of file name. If none is specified, a default
  time-based name is given.

## Value

The path to the downloaded NetCDF file.

## Examples

``` r
if (FALSE) { # \dontrun{
save_era5(years = c(2000, 2001), path = "./data")
} # }
```
