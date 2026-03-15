# Extract Temperature Data from ERA5 File

Extracts temperature data for a specified longitude and latitude from a
NetCDF ERA5 file.

## Usage

``` r
extract_era5(lon, lat, path_to_file)
```

## Arguments

- lon:

  A numeric value specifying the longitude of the location.

- lat:

  A numeric value specifying the latitude of the location.

- path_to_file:

  A string specifying the path to the NetCDF file containing ERA5 data.

## Value

A dataframe containing two columns: `Date` (POSIXct) and `Temperature`
(in Celsius).

## Examples

``` r
if (FALSE) { # \dontrun{
extract_era5(lon = 10.5, lat = -25.3, path_to_file = "./data/era5_011220.nc")
} # }
```
