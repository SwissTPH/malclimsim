# Saving Climate Data From ERA5 and CHIRTSdaily

Saving Climate Data From ERA5 and CHIRTSdaily

## Usage

``` r
save_climate_data(
  lon,
  lat,
  years,
  path_to_data,
  rain = TRUE,
  temp = TRUE,
  temp_file_name = NULL,
  rain_file_name = NULL
)
```

## Arguments

- lon:

  - longitude coordinate

- lat:

  - latitude coordinate

- years:

  - vector containing each year for data to be downloaded

- path_to_data:

  - folder where data is to be saved

- rain:

  - boolean for if CHIRPS data should be downloaded and saved

- temp:

  - boolean for if ERA5 data should be downloaded and saved

- temp_file_name:

  - character file name for stored temperature data (ERA5)

- rain_file_name:

  - character file name for stored rainfall data (CHIRPS)

## Value

- nothing

## Examples

``` r
if (FALSE) { # \dontrun{
years <- 2014:2023
lon <- 17.9
lat <- 8.3
path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"
save_climate_data(lon, lat, years, path_to_data)
} # }
```
