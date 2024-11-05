library(lubridate)
library(ggplot2)
library(ecmwfr)
library(ncdf4)
library(sf)
library(chirps)
library(dplyr)
library(roxygen2)
library(pak)

library(malclimsim)
years <- 2014:2023
lat <- 8.3
lon <- 17.9
path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"

save_climate_data(lon, lat, years, path_to_data)

temp_path <- paste0(path_to_data, "era5_11051417.nc")
rain_path <- paste0(path_to_data, "chirps_11051418.rds")

met <- process_climate_data(lon, lat, years, D = 30, temp_path = temp_path,
                            rain_path = rain_path, path_to_data)

