avg_across_loc <- function(dat, R = TRUE){
  if(R){
    avg_dat <- dat %>% group_by(date) %>% summarize(rainfall = mean(chirps))  # Average daily rainfall
  } else {
    avg_dat <- dat %>% group_by(date) %>% summarize(temp = mean(chirts))  # Average daily temperature
  }
  return(avg_dat)
}

date_bounds_from_years <- function(years){
  # Determine the first and last year from the vector
  start_year <- min(years)  # This will give 2014
  end_year <- max(years)    # This will give 2023

  # Create the first and last dates
  first_date <- as.Date(paste0(start_year, "-01-01"))  # January 1, 2014
  last_date <- as.Date(paste0(end_year, "-12-31"))     # December 31, 2023

  return(c(first_date, last_date))
}

save_CHIRPS <- function(lon, lat, years, save = TRUE, path_to_data = NULL) {
  # Create a data frame for the single coordinate
  coordinates <- data.frame(longitude = lon, latitude = lat)

  # Convert the coordinates to a simple feature object
  tp_point <- st_as_sf(coordinates, coords = c("longitude", "latitude"), crs = 4326)  # WGS84 CRS

  dates_for_rain <- as.character(date_bounds_from_years(years))

  # Retrieve CHIRPS rainfall data for the sampled points and date range
  dat <- chirps::get_chirps(tp_point, dates = dates_for_rain, server = "ClimateSERV")

  # Calculate the average rainfall across the sampled locations
  avg_rain <- avg_across_loc(dat)

  # Optionally save the average rainfall data
  if (save) {
    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%m%d%H%M")
    file_name <- paste0("chirps_", formatted_datetime, ".rds")
    saveRDS(avg_rain, paste0(path_to_data, file_name))
  }

  return(avg_rain)  # Return the average rainfall data
}

rolling_average_D_days <- function(avg_rain, D, save = TRUE, file = "", rain = TRUE) {

  # Extract date and rainfall columns
  c1 <- as.Date(as.matrix(avg_rain[1]))  # Date column
  c2 <- as.numeric(as.matrix(avg_rain[2]))  # Rainfall column

  # Create a time-series object using zoo
  x <- zoo::zoo(c2, c1)

  # Calculate rolling average with a window of D days
  avg_rain$rollrain <- zoo::rollmean(x, D, fill = NA, align = "right")

  # Optionally save the updated data with rolling average
  if (save) {
    saveRDS(avg_rain, file = paste(dir, file, sep = ""))
  }

  return(avg_rain)  # Return the updated data frame
}

rolling_average_temp_D_days <- function(avg_temp, D, save = TRUE, file = "") {

  # Extract date and temperature columns
  c1 <- as.Date(as.matrix(avg_temp[1]))  # Date column
  c2 <- as.numeric(as.matrix(avg_temp[2]))  # Temperature column

  # Create a time-series object using zoo
  x <- zoo::zoo(c2, c1)

  # Calculate rolling average with a window of D days
  avg_temp$rolltemp <- zoo::rollmean(x, D, fill = NA, align = "right")

  # Optionally save the updated data with rolling average
  if (save) {
    saveRDS(avg_temp, file = file)
  }

  return(avg_temp)  # Return the updated data frame
}



standardize_rainfall <- function(avg_rain, save = TRUE, file = "") {

  # Calculate z-scores for rainfall anomalies
  avg_rain$anom <- (avg_rain$rollrain - mean(avg_rain$rollrain, na.rm = TRUE)) /
    sd(avg_rain$rollrain, na.rm = TRUE)

  # Extract date components (month, week, day) for further analysis
  avg_rain$month <- month(avg_rain$date)
  avg_rain$week <- week(avg_rain$date)
  avg_rain$day <- day(avg_rain$date)

  # Optionally save the updated data
  if (save) {
    saveRDS(avg_rain, file = paste(dir, file, sep = ""))
  }

  return(avg_rain)  # Return the updated data frame
}

#' Save ERA5 Reanalysis Data
#'
#' Downloads ERA5 reanalysis data for specified years and saves it to the specified path.
#'
#' @param years A vector of years for which the data is to be downloaded (e.g., \code{c(2000, 2001)}).
#' @param path A string specifying the path where the downloaded data file should be stored.
#'
#' @return The path to the downloaded NetCDF file.
#' @export
#'
#' @examples
#' \dontrun{
#' save_era5(years = c(2000, 2001), path = "./data")
#' }
save_era5 <- function(years, path){
  current_datetime <- Sys.time()
  formatted_datetime <- format(current_datetime, "%m%d%H%M")
  target <- paste0("era5_", formatted_datetime, ".nc")

  request <- list(
    dataset_short_name = "reanalysis-era5-land-monthly-means",
    product_type = "monthly_averaged_reanalysis",
    variable = "2m_temperature",
    year = as.character(years),
    month = c("01", "02", "03",
              "04", "05", "06",
              "07", "08", "09",
              "10", "11", "12"),
    time = c("00:00"),
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(37, -20, -35, 52),  # Updated bounds: (North, West, South, East)
    target = target
  )

  file <- wf_request(
    request  = request,  # the request
    transfer = TRUE,     # download the file
    path = path       # store data in current working directory
  )
}

#' Extract Temperature Data from ERA5 File
#'
#' Extracts temperature data for a specified longitude and latitude from a NetCDF ERA5 file.
#'
#' @param lon A numeric value specifying the longitude of the location.
#' @param lat A numeric value specifying the latitude of the location.
#' @param path_to_file A string specifying the path to the NetCDF file containing ERA5 data.
#'
#' @return A dataframe containing two columns: \code{Date} (POSIXct) and \code{Temperature} (in Celsius).
#' @export
#'
#' @examples
#' \dontrun{
#' extract_era5(lon = 10.5, lat = -25.3, path_to_file = "./data/era5_011220.nc")
#' }
extract_era5 <- function(lon, lat, path_to_file){
  nc <- ncdf4::nc_open(path_to_file)
  temp_values <- ncdf4::ncvar_get(nc, "t2m")
  lat_values <- ncdf4::ncvar_get(nc, "latitude")
  lon_values <- ncdf4::ncvar_get(nc, "longitude")
  valid_time <- ncdf4::ncvar_get(nc, "valid_time")  # Extract valid_time

  dim1_ind <- which.min(abs(lon_values - lon))  # Closest longitude index
  dim2_ind <- which.min(abs(lat_values - lat))  # Closest latitude index

  temp <- temp_values[dim1_ind, dim2_ind,]

  # Convert 'valid_time' to readable dates
  dates <- as.POSIXct(valid_time, origin = "1970-01-01", tz = "UTC")

  # Create a dataframe with dates and temperature
  temp_df <- data.frame(
    Date = dates,
    Temperature = temp - 273.15
  )
  return(temp_df)
}


#' Daily Smoothed Temperature Function
#'
#' This function takes a data frame with monthly temperature data and produces a daily smoothed
#' temperature series using a smoothing spline. The resulting data frame includes predicted temperatures
#' for every day in the range from the earliest to the latest date in the original data.
#'
#' @param temp_df A data frame containing two columns:
#'   \describe{
#'     \item{Date}{A character or numeric column representing the year and month (e.g., "2024-11").}
#'     \item{Temperature}{A numeric column representing the temperature for the corresponding month.}
#'   }
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{Date}{A Date column representing daily dates.}
#'     \item{Temperature}{A numeric column with the smoothed daily temperature predictions.}
#'   }
#'
#' @examples
#' # Example usage:
#' temp_data <- data.frame(
#'   Date = c("2024-01", "2024-02", "2024-03"),
#'   Temperature = c(30.5, 28.0, 25.3)
#' )
#' daily_smooth_temp(temp_data)
#'
#' @importFrom stats smooth.spline predict
#' @export
daily_smooth_temp <- function(temp_df) {
  # Convert the Date column to Date type
  temp_df$Date <- as.Date(paste0(temp_df$Date, "-01"))

  # Convert dates to numeric format (number of days since 1970-01-01)
  temp_df$Date_numeric <- as.numeric(temp_df$Date)

  # Fit a smoothing spline using the numeric date values
  spline_fit <- smooth.spline(x = temp_df$Date_numeric, y = temp_df$Temperature)

  # Find the last date in temp_df
  last_date <- max(temp_df$Date)

  # Determine the end of the month for the last_date
  end_of_month <- as.Date(format(last_date, "%Y-%m-01")) + months(1) - 1

  # Create a numeric sequence from the first date to the end of the last month
  predicted_dates_numeric <- seq(min(temp_df$Date_numeric), as.numeric(end_of_month), by = 1)

  # Predict temperatures using the smoothing spline
  predicted_temps <- predict(spline_fit, x = predicted_dates_numeric)

  # Convert numeric predictions back to Date format
  predicted_dates <- as.Date(predicted_dates_numeric, origin = "1970-01-01")

  # Combine the predicted data into a data frame
  temp_df <- data.frame(Date = predicted_dates, Temperature = predicted_temps$y)

  return(temp_df)
}



#' Saving Climate Data From ERA5 and CHIRTSdaily
#'
#' @param lon - longitude coordinate
#' @param lat - latitude coordinate
#' @param years - vector containing each year for data to be downloaded
#' @param path_to_data - folder where data is to be saved
#' @param rain - boolean for if CHIRPS data should be downloaded and saved
#' @param temp - boolean for if ERA5 data should be downloaded and saved
#'
#' @return - nothing
#' @export
#'
#' @examples
#' years <- 2014:2023
#' lon <- 17.9
#' lat <- 8.3
#' path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"
#' save_climate_data(lon, lat, years, path_to_data)
save_climate_data <- function(lon, lat, years, path_to_data, rain = TRUE, temp = TRUE){
  if(temp){
    save_era5(years, path = path_to_data)
  }
  if(rain){
    save_CHIRPS(lon, lat, years, save = TRUE, path_to_data = path_to_data)
  }
}

#' Title
#'
#' @param lon - longitude coordinate
#' @param lat - latitude coordinate
#' @param years - vector containing each year for the analysis
#' @param D1 - the number of previous days to take an average over for the rainfall rolling mean
#' @param D2 - the number of previous days to take an average over for the temperature rolling mean
#' @param temp_path - path to temperature data saved by save_climate_data
#' @param rain_path path to rainfall data saved by save_climate_data
#' @param path_to_data path to directory where the processed climate data will be stored
#'
#' @return a data frame
#' @export
#'
#' @examples
#' years <- 2014:2023
#' lon <- 17.9
#' lat <- 8.3
#' path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"
#' save_climate_data(lon, lat, years, path_to_data)
#' temp_path <- paste0(path_to_data, "era5_11051417.nc")
#' rain_path <- paste0(path_to_data, "chirps_11051418.rds")
#' met <- process_climate_data(lon, lat, years, D1 = 30, D2 = 30, temp_path = temp_path,
#' rain_path = rain_path, path_to_data
process_climate_data <- function(lon, lat, years, D1, D2, temp_path, rain_path, path_to_data, months_30_days = TRUE){
  temp_df <- extract_era5(lat = lat, lon = lon, path_to_file = temp_path)
  temp <- daily_smooth_temp(temp_df)

  # Rainfall
  avg_rain <- readRDS(rain_path)

  temp <- temp[which(temp$Date %in% as.Date(intersect(as.character(temp$Date), as.character(avg_rain$date)))),]

  # - `rolling_average_D_days` calculates a rolling average of rainfall with a window size of `D1` days (30 days here)
  # - `rolling_avg_rwa$rollmean` contains the calculated rolling averages
  rolling_rain <- rolling_average_D_days(avg_rain, D1, save = FALSE)
  rolling_temp <- rolling_average_temp_D_days(temp, D2, save = FALSE)
  # Fill missing values in the rolling averages by extending the last available value forward
  # - `na.fill(..., c("extend"))` fills any NA values by extending the last non-NA value forward
  rolling_rain$rollrain <- zoo::na.fill(rolling_rain$rollrain, c("extend"))
  rolling_temp$rolltemp <- zoo::na.fill(rolling_temp$rolltemp, c("extend"))

  # Standardize the rainfall data using z-scores (anomalies from the mean)
  # - `standardize_rainfall` calculates the anomaly (z-score) for rainfall to identify deviations from the norm
  met <- standardize_rainfall(rolling_rain, save = FALSE)
  # Combine rainfall and temperature
  # Overwrite the `rwa_met` object with average temperature data from `temp_rwa`
  # - `temp_rwa$tavg` contains the average daily temperature extracted earlier
  met$temp <- temp$Temperature
  met$rolltemp <- rolling_temp$rolltemp
  if(months_30_days){
    met <- climate_to_30_day_months(met, start_year = years[1], end_year = years[length(years)])
  }
  # Save the processed Rwanda climate data (rainfall and temperature) to the specified directory
  # - The data is saved as an RDS file, which is an efficient format for storing R objects
  current_datetime <- Sys.time()
  formatted_datetime <- format(current_datetime, "%m%d%H%M")
  saveRDS(met, paste0(path_to_data, "met_", D1, "_", formatted_datetime, ".rds"))

  return(met)
}
