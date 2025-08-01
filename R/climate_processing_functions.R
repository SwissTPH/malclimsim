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

save_CHIRPS <- function(lon, lat, years, save = TRUE, path_to_data = NULL, file_name = NULL) {
  # Create a data frame for the single coordinate
  coordinates <- data.frame(longitude = lon, latitude = lat)

  # Convert the coordinates to a simple feature object
  tp_point <- st_as_sf(coordinates, coords = c("longitude", "latitude"), crs = 4326)  # WGS84 CRS

  dates_for_rain <- as.character(date_bounds_from_years(years))

  # Retrieve CHIRPS rainfall data for the sampled points and date range
  dat <- chirps::get_chirps(tp_point, dates = dates_for_rain, server = "ClimateSERV")

  # Calculate the average rainfall across the sampled locations
  avg_rain <- avg_across_loc(dat)

  if(is.null(file_name)){
    file_name <- paste0("chirps_", formatted_datetime, ".rds")
  }
  # Optionally save the average rainfall data
  if (save) {
    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%m%d%H%M")
    file_name <- file_name
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


standardize_rainfall <- function(cum_rain, save = TRUE, file = "") {

  # Calculate z-scores for rainfall anomalies
  cum_rain$anom <- (cum_rain$CumulativeRainfall - mean(cum_rain$CumulativeRainfall, na.rm = TRUE)) /
    sd(cum_rain$CumulativeRainfall, na.rm = TRUE)

  # Extract date components (month, week, day) for further analysis
  cum_rain$month <- month(cum_rain$date)
  cum_rain$week <- week(cum_rain$date)
  cum_rain$day <- day(cum_rain$date)

  # Optionally save the updated data
  if (save) {
    saveRDS(cum_rain, file = paste(dir, file, sep = ""))
  }

  return(cum_rain)  # Return the updated data frame
}


#' Save ERA5 Reanalysis Data
#'
#' Downloads ERA5 reanalysis data for specified years and saves it to the specified path.
#'
#' @param years A vector of years for which the data is to be downloaded (e.g., \code{c(2000, 2001)}).
#' @param path A string specifying the path where the downloaded data file should be stored.
#' @param file_name Optional string of file name. If none is specified, a default time-based name is given.
#'
#' @return The path to the downloaded NetCDF file.
#' @export
#'
#' @examples
#' \dontrun{
#' save_era5(years = c(2000, 2001), path = "./data")
#' }
save_era5 <- function(years, path, file_name = NULL){
  if(is.null(file_name)){
    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%m%d%H%M")
    file_name <- paste0("era5_", formatted_datetime, ".nc")
  }

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
    target = file_name
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
#' temperature series using a smoothing spline. The output includes predicted temperatures
#' for each day from the earliest to the latest date in the original data.
#'
#' @param temp_df A data frame with two columns:
#' \itemize{
#'   \item \code{Date}: A character vector in "YYYY-MM" format (e.g., "2024-11") or a Date object representing the first day of the month.
#'   \item \code{Temperature}: A numeric vector giving the temperature for each month.
#' }
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{Date}: A Date vector representing daily dates.
#'   \item \code{Temperature}: A numeric vector with the smoothed daily temperature predictions.
#' }
#'
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

calculate_cumulative_rainfall <- function(rain_path, time_period = c("monthly", "weekly")) {
  # Read the RDS file
  avg_rain <- readRDS(rain_path)

  # Ensure the date column is in Date format
  avg_rain$date <- as.Date(avg_rain$date)

  # Determine the aggregation period and calculate cumulative rainfall
  if (time_period == "monthly") {
    cum_rain <- avg_rain %>%
      mutate(period = floor_date(date, "month")) %>%
      group_by(period) %>%
      summarise(cumulative_rainfall = sum(rainfall, na.rm = TRUE)) %>%
      ungroup()
  } else if (time_period == "weekly") {
    cum_rain <- avg_rain %>%
      mutate(period = floor_date(date, "week", week_start = 1)) %>%  # Week starts on Monday
      group_by(period) %>%
      summarise(cumulative_rainfall = sum(rainfall, na.rm = TRUE)) %>%
      ungroup()
  } else {
    stop("Invalid time_period. Choose 'monthly' or 'weekly'.")
  }
  return(cum_rain)
}

#' Smooth Daily Rainfall Data
#'
#' This function uses a smoothing spline to interpolate daily cumulative rainfall data from a dataset containing
#' weekly or monthly cumulative rainfall values.
#'
#' @param cum_rain A data frame with two columns:
#' \itemize{
#'   \item \code{period}: A column representing the time period (e.g., weekly dates), of class \code{Date}.
#'   \item \code{cumulative_rainfall}: A column representing the cumulative rainfall for each period.
#' }
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{Date}: A column with daily dates, interpolated from the original period.
#'   \item \code{CumulativeRainfall}: A column with smoothed daily cumulative rainfall values.
#' }
#'
#' @examples
#' # Example dataset with weekly cumulative rainfall
#' weekly_dates <- seq(as.Date("2023-01-01"), as.Date("2023-12-24"), by = "7 days")
#' cum_rain <- data.frame(
#'   period = weekly_dates,
#'   cumulative_rainfall = cumsum(runif(length(weekly_dates), min = 0, max = 20))
#' )
#'
#' # Smooth to daily cumulative rainfall
#' daily_rain <- daily_smooth_rain(cum_rain)
#'
#' # Plot the results
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(daily_rain, ggplot2::aes(x = Date, y = CumulativeRainfall)) +
#'     ggplot2::geom_line(color = "blue") +
#'     ggplot2::labs(title = "Daily Smoothed Cumulative Rainfall",
#'                   x = "Date", y = "Cumulative Rainfall") +
#'     ggplot2::theme_minimal()
#' }
#'
#' @export
daily_smooth_rain <- function(cum_rain) {
  # Convert the period to numeric format
  cum_rain$Date_numeric <- as.numeric(cum_rain$period)

  # Fit a smoothing spline using the numeric date values
  spline_fit <- smooth.spline(x = cum_rain$Date_numeric, y = cum_rain$cumulative_rainfall)

  # Create a sequence of daily dates
  predicted_dates_numeric <- seq(min(cum_rain$Date_numeric), max(cum_rain$Date_numeric), by = 1)

  # Predict daily cumulative rainfall using the smoothing spline
  predicted_rain <- predict(spline_fit, x = predicted_dates_numeric)

  # Convert numeric predictions back to Date format
  predicted_dates <- as.Date(predicted_dates_numeric, origin = "1970-01-01")

  # Combine the predicted data into a data frame
  daily_rain_df <- data.frame(Date = predicted_dates, CumulativeRainfall = predicted_rain$y)

  return(daily_rain_df)
}

#' Saving Climate Data From ERA5 and CHIRTSdaily
#'
#' @param lon - longitude coordinate
#' @param lat - latitude coordinate
#' @param years - vector containing each year for data to be downloaded
#' @param path_to_data - folder where data is to be saved
#' @param rain - boolean for if CHIRPS data should be downloaded and saved
#' @param temp - boolean for if ERA5 data should be downloaded and saved
#' @param temp_file_name - character file name for stored temperature data (ERA5)
#' @param rain_file_name - character file name for stored rainfall data (CHIRPS)
#'
#' @return - nothing
#' @export
#'
#' @examples
#' \dontrun{
#' years <- 2014:2023
#' lon <- 17.9
#' lat <- 8.3
#' path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"
#' save_climate_data(lon, lat, years, path_to_data)
#' }
save_climate_data <- function(lon, lat, years, path_to_data, rain = TRUE, temp = TRUE,
                              temp_file_name = NULL, rain_file_name = NULL){
  if(temp){
    save_era5(years, path = path_to_data, file_name = temp_file_name)
  }
  if(rain){
    save_CHIRPS(lon, lat, years, save = TRUE, path_to_data = path_to_data, file_name = rain_file_name)
  }
}


#' Process and Format Climate Data for Modeling
#'
#' Reads and processes daily temperature and rainfall data, applies smoothing, standardizes rainfall,
#' aligns dates, and optionally transforms the data to use 30-day months. Returns a clean climate dataset
#' ready for use in modeling workflows.
#'
#' @param lon Numeric. Longitude coordinate for the location of interest.
#' @param lat Numeric. Latitude coordinate for the location of interest.
#' @param years Integer vector. Years to include in the climate dataset.
#' @param temp_path Character. Path to the saved ERA5 temperature data file (NetCDF format).
#' @param rain_path Character. Path to the saved CHIRPS rainfall data file (RDS format).
#' @param months_30_days Logical. If TRUE (default), transforms calendar months into uniform 30-day months.
#' @param save Logical. If TRUE, saves the processed data as an RDS file (default is FALSE).
#' @param path_to_data Character. Directory where the processed climate data will be saved (required if `save = TRUE`).
#'
#' @return A data frame with daily or monthly climate data including standardized rainfall and smoothed temperature.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' years <- 2014:2023
#' lon <- 17.9
#' lat <- 8.3
#' path_to_data <- "C:/Users/putnni/Documents/r-packages/data/"
#' save_climate_data(lon, lat, years, path_to_data)
#' temp_path <- paste0(path_to_data, "era5_11051417.nc")
#' rain_path <- paste0(path_to_data, "chirps_11051418.rds")
#' met <- process_climate_data(lon, lat, years, D1 = 30, D2 = 30, temp_path = temp_path,
#' rain_path = rain_path, path_to_data)
#' }
process_climate_data <- function(lon, lat, years, temp_path, rain_path, months_30_days = TRUE, save = FALSE, path_to_data = NULL){
  temp_df <- extract_era5(lat = lat, lon = lon, path_to_file = temp_path)
  temp <- daily_smooth_temp(temp_df)

  # Rainfall
  avg_rain <- readRDS(rain_path)
  #cum_rain <- calculate_cumulative_rainfall(rain_path, time_period = "monthly") %>% daily_smooth_rain()
  cum_rain <- calculate_cumulative_rainfall(rain_path, time_period = "monthly") %>% daily_smooth_rain() %>% filter_by_year("Date", years[1]:years[length(years)])
  colnames(cum_rain)[[1]] <- "date"
  temp <- temp[which(temp$Date %in% as.Date(intersect(as.character(temp$Date), as.character(cum_rain$date)))),]
  met <- standardize_rainfall(cum_rain, save = FALSE)
  met$temp <- temp$Temperature

  avg_rain <- avg_rain %>% filter(date %in% as.Date(met$date[1] : met$date[nrow(met)]))
  met$rainfall <- avg_rain$rainfall


  if(months_30_days){
    met <- climate_to_30_day_months(met, start_year = years[1], end_year = years[length(years)])
  }


  colnames(met)[1] <- "dates"
  # Save the processed Rwanda climate data (rainfall and temperature) to the specified directory
  # - The data is saved as an RDS file, which is an efficient format for storing R objects
  current_datetime <- Sys.time()
  formatted_datetime <- format(current_datetime, "%m%d%H%M")
  if(save){saveRDS(met, paste0(path_to_data, "met_", "_", formatted_datetime, ".rds"))}

  return(met)
}
