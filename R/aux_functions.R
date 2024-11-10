# Function to convert a date range into a vector of week-year format
# Takes in a start date and an end date (e.g., "2014-01-01") and returns
# weeks in the format "YYYY-WW" (ISO week number). Handles year transitions.
date_to_weeks <- function(start_date, end_date){
  weeks <- c()  # Initialize an empty vector to store weeks

  # Loop through each week in the date range
  while (start_date <= end_date) {
    # Extract year and ISO week number
    year <- format(start_date, "%Y")
    week <- format(start_date, "%V")  # %V gives the ISO-8601 week number

    # Combine year and week to form "YYYY-WW"
    year_week <- paste(year, week, sep = "-")

    # Append this week to the list of weeks
    weeks <- c(weeks, year_week)

    # Move to the next week
    start_date <- start_date + 7
  }

  # Remove any potential duplicate weeks (e.g., at year boundaries) and append the last week
  weeks <- c(unique(weeks), "2022-53")
  return(weeks)
}

# Function to convert a date range into a sequence of months
# It takes start and end dates and returns a list of months in the format "YYYY-MM"
date_to_months <- function(start_date, end_date) {
  # Convert start and end dates from strings to Date objects
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  # Generate a sequence of dates with monthly intervals
  all_months <- seq(start_date, end_date, by = "month")

  # Extract the year and month as "YYYY-MM"
  year_month <- format(all_months, "%Y-%m")

  # Create a dataframe with the monthly dates
  result_df <- data.frame(Date = all_months)

  # Return the vector of dates (first day of each month)
  return(c(result_df$Date))
}

################################################################################
### ------- FUNCTION FOR CONVERTING CLIMATE DATA TO 360 DAY YEARS ---------- ###
################################################################################
# given the climate dataframe produced in the climate processing file,
# creates 360 day years by forcing each month to be 30 days by either
# removing or repeating values
climate_to_30_day_months <- function(clim_df, start_year = 2014, end_year = 2022){
  dates_30_day_months <- generate_360_day_dates(start_year = start_year, end_year = end_year)
  n_years <- max(year(clim_df$date)) - min(year(clim_df$date)) + 1
  anom_360 <- rep(NA, 360 * n_years)
  temp_360 <- rep(NA, 360 * n_years)
  rollmean_360 <- rep(NA, 360 * n_years)
  vec_i = 1 # keep track of vector position
  months_31_days <- c(1, 3, 5, 7, 8, 10, 12)
  for(date_i in 1:nrow(clim_df)){ # keep track of date position
    if((clim_df$month[date_i] == 2) & (clim_df$day[date_i] == 28)){
      if(leap_year(year(clim_df$date[date_i]))){
        anom_360[vec_i:(vec_i+1)] <- clim_df$anom[date_i]
        temp_360[vec_i:(vec_i+1)] <- clim_df$temp[date_i]
        rollmean_360[vec_i:(vec_i+1)] <- clim_df$rollmean[date_i]
        vec_i = vec_i + 2
      }else{
        anom_360[vec_i:(vec_i+2)] <- clim_df$anom[date_i]
        temp_360[vec_i:(vec_i+2)] <- clim_df$temp[date_i]
        rollmean_360[vec_i:(vec_i+2)] <- clim_df$rollmean[date_i]
        vec_i = vec_i + 3
      }
    }else if((clim_df$month[date_i] %in% months_31_days) & (clim_df$day[date_i] == 31)){
      next
    }else{
      anom_360[vec_i] <- clim_df$anom[date_i]
      temp_360[vec_i] <- clim_df$temp[date_i]
      rollmean_360[vec_i] <- clim_df$rollmean[date_i]
      vec_i = vec_i + 1
    }
  }
  clim_df_360_day_years <- data.frame(dates = dates_30_day_months, anom = anom_360, temp = temp_360, rollmean = rollmean_360)
  return(clim_df_360_day_years)
}

# Function to generate a sequence of dates assuming a 360-day year (12 months of 30 days each).
# This is used to simulate a calendar without leap years and 31-day months.
# Args:
#   - start_year: The starting year of the date sequence.
#   - end_year: The ending year of the date sequence.
# Returns:
#   - A vector of dates formatted as "YYYY-MM-DD" for the entire period.
generate_360_day_dates <- function(start_year, end_year) {
  dates <- vector("list", length = (end_year - start_year + 1) * 12 * 30)  # Initialize list to store dates
  counter <- 1  # Counter to track the number of generated dates

  # Loop through each year
  for (year in start_year:end_year) {
    # Loop through each month
    for (month in 1:12) {
      # Loop through each day (assuming 30 days per month)
      for (day in 1:30) {
        # Handle February (only 28 days in this simplified calendar)
        if ((month == 2) & (day > 28)) {
          dates[[counter]] <- paste(year, sprintf("%02d", month), day, sep = "-")
        } else {
          dates[[counter]] <- format(as.Date(paste(year, month, day, sep = "-")), "%Y-%m-%d")
        }
        counter <- counter + 1  # Increment the counter
      }
    }
  }

  return(unlist(dates))  # Return the list of dates as a vector
}

#' Calculate the Day Difference Assuming a 360-Day Year
#'
#' This function calculates the difference between two dates, assuming a 360-day year
#' where each month is treated as having 30 days. This approach is commonly used in
#' financial calculations for simplicity.
#'
#' @param date1 The start date as a character string or `Date` object.
#' @param date2 The end date as a character string or `Date` object.
#'
#' @return The difference in days between `date1` and `date2` assuming a 360-day year.
#' @export
#'
#' @examples
#' # Example 1: Calculate difference between two dates
#' calculate_360_day_difference("2023-01-01", "2023-12-31")
#'
#' # Example 2: Use Date class objects
#' calculate_360_day_difference(as.Date("2023-01-01"), as.Date("2024-01-01"))
calculate_360_day_difference <- function(date1, date2) {
  # Convert dates to Date class if they aren't already
  date1 <- as.Date(date1)
  date2 <- as.Date(date2)

  # Extract year, month, and day components for each date
  y1 <- as.numeric(format(date1, "%Y"))
  m1 <- as.numeric(format(date1, "%m"))
  d1 <- as.numeric(format(date1, "%d"))

  y2 <- as.numeric(format(date2, "%Y"))
  m2 <- as.numeric(format(date2, "%m"))
  d2 <- as.numeric(format(date2, "%d"))

  # Calculate the difference in years, months, and days
  year_diff <- y2 - y1
  month_diff <- m2 - m1
  day_diff <- d2 - d1

  # Adjust month and day differences
  # Each month is treated as 30 days, and each year as 360 days
  total_days <- year_diff * 360 + month_diff * 30 + day_diff

  return(total_days)
}

