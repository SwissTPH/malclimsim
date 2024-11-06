# Function to generate an SMC deployment schedule and calculate efficacy decay.
# This function simulates an SMC schedule based on given parameters and calculates efficacy decay.
# Args:
#   - start_date: The start date of the SMC deployment (format: "YYYY-MM-DD").
#   - end_date: The end date of the SMC deployment (format: "YYYY-MM-DD").
#   - years: A vector of years for which SMC will be deployed.
#   - months_active: A matrix where each row corresponds to a year and columns represent months (1 = active, 0 = inactive).
#   - months_30_days: A boolean flag to simulate a 360-day calendar (default = FALSE).
#   - coverage: The coverage rate of SMC (default = 0.90).
# Returns:
#   - A dataframe with dates, SMC schedule, and efficacy decay over time.
gen_smc_schedule <- function(start_date, end_date, years, months_active, months_30_days = FALSE, coverage = 0.90) {
  # Generate the sequence of dates
  if (months_30_days) {
    dates <- generate_360_day_dates(years[1], years[length(years)])
    smc_df <- data.frame(dates = dates, SMC = numeric(length(dates)), cov = rep(coverage, length(dates)))
  } else {
    dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
    smc_df <- data.frame(dates = format(dates, "%Y-%m-%d"), SMC = numeric(length(dates)), cov = rep(coverage, length(dates)))
  }

  # Loop through each year to assign SMC deployment to the active months
  for (i in 1:length(years)) {
    if (sum(months_active[i, ]) > 0) {  # Check if there are active months in the year
      active_months <- which(months_active[i, ] == 1)
      smc_days <- format(as.Date(paste(years[i], active_months, "01", sep = "-"), "%Y-%m-%d"))
      smc_df[match(smc_days, smc_df$dates), ]$SMC <- 1  # Set SMC to 1 on the first day of active months
    }
  }

  # Adjust the decay to start on the first day of each month (reset decay at the beginning of each month)
  decay_df <- smc_df
  start_SMC <- which(day(as.Date(decay_df$dates)) != 1)  # Get all dates that are not the first day of the month
  decay_df$SMC[start_SMC] <- 0  # Set SMC to 0 for non-first-day dates

  # Calculate the efficacy decay over time
  smc_df$decay <- calc_decay_arr(decay_df$SMC, const = -0.1806)

  return(smc_df)  # Return the SMC schedule with decay
}




smc_schedule_from_data <- function(smc_cov, months_30_days, years, const = -0.1806) {
  start_date <- as.Date(paste0(years[1], "-01-01"))
  end_date <- as.Date(paste0(years[length(years)], "-12-31"))
  smc_cov$date_start <- as.character(smc_cov$date_start)

  # Create a sequence of all days from 2014-01-01 to 2022-12-31
  all_dates <- data.frame(date_start = as.character(seq(as.Date(start_date),
                                                        as.Date(end_date), by = "day")))
  if(month_30_days){
    all_dates <- data.frame(date_start = generate_360_day_dates(years[1], years[length(years)]))
  }

  # Join the sequence of dates with the original data
  # NA values will appear for dates that were not in the original dataset
  smc_cov_day <- all_dates %>%
    left_join(smc_cov, by = "date_start") %>%
    # Replace NA values with zero initially
    mutate(coverage = ifelse(is.na(coverage), 0, coverage))

  # Fill in missing coverage values with the previous non-zero value
  smc_cov_exp <- smc_cov_day %>%
    mutate(coverage = ifelse(coverage == 0, NA, coverage)) %>%
    fill(coverage, .direction = "down") %>%
    # Replace any remaining NA values with zero
    replace_na(list(coverage = 0))

  smc_df <- data.frame(dates = all_dates$date_start, SMC = 0,
                       cov = smc_cov_exp$coverage)
  smc_df$SMC[which(smc_cov_day$coverage > 0)] <- 1
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = const)

  return(smc_df)
}

