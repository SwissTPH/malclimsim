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


#' Generate a schedule of coverage values for a given period with decay
#'
#' This function generates a schedule of coverage values for the given time period,
#' filling missing values with previous non-zero coverage and applying a decay calculation.
#'
#' @param smc_cov A data frame containing the coverage data. It should have at least two columns:
#'   - `date_start`: The date of the coverage value.
#'   - `coverage`: The coverage value corresponding to that date.
#' @param months_30_days Logical; whether to generate dates in 30-day month format (defaults to `FALSE`).
#'   If `TRUE`, the sequence of dates will follow a 360-day year with 12 months of 30 days each.
#' @param years A numeric vector containing the range of years to generate the schedule for (e.g., `c(2014, 2022)`).
#' @param const A numeric value representing the decay constant used in the decay calculation (default is `-0.1806`).
#'
#' @return A data frame containing the following columns:
#'   - `dates`: The sequence of dates from the beginning to the end of the time period.
#'   - `SMC`: A binary value indicating the presence (1) or absence (0) of coverage on each date.
#'   - `cov`: The coverage value for each date.
#'   - `decay`: The decay values calculated for each date based on the SMC and the decay constant.
#'
#' @export
#'
#' @examples
#' smc_schedule_from_data(smc_cov = smc_cov_data, months_30_days = FALSE, years = c(2014, 2022))
smc_schedule_from_data <- function(smc_cov, months_30_days, years, const = -0.1806) {

  # Set start and end dates based on the input years
  start_date <- as.Date(paste0(years[1], "-01-01"))  # Start date of the time range (Jan 1 of the first year)
  end_date <- as.Date(paste0(years[length(years)], "-12-31"))  # End date of the time range (Dec 31 of the last year)

  # Ensure that `date_start` in `smc_cov` is in Date format
  smc_cov$date_start <- as.Date(smc_cov$date_start)

  # Create a sequence of all dates from the start date to the end date (by default 1 day interval)
  all_dates <- data.frame(date_start = seq(start_date, end_date, by = "day"))

  # If `months_30_days` is TRUE, generate dates using the 30-day month format (360-day year)
  if (months_30_days) {
    all_dates <- data.frame(date_start = generate_360_day_dates(years[1], years[length(years)]))
  }

  # Convert `date_start` in `all_dates` to Date format to ensure compatibility with `smc_cov`
  all_dates$date_start <- as.Date(all_dates$date_start)

  # Merge the sequence of dates with the original `smc_cov` data frame
  smc_cov_day <- all_dates %>%
    left_join(smc_cov, by = "date_start") %>%
    # Replace NA values in the coverage column with zero initially (indicating no coverage)
    mutate(coverage = ifelse(is.na(coverage), 0, coverage))

  # Here, we ensure that coverage is zero for any year/month without SMC intervention
  # We do NOT fill missing values with previous non-zero coverage.
  # Any NA after left_join is set to zero, and we do not propagate coverage forward.

  # Create the final schedule data frame, including a binary `SMC` indicator and the coverage values
  smc_df <- data.frame(dates = smc_cov_day$date_start, SMC = 0, cov = smc_cov_day$coverage)

  # Set the `SMC` column to 1 for dates where there is coverage (coverage > 0)
  smc_df$SMC[smc_df$cov > 0] <- 1

  # Calculate the decay values for each date using the `calc_decay_arr` function
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = const)

  # Return the final data frame with dates, SMC, coverage, and decay
  return(smc_df)
}
