#' Convert Date Range to Weekly Format (Sundays)
#'
#' This function takes a start date and an end date, and returns a vector of Dates
#' representing the Sunday of each week in the range.
#'
#' @param start_date The start date as a character string or `Date` object.
#' @param end_date The end date as a character string or `Date` object.
#'
#' @return A vector of `Date` objects representing each Sunday between the dates.
#' @export
#'
#' @examples
#' date_to_weeks("2022-01-01", "2022-03-01")
date_to_weeks <- function(start_date, end_date) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  # Align start_date to the nearest previous Sunday
  start_sunday <- lubridate::floor_date(start_date, unit = "week", week_start = 7)

  # Generate sequence of Sundays
  all_weeks <- seq(start_sunday, end_date, by = "1 week")

  return(all_weeks)
}


#' Convert Date Range to Monthly Format
#'
#' This function takes a start date and an end date, and converts the range into a vector of the first day
#' of each month in the range, formatted as "YYYY-MM".
#'
#' @param start_date The start date as a character string or `Date` object.
#' @param end_date The end date as a character string or `Date` object.
#'
#' @return A vector of `Date` objects representing the first day of each month within the given date range.
#' @export
#'
#' @examples
#' date_to_months("2022-01-01", "2022-06-01")
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

#' Convert Climate Data to 360-Day Years
#'
#' This function converts climate data to a 360-day calendar, with each month forced to 30 days.
#' Data is either removed or repeated to ensure the adjustment, useful for simplifying climate simulations.
#'
#' @param clim_df A data frame of climate data, which must include columns for `date`, `month`, `day`, `anom`, `temp`, and `rollrain`.
#' @param start_year The starting year of the conversion (default is 2014).
#' @param end_year The ending year of the conversion (default is 2022).
#'
#' @return A data frame with 360-day years, where each month has exactly 30 days.
#' @export
#'
#' @examples
#' climate_to_30_day_months(clim_df, start_year = 2014, end_year = 2022)
climate_to_30_day_months <- function(clim_df, start_year = 2014, end_year = 2022){
  dates_30_day_months <- generate_360_day_dates(start_year = start_year, end_year = end_year)
  n_years <- max(year(clim_df$date)) - min(year(clim_df$date)) + 1
  rain_360 <- rep(NA, 360 * n_years)
  anom_360 <- rep(NA, 360 * n_years)
  temp_360 <- rep(NA, 360 * n_years)
  rollrain_360 <- rep(NA, 360 * n_years)
  vec_i <- 1 # track vector position
  months_31_days <- c(1, 3, 5, 7, 8, 10, 12)

  for(date_i in 1:nrow(clim_df)){ # track date position
    if((clim_df$month[date_i] == 2) & (clim_df$day[date_i] == 28)){
      if(leap_year(year(clim_df$date[date_i]))){
        anom_360[vec_i:(vec_i+1)] <- clim_df$anom[date_i]
        temp_360[vec_i:(vec_i+1)] <- clim_df$temp[date_i]
        rollrain_360[vec_i:(vec_i+1)] <- clim_df$CumulativeRainfall[date_i]
        vec_i <- vec_i + 2
      } else {
        anom_360[vec_i:(vec_i+2)] <- clim_df$anom[date_i]
        temp_360[vec_i:(vec_i+2)] <- clim_df$temp[date_i]
        rollrain_360[vec_i:(vec_i+2)] <- clim_df$CumulativeRainfall[date_i]
        vec_i <- vec_i + 3
      }
    } else if((clim_df$month[date_i] %in% months_31_days) & (clim_df$day[date_i] == 31)){
      next
    } else {
      anom_360[vec_i] <- clim_df$anom[date_i]
      temp_360[vec_i] <- clim_df$temp[date_i]
      rollrain_360[vec_i] <- clim_df$CumulativeRainfall[date_i]
      vec_i <- vec_i + 1
    }
  }

  clim_df_360_day_years <- data.frame(date = dates_30_day_months,
                                      anom = anom_360,
                                      temp = temp_360,
                                      cumrain = rollrain_360)
  return(clim_df_360_day_years)
}




#' Generate Dates for a 360-Day Year Calendar
#'
#' This function generates a sequence of dates assuming a 360-day year, with each month consisting of 30 days.
#' Useful for simplified models where leap years and 31-day months are ignored.
#'
#' @param start_year The starting year of the date sequence.
#' @param end_year The ending year of the date sequence.
#'
#' @return A vector of dates formatted as "YYYY-MM-DD" for each day within the given range, assuming a 360-day year.
#' @export
#'
#' @examples
#' generate_360_day_dates(2014, 2022)
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
#' @param date1 The start date as a character string or Date object.
#' @param date2 The end date as a character string or Date object.
#'
#' @return The difference in days between date1 and date2 assuming a 360-day year.
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

#' Filter a dataset by a range of years
#'
#' This function filters a dataset to include only rows where the year extracted
#' from a specified date column falls within a given range of years. The date column
#' should be in the format "YYYY-MM-DD" as a string.
#'
#' @param data A data frame containing the dataset to filter.
#' @param date_column A string specifying the name of the column that contains dates as strings.
#' @param years_range A numeric vector specifying the range of years to include in the filtered dataset.
#'
#' @return A filtered data frame containing only rows within the specified year range.
#'
#' @examples
#' # Example dataset
#' met_360 <- data.frame(
#'   date = c("2012-01-01", "2015-06-15", "2020-03-22"),
#'   temp = c(26.3, 27.5, 28.1)
#' )
#'
#' # Filter data for the years 2014 to 2023
#' years_analysis <- 2014:2023
#' met_360_filtered <- filter_by_year(met_360, "date", years_analysis)
#'
#' @export
filter_by_year <- function(data, date_column, years_range) {
  # Extract the year from the specified date column
  data$year <- as.numeric(substr(data[[date_column]], 1, 4))

  # Filter the dataset to include only rows within the specified years range
  filtered_data <- subset(data, year %in% years_range)

  # Remove the temporary 'year' column
  filtered_data$year <- NULL

  # Return the filtered dataset
  return(filtered_data)
}

# edit_observation_function(){
#   package_code_path <- paste0(find.package("malclimsim"), "/R/")
#   utils::file.edit(paste0(package_code_path, "observation_functions.R"))
# }




#' Export Estimated and Fixed Parameters to LaTeX
#'
#' Extracts posterior summaries for estimated parameters and writes both estimated
#' and fixed parameters as LaTeX tables.
#'
#' @param results A list containing inference results with `coda_pars` and `param_inputs`.
#' @param params_to_estimate Character vector of parameter names that were inferred.
#' @param out_dir Directory to save the LaTeX files.
#' @param suffix Optional suffix to append to filenames (e.g., dataset name).
#' @param sigfig Number of significant figures to use (default = 3).
#'
#' @return No return value; writes LaTeX tables to disk.
#' @export
export_param_table_tex <- function(results, params_to_estimate, out_dir, suffix = NULL, sigfig = 3) {
  # Helper to build filenames with optional suffix
  suffix_str <- if (!is.null(suffix)) paste0("_", suffix) else ""

  # ---- Estimated Parameters ----
  post_samples <- results$coda_pars[, params_to_estimate, drop = FALSE]

  est_summary <- apply(post_samples, 2, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  })

  est_df <- data.frame(
    Parameter = colnames(est_summary),
    `Lower 2.5%` = signif(est_summary[1, ], sigfig),
    Median = signif(est_summary[2, ], sigfig),
    `Upper 97.5%` = signif(est_summary[3, ], sigfig),
    row.names = NULL
  )

  est_out_path <- file.path(out_dir, paste0("estimated_params", suffix_str, ".tex"))
  print(
    xtable::xtable(est_df, caption = "Estimated Parameters with 95\\% Credible Intervals",
                   label = paste0("tab:estimated_params", suffix_str)),
    file = est_out_path, include.rownames = FALSE
  )

  # ---- Fixed Parameters ----
  param_inputs <- results$param_inputs
  param_scalar <- param_inputs[sapply(param_inputs, function(x) is.numeric(x) && length(x) == 1)]
  fixed_params <- param_scalar[!names(param_scalar) %in% params_to_estimate]

  if (length(fixed_params) > 0) {
    fix_df <- data.frame(
      Parameter = names(fixed_params),
      Value = signif(unlist(fixed_params), sigfig),
      row.names = NULL
    )
    fix_out_path <- file.path(out_dir, paste0("fixed_params", suffix_str, ".tex"))
    print(
      xtable::xtable(fix_df, caption = "Fixed Model Parameters",
                     label = paste0("tab:fixed_params", suffix_str)),
      file = fix_out_path, include.rownames = FALSE
    )
  }
}


#' Save LaTeX Table of Scenario Estimates
#'
#' Writes a LaTeX-formatted table summarizing scenario estimates to file.
#'
#' @param summary_estimates A data frame with columns: Scenario, Lower_2.5, Median, Upper_97.5.
#' @param out_dir Directory to save the LaTeX file.
#' @param file_name Name of the .tex file (default = "scenario_summary.tex").
#'
#' @return No return value. Side effect: saves LaTeX file to disk.
#' @export
save_scenario_summary_tex <- function(summary_estimates, out_dir, file_name = "scenario_summary.tex") {
  out_path <- file.path(out_dir, file_name)
  latex_tbl <- xtable::xtable(
    summary_estimates,
    caption = "Estimated Cases Averted per Scenario (95\\% CrI)",
    label = "tab:scenario_summary"
  )
  print(latex_tbl, file = out_path, include.rownames = FALSE)
}

#' Save ggplot with Dynamic Filename
#'
#' Utility to save plots with a consistent naming scheme.
#'
#' @param plot_obj The ggplot object.
#' @param filename_stub Name without extension.
#' @param out_dir Directory to save to.
#' @param width Width in inches.
#' @param height Height in inches.
#'
#' @return Path to saved file (invisible).
#' @export
save_plot_dynamic <- function(plot_obj, filename_stub, out_dir, width = 10, height = 6) {
  fname <- file.path(out_dir, paste0(filename_stub, ".png"))
  ggplot2::ggsave(fname, plot = plot_obj, width = width, height = height, dpi = 300)
  invisible(fname)
}


#' Generate Vectors of Population Growth Scaling (r_C and r_A)
#'
#' Computes population multipliers over time for under-5 (r_C) and over-5 (r_A) groups,
#' using fixed daily exponential growth. Time steps can be monthly (30 days) or weekly (7 days).
#'
#' @param n Integer. Number of time points to generate (e.g. months or weeks).
#' @param month Logical. If `TRUE`, use 30-day months; if `FALSE`, use 7-day weeks.
#' @param growth_rate_C Daily growth rate for under-5 population. Default corresponds to 2.6% per year.
#' @param growth_rate_A Daily growth rate for over-5 population. Default corresponds to 2.6% per year.
#' @param init_r_C Initial multiplier for r_C (default = 1).
#' @param init_r_A Initial multiplier for r_A (default = 1).
#'
#' @return A data frame with columns: `timestep`, `r_C`, `r_A`
#'
#' @examples
#' # Monthly population growth over 60 months
#' get_population_scaling(n = 60, month = TRUE)
#'
#' # Weekly population growth over 120 weeks with custom growth rates
#' get_population_scaling(n = 120, month = FALSE, growth_rate_C = 1.00005, growth_rate_A = 1.00003)
#'
#' @export
get_population_scaling <- function(n,
                                   month = TRUE,
                                   growth_rate_C = 1.026^(1 / 360),
                                   growth_rate_A = 1.026^(1 / 360),
                                   init_r_C = 1,
                                   init_r_A = 1) {
  if (n <= 0) stop("n must be a positive integer")

  # Time steps in model days
  step_size <- if (month) 30 else 7
  model_days <- seq(0, by = step_size, length.out = n)

  # Compute exponential growth
  r_C <- init_r_C * growth_rate_C ^ model_days
  r_A <- init_r_A * growth_rate_A ^ model_days

  data.frame(
    timestep = seq_len(n),
    model_days = model_days,
    r_C = r_C,
    r_A = r_A
  )
}

#' Extend Time-Varying Inputs to Match Simulation Duration
#'
#' Ensures that all relevant time-varying parameters (e.g., temperature, SMC coverage)
#' are extended to match the required number of simulation days by prepending repeated
#' values from the first year of data.
#'
#' @param param_inputs A named list of model input vectors. Some entries should be time-varying inputs.
#' @param days_needed Integer specifying the total number of days required for simulation (e.g., `n_days`).
#' @param days_per_year Integer number of days per model year (default is 360).
#'
#' @return A list of parameter inputs where time-varying entries have been extended (and trimmed) to exactly `days_needed` days.
#'
#' @details
#' The following keys are treated as time-varying: `"cov_SMC"`, `"SMC"`, `"decay"`, `"c_R_D"`, `"temp"`.
#' These are extended by repeatedly prepending the first year of values until the total length is sufficient.
#' Non-time-varying keys are returned unchanged.
#'
#' @examples
#' inputs <- list(temp = rep(25, 360), constant = 0.1)
#' extended <- extend_time_varying_inputs_to_length(inputs, days_needed = 720)
#' length(extended$temp)  # should be 720
#'
#' @export
extend_time_varying_inputs_to_length <- function(param_inputs, days_needed, days_per_year = 360) {
  time_varying_keys <- c("cov_SMC", "SMC", "decay", "c_R_D", "temp")
  extended_inputs <- list()

  for (key in names(param_inputs)) {
    x <- param_inputs[[key]]
    if (key %in% time_varying_keys) {
      # Extract first year of data
      first_year <- head(x, days_per_year)
      # Repeat until we reach at least days_needed
      while (length(x) < days_needed) {
        x <- c(first_year, x)
      }
      # Trim to exact length
      extended_inputs[[key]] <- x[1:days_needed]
    } else {
      extended_inputs[[key]] <- x
    }
  }

  return(extended_inputs)
}


#' Prepend Prewarm Climate to Your True Analysis Series
#'
#' Takes your genuine, day‑1–to‑day‑N climate vectors (temp, rain, etc.)
#' and builds a single vector of length prewarm_years*360 + analysis_days,
#' where the first prewarm_years*360 days are just repeats of the first
#' 360 days of your real series, and the last analysis_days days are the
#' real data you passed in.
#'
#' @param param_inputs Named list of inputs, where time‑varying elements
#'   include "cov_SMC", "SMC", "decay", "c_R_D", "temp".
#'   Each of those should be a numeric vector of length = raw_analysis_days.
#' @param start_date Date at which your real series begins (e.g. "2018-01-01").
#' @param end_date   Date at which your real series ends   (e.g. "2023-12-31").
#' @param prewarm_years Integer number of 360‑day years to prepend.
#' @param days_per_year Number of days per model year (default: 360).
#'
#' @return A list like `param_inputs` but with each time‑varying vector
#'   extended to length = prewarm_years*360 + floor(analysis_days/7)*7.
#' @export
extend_inputs_with_prewarm <- function(param_inputs,
                                       start_date, end_date,
                                       prewarm_years,
                                       days_per_year = 360) {
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)

  # how many real days you need?
  raw_analysis_days <- calculate_360_day_difference(start_date, end_date) + 1
  analysis_days_use <- floor(raw_analysis_days / 7) * 7

  days_prewarm <- prewarm_years * days_per_year
  total_days   <- days_prewarm + analysis_days_use

  tv_keys <- c("cov_SMC","SMC","decay","c_R_D","temp")
  out <- list()

  for (k in names(param_inputs)) {
    x <- param_inputs[[k]]
    if (k %in% tv_keys) {
      if (length(x) < analysis_days_use) {
        stop(sprintf("Your '%s' vector is only length %d; you need at least %d days of real data.",
                     k, length(x), analysis_days_use))
      }
      real_segment <- x[1:analysis_days_use]
      baseline     <- real_segment[1:days_per_year]
      prewarm_seg  <- rep(baseline, length.out = days_prewarm)
      full_vec     <- c(prewarm_seg, real_segment)
      out[[k]]     <- full_vec[1:total_days]
    } else {
      out[[k]] <- x
    }
  }

  return(out)
}


#' Given a start-year and a 360-day offset, return (Year,Mon,Day) in 30‑day months
#'
#' @param day_index
#' @param start_year
#' @export
model_date_360 <- function(day_index, start_year) {
  # day_index: 1..N
  year_num     <- start_year + (day_index - 1) %/% 360
  day_in_year  <- ((day_index - 1) %% 360) + 1
  month_num    <- ((day_in_year - 1) %/% 30) + 1
  day_in_month <- ((day_in_year - 1) %% 30) + 1
  data.frame(
    year  = year_num,
    month = month_num,
    day   = day_in_month
  )
}

#' Generate weekly model dates (by 7-day blocks) in 360‑day calendar
#'
#' @param start_year integer, year of day 1
#' @param n_days     integer, total days (must be divisible by 7)
#'
#' @return A data.frame with one row per week: year, month, day (all numeric)
#' @export
date_to_weeks_360 <- function(start_year, n_days) {
  if (n_days %% 7 != 0) stop("n_days must be multiple of 7")
  wk_idxs <- seq(1, n_days, by = 7)
  md      <- model_date_360(wk_idxs, start_year)
  md$week_no <- seq_along(wk_idxs) - 1
  md
}

#' Impute Missing Climate Data Using Climatological Daily Means
#'
#' This function imputes missing dates in a climate time series using the average
#' value for each (month, day) combination across all prior years in the dataset.
#' It is useful for filling in short gaps at the end of a time series, such as
#' imputing values for December 2023 using earlier years' averages.
#'
#' @param data A data frame with at least the following columns:
#'   - `date` (Date): The observation date.
#'   - `CumulativeRainfall`, `anom`, `temp`: Numeric climate variables to impute.
#' @param end_date A `Date` object indicating the final date that should be present.
#'   The function will impute values from the last available date up to `end_date`.
#'   Assumes `date` column ends before this value.
#'
#' @return A data frame with the original data and appended imputed values, sorted by date.
#' @export
#'
#' @examples
#' # Suppose your climate data ends at 2023-12-01
#' # met_360_completed <- impute_climate_to_end_date(met_360, as.Date("2023-12-31"))
impute_climate_to_end_date <- function(data, end_date) {
  library(dplyr)
  library(lubridate)

  # Ensure required columns are present
  required_vars <- c("date", "CumulativeRainfall", "anom", "temp")
  missing_vars <- setdiff(required_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Missing required columns: ", paste(missing_vars, collapse = ", "))
  }

  # Add date components
  data <- data %>%
    mutate(month = month(date),
           day   = day(date),
           year  = year(date))

  # Find last available date
  last_date <- max(data$date, na.rm = TRUE)

  if (last_date >= end_date) {
    warning("No missing dates to impute; last date is already >= end_date.")
    return(data)
  }

  # Sequence of dates to impute
  missing_dates <- seq.Date(from = last_date + 1, to = end_date, by = "day")

  # Compute daily climatology (excluding the target year)
  clim_means <- data %>%
    filter(year < year(end_date)) %>%
    group_by(month, day) %>%
    summarise(
      CumulativeRainfall = mean(CumulativeRainfall, na.rm = TRUE),
      anom               = mean(anom, na.rm = TRUE),
      temp               = mean(temp, na.rm = TRUE),
      .groups = "drop"
    )

  # Create imputed data frame
  imputed_data <- data.frame(date = missing_dates) %>%
    mutate(month = month(date),
           day   = day(date),
           year  = year(date),
           week  = isoweek(date)) %>%
    left_join(clim_means, by = c("month", "day")) %>%
    select(date, CumulativeRainfall, anom, month, week, day, temp)

  # Bind and sort
  completed_data <- bind_rows(data, imputed_data) %>%
    arrange(date)

  return(completed_data)
}
