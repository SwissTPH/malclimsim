#' #' Generate an SMC Deployment Schedule and Calculate Efficacy Decay
#' #'
#' #' This function simulates an SMC (Seasonal Malaria Chemoprevention) schedule based on given parameters
#' #' and calculates efficacy decay over time.
#' #'
#' #' @param start_date A string specifying the start date of the SMC deployment (format: "YYYY-MM-DD").
#' #' @param end_date A string specifying the end date of the SMC deployment (format: "YYYY-MM-DD").
#' #' @param years A numeric vector of years for which SMC will be deployed.
#' #' @param months_active A matrix where each row corresponds to a year and columns represent months
#' #' (1 = active, 0 = inactive).
#' #' @param months_30_days A logical flag to simulate a 360-day calendar. Default is \code{FALSE}.
#' #' @param coverage A numeric value specifying the coverage rate of SMC. Default is \code{0.90}.
#' #' @param smc_day_of_month Determines which day of the month the SMC round begins.
#' #'
#' #' @return A dataframe with columns:
#' #' \item{dates}{The sequence of dates for the given period.}
#' #' \item{SMC}{Binary values indicating SMC deployment (1 = deployed, 0 = not deployed).}
#' #' \item{cov}{The coverage rate of SMC over the period.}
#' #' \item{decay}{The calculated efficacy decay over time.}
#' #'
#' #' @details
#' #' The function generates a sequence of dates and assigns SMC deployments based on active months
#' #' provided for each year. It also calculates the decay in efficacy over time for the deployed SMC.
#' #'
#' #' @examples
#' #' # Example usage:
#' #' years <- 2023:2024
#' #' months_active <- matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0,  # Year 2023
#' #'                           1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0), # Year 2024
#' #'                         nrow = 2, byrow = TRUE)
#' #' start_date <- "2023-01-01"
#' #' end_date <- "2024-12-31"
#' #' schedule <- gen_smc_schedule(start_date, end_date, years, months_active, coverage = 0.85)
#' #' head(schedule)
#' #'
#' #' @export
#' gen_smc_schedule <- function(start_date, end_date, years, months_active,
#'                              months_30_days = FALSE, coverage = 0.90,
#'                              smc_day_of_month = 1) {
#'
#'   if (months_30_days) {
#'     dates <- generate_360_day_dates(years[1], years[length(years)])
#'   } else {
#'     dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
#'   }
#'
#'   smc_df <- data.frame(
#'     dates = dates,
#'     SMC = 0,
#'     cov = 0
#'   )
#'
#'   # --- Assign SMC = 1 to selected day of active months ---
#'   for (i in seq_along(years)) {
#'     active_months <- which(months_active[i, ] == 1)
#'     year_i <- years[i]
#'     for (month in active_months) {
#'       # Ensure smc_day_of_month does not exceed the number of days in that month
#'       first_day <- as.Date(sprintf("%04d-%02d-01", year_i, month))
#'       days_in_month <- as.integer(format(first_day + months(1) - 1, "%d"))
#'       actual_day <- min(smc_day_of_month, days_in_month)
#'       smc_date <- as.Date(sprintf("%04d-%02d-%02d", year_i, month, actual_day))
#'
#'       idx <- which(smc_df$dates == smc_date)
#'       if (length(idx) == 1) {
#'         smc_df$SMC[idx] <- 1
#'       }
#'     }
#'   }
#'
#'   # --- Apply coverage for each SMC round (30-day window) ---
#'   smc_rounds <- which(smc_df$SMC == 1)
#'   for (idx in smc_rounds) {
#'     window_end <- min(idx + 120, nrow(smc_df))  # Avoid going past end of data
#'     smc_df$cov[idx:window_end] <- coverage
#'   }
#'
#'   # --- Decay is only calculated from SMC = 1 onward ---
#'   smc_df$decay <- calc_decay_arr(smc_df$SMC, const = -0.1806)
#'
#'   return(smc_df)
#' }

#' Generate an SMC deployment schedule with optional month-specific coverages
#'
#' @param start_date Character. First date of the schedule ("YYYY-MM-DD").
#' @param end_date   Character. Last  date of the schedule ("YYYY-MM-DD").
#' @param years      Integer vector. Years included in the schedule.
#' @param months_active Matrix (n_years × 12) with 1 = SMC month, 0 = not.
#' @param months_30_days Logical. Use 360-day calendar (every month = 30 days)?
#' @param coverage   Either **(i)** a single numeric (legacy behaviour) or
#'                   **(ii)** a numeric vector of length 12 *or* a named vector
#'                   whose names are month numbers (1–12) / abbreviations
#'                   ("Jan", "Feb", …).  Values give coverage for each month.
#' @param smc_day_of_month Integer. Day of month SMC is administered.
#'
#' @return data.frame with columns:
#'   * dates – sequence of dates in the period
#'   * SMC   – 1 on round day, 0 otherwise
#'   * cov   – coverage applied (piecewise-constant blocks)
#'   * decay – efficacy decay computed from round days
#' @export
gen_smc_schedule <- function(start_date, end_date, years, months_active,
                             months_30_days = FALSE,
                             coverage = 0.90,
                             smc_day_of_month = 1) {

  ## -------- 1. Build date skeleton ----------------------------------------
  dates <- if (months_30_days) {
    generate_360_day_dates(years[1], tail(years, 1))
  } else {
    seq(as.Date(start_date), as.Date(end_date), by = "day")
  }

  smc_df <- tibble::tibble(
    dates = dates,
    SMC   = 0,
    cov   = 0
  )

  ## -------- 2. Helper: resolve month-specific coverage --------------------
  # Convert single numeric → 12-length vector
  if (length(coverage) == 1) coverage <- rep(coverage, 12)

  # If vector is *named*, match against month number / abbrev
  if (!is.null(names(coverage)) && any(names(coverage) != "")) {
    # Re-index so that coverage_vec[1] = Jan, …, [12] = Dec
    cov_vec <- numeric(12)
    month_names <- c(month.abb, month.name)
    for (m in 1:12) {
      # accepted names: "1", "Jan", "January"
      name_hits <- c(
        as.character(m),
        month.abb[m],
        month.name[m]
      )
      hit <- intersect(name_hits, names(coverage))
      cov_vec[m] <- if (length(hit)) coverage[hit[1]] else NA_real_
    }
    if (anyNA(cov_vec))
      stop("Coverage missing for some months; supply all 12 values.")
    coverage <- cov_vec
  } else if (length(coverage) != 12) {
    stop("When coverage is a vector it must have length 12 or names for all months.")
  }

  ## -------- 3. Mark round days -------------------------------------------
  for (i in seq_along(years)) {
    active <- which(months_active[i, ] == 1)
    for (m in active) {
      y <- years[i]
      first_day_m <- as.Date(sprintf("%04d-%02d-01", y, m))
      n_days_m    <- as.integer(format(first_day_m + months(1) - 1, "%d"))
      round_day   <- as.Date(sprintf("%04d-%02d-%02d",
                                     y, m, min(smc_day_of_month, n_days_m)))
      idx <- match(round_day, smc_df$dates, nomatch = 0L)
      if (idx) smc_df$SMC[idx] <- 1
    }
  }

  ## -------- 4. Apply (month-specific) coverage in a 120-day window --------
  round_idx <- which(smc_df$SMC == 1)
  for (idx in round_idx) {
    this_month <- lubridate::month(smc_df$dates[idx])
    cov_value  <- coverage[this_month]
    window_end <- min(idx + 120, nrow(smc_df))
    smc_df$cov[idx:window_end] <- cov_value
  }

  ## -------- 5. Compute decay ---------------------------------------------
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = -0.1806)

  smc_df
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

  # Ensure that `date_start` in `smc_cov` is a character vector (for compatibility in joining)
  smc_cov$date_start <- as.character(smc_cov$date_start)

  # Create a sequence of all dates from the start date to the end date (by default 1 day interval)
  all_dates <- data.frame(date_start = as.character(seq(as.Date(start_date),
                                                        as.Date(end_date), by = "day")))

  # If `months_30_days` is TRUE, generate dates using the 30-day month format (360-day year)
  if(months_30_days) {
    all_dates <- data.frame(date_start = generate_360_day_dates(years[1], years[length(years)]))
  }

  # Merge the sequence of dates with the original `smc_cov` data frame
  smc_cov_day <- all_dates %>%
    left_join(smc_cov, by = "date_start") %>%
    # Replace NA values in the coverage column with zero initially (indicating no coverage)
    mutate(coverage = ifelse(is.na(coverage), 0, coverage))

  # Fill missing coverage values with the previous non-zero value using the `fill` function
  smc_cov_exp <- smc_cov_day %>%
    mutate(coverage = ifelse(coverage == 0, NA, coverage)) %>%
    fill(coverage, .direction = "down") %>%  # Propagate the last non-zero coverage forward
    # Replace any remaining NA values (those that had no coverage) with zero
    replace_na(list(coverage = 0))

  # Create the final schedule data frame, including a binary `SMC` indicator and the coverage values
  smc_df <- data.frame(dates = all_dates$date_start, SMC = 0, cov = smc_cov_exp$coverage)

  # Set the `SMC` column to 1 for dates where there is coverage (coverage > 0)
  smc_df$SMC[which(smc_cov_day$coverage > 0)] <- 1

  # New Step: Set coverage to zero for months without an SMC round
  # Convert dates to Date type for easier manipulation
  smc_df$dates <- as.Date(smc_df$dates)
  smc_cov$date_start <- as.Date(smc_cov$date_start)

  # Extract unique months where SMC occurs from the `smc_cov` data frame
  smc_months <- unique(format(smc_cov$date_start, "%Y-%m"))

  # Include the month after each SMC round as an extended coverage month
  smc_months_extended <- c(smc_months,
                           unique(format(smc_cov$date_start + months(1), "%Y-%m"))) %>%
    unique()  # Avoid duplicates

  # Set coverage to zero for dates not in `smc_months_extended`
  smc_df <- smc_df %>%
    mutate(month = format(dates, "%Y-%m")) %>%
    mutate(cov = ifelse(!(month %in% smc_months_extended), 0, cov)) %>%
    select(-month)  # Drop the month column for clarity

  # Calculate the decay values for each date using the `calc_decay_arr` function
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = const)

  # Return the final data frame with dates, SMC, coverage, and decay
  return(smc_df)
}

#' Generate SMC Coverage Matrix for Simulation
#'
#' Creates a monthly covariate matrix for SMC coverage from a deployment pattern.
#'
#' @param start_date Simulation start date ("YYYY-MM-DD").
#' @param end_date Simulation end date ("YYYY-MM-DD").
#' @param years Numeric vector of years to apply the SMC schedule.
#' @param month_pattern A 12-length binary vector representing monthly SMC rounds.
#' @param avg_cov Numeric average SMC coverage to apply.
#' @param exclude_years Optional vector of years to exclude from coverage.
#'
#' @return A data frame with `date_ymd` and `cov_SMC` columns.
#' @export
generate_smc_coverage_matrix <- function(start_date, end_date, years, month_pattern, avg_cov, exclude_years = NULL) {
  n_years <- (year(end_date) + 1) - year(start_date)
  months_active <- matrix(rep(month_pattern, n_years + 1), nrow = n_years + 1, byrow = TRUE)

  smc_schedule <- gen_smc_schedule(start_date, end_date, years, months_active, coverage = avg_cov)
  smc_schedule_monthly <- calculate_monthly_metrics(smc_schedule, exclude_years = exclude_years)

  data.frame(date_ymd = as.Date(smc_schedule_monthly$month), cov_SMC = smc_schedule_monthly$cov)
}

#' Compute monthly metrics from SMC schedule
#'
#' @param schedule Data frame with SMC schedule (must include `dates`, `SMC`, `cov`, `decay`)
#' @param exclude_years Vector of years to exclude (e.g., c(2023))
#'
#' @return Monthly summarized schedule with columns: month, SMC, cov, decay
#' @export
calculate_monthly_metrics <- function(schedule, exclude_years = NULL) {
  schedule %>%
    dplyr::mutate(date_ymd = format(as.Date(dates), "%Y-%m-01")) %>%
    dplyr::group_by(date_ymd) %>%
    dplyr::summarise(
      SMC = ifelse(sum(SMC, na.rm = TRUE) > 0, 1, 0),
      cov = sum(cov, na.rm = TRUE),
      decay = sum(decay, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!(lubridate::year(as.Date(date_ymd)) %in% exclude_years))
}



#' Compute weekly metrics from SMC schedule (epidemiological weeks, starting Sunday)
#'
#' @param schedule Data frame with SMC schedule (must include `dates`, `SMC`, `cov`, `decay`)
#' @param exclude_years Vector of years to exclude (e.g., c(2023))
#'
#' @return Weekly summarized schedule with columns: week_start (Date), SMC, cov, decay
#' @export
calculate_weekly_metrics <- function(schedule, exclude_years = NULL) {
  schedule %>%
    dplyr::mutate(
      date_ymd = lubridate::floor_date(as.Date(dates), unit = "week", week_start = 7)
    ) %>%
    dplyr::group_by(date_ymd) %>%
    dplyr::summarise(
      SMC = ifelse(sum(SMC, na.rm = TRUE) > 0, 1, 0),
      cov = sum(cov, na.rm = TRUE),
      decay = sum(decay, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!(lubridate::year(date_ymd) %in% exclude_years))
}
