#' Generate an SMC deployment schedule with optional month-specific coverages
#'
#' @param start_date Character. First date of the schedule ("YYYY-MM-DD").
#' @param end_date   Character. Last  date of the schedule ("YYYY-MM-DD").
#' @param years      Integer vector. Years included in the schedule.
#' @param months_active Matrix (n_years x 12) with 1 = SMC month, 0 = not.
#' @param months_30_days Logical. Use 360-day calendar (every month = 30 days)?
#' @param coverage   Either **(i)** a single numeric (legacy behaviour) or
#'                   **(ii)** a numeric vector of length 12 *or* a named vector
#'                   whose names are month numbers (1-12) / abbreviations
#'                   ("Jan", "Feb", ...).  Values give coverage for each month.
#' @param smc_day_of_month Integer. Day of month SMC is administered.
#'
#' @return data.frame with columns:
#'   * dates - sequence of dates in the period
#'   * SMC   - 1 on round day, 0 otherwise
#'   * cov   - coverage applied (piecewise-constant blocks)
#'   * decay - efficacy decay computed from round days
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
  if (length(coverage) == 1) coverage <- rep(coverage, 12)

  if (!is.null(names(coverage)) && any(names(coverage) != "")) {
    cov_vec <- numeric(12)
    month_names <- c(month.abb, month.name)
    for (m in 1:12) {
      name_hits <- c(as.character(m), month.abb[m], month.name[m])
      hit <- intersect(name_hits, names(coverage))
      cov_vec[m] <- if (length(hit)) coverage[hit[1]] else NA_real_
    }
    if (anyNA(cov_vec))
      stop("Coverage missing for some months; supply all 12 values.")
    coverage <- cov_vec
  } else if (length(coverage) != 12) {
    stop("When coverage is a vector it must have length 12 or names for all months.")
  }

  ## -------- 3. Handle flexible months_active ------------------------------
  get_months_row <- function(y) {
    if (is.matrix(months_active)) {
      if (nrow(months_active) != length(years)) {
        stop("Number of rows in 'months_active' matrix must match length of 'years'")
      }
      return(months_active[which(years == y), ])
    } else if (is.list(months_active)) {
      y_str <- as.character(y)
      if (!y_str %in% names(months_active)) {
        stop(sprintf("Year %s not found in months_active list.", y_str))
      }
      row <- months_active[[y_str]]
      if (length(row) != 12) stop("Each months_active entry must have length 12.")
      return(row)
    } else {
      stop("'months_active' must be either a matrix or a named list.")
    }
  }

  ## -------- 4. Mark round days -------------------------------------------
  for (y in years) {
    active <- which(get_months_row(y) == 1)
    for (m in active) {
      first_day_m <- as.Date(sprintf("%04d-%02d-01", y, m))
      n_days_m <- as.integer(format(first_day_m + months(1) - 1, "%d"))
      round_day <- as.Date(sprintf("%04d-%02d-%02d",
                                   y, m, min(smc_day_of_month, n_days_m)))
      idx <- match(round_day, smc_df$dates, nomatch = 0L)
      if (idx) smc_df$SMC[idx] <- 1
    }
  }

  ## -------- 5. Apply coverage in 120-day window ---------------------------
  round_idx <- which(smc_df$SMC == 1)
  for (idx in round_idx) {
    this_month <- lubridate::month(smc_df$dates[idx])
    cov_value  <- coverage[this_month]
    window_end <- min(idx + 120, nrow(smc_df))
    smc_df$cov[idx:window_end] <- cov_value
  }

  ## -------- 6. Compute decay ----------------------------------------------
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = -0.1806)

  smc_df
}


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

  ## -------- 2. Resolve coverage per month --------------------------------
  if (length(coverage) == 1) coverage <- rep(coverage, 12)
  if (!is.null(names(coverage)) && any(names(coverage) != "")) {
    # named vector case
    cov_vec <- numeric(12)
    for (m in 1:12) {
      hits <- intersect(c(as.character(m), month.abb[m], month.name[m]),
                        names(coverage))
      cov_vec[m] <- if (length(hits)) coverage[hits[1]] else NA_real_
    }
    if (anyNA(cov_vec))
      stop("Coverage missing for some months; supply all 12 values.")
    coverage <- cov_vec
  } else if (length(coverage) != 12) {
    stop("When coverage is a vector it must have length 12 or have names for all months.")
  }

  ## -------- 3a. Helper: resolve months_active per year  -----------------
  get_months_row <- function(y) {
    if (is.matrix(months_active)) {
      if (nrow(months_active) != length(years))
        stop("Rows of 'months_active' must match length(years).")
      return(months_active[which(years == y), ])
    }
    if (is.list(months_active)) {
      y_str <- as.character(y)
      if (!y_str %in% names(months_active))
        stop("Year ", y, " not found in months_active list.")
      row <- months_active[[y_str]]
      if (length(row) != 12)
        stop("Each element of 'months_active' list must have length 12.")
      return(row)
    }
    stop("'months_active' must be either a matrix or a named list.")
  }

  ## -------- 3b. Helper: resolve smc_day_of_month per year/month ----------
  get_smc_day <- function(y, m) {
    # scalar
    if (length(smc_day_of_month) == 1) return(as.integer(smc_day_of_month))
    # length-12 vector
    if (is.numeric(smc_day_of_month) && length(smc_day_of_month) == 12) {
      return(as.integer(smc_day_of_month[m]))
    }
    # matrix of shape years x 12
    if (is.matrix(smc_day_of_month)) {
      if (nrow(smc_day_of_month) != length(years) ||
          ncol(smc_day_of_month) != 12) {
        stop("If matrix, 'smc_day_of_month' must be dim(years) x 12.")
      }
      return(as.integer(smc_day_of_month[which(years == y), m]))
    }
    # named list per year
    if (is.list(smc_day_of_month)) {
      y_str <- as.character(y)
      if (!y_str %in% names(smc_day_of_month))
        stop("Year ", y, " not found in smc_day_of_month list.")
      vec <- smc_day_of_month[[y_str]]
      if (length(vec) != 12)
        stop("Each element of 'smc_day_of_month' list must have length 12.")
      return(as.integer(vec[m]))
    }
    stop("'smc_day_of_month' must be a single value, length-12 vector, matrix, or named list.")
  }

  ## -------- 4. Mark SMC rounds --------------------------------------------
  for (y in years) {
    active_months <- which(get_months_row(y) == 1)
    for (m in active_months) {
      # first day of that month
      first_of_m <- as.Date(sprintf("%04d-%02d-01", y, m))
      # last day of that month
      last_of_m  <- as.Date(first_of_m + months(1) - days(1))
      # choose day within [1, n_days]
      desired_day <- get_smc_day(y, m)
      day_in_m    <- pmin(desired_day, as.integer(format(last_of_m, "%d")))
      round_day   <- as.Date(sprintf("%04d-%02d-%02d", y, m, day_in_m))
      idx <- match(round_day, smc_df$dates, nomatch = 0L)
      if (idx > 0L) smc_df$SMC[idx] <- 1
    }
  }

  ## -------- 5. Spread coverage over 120-day window -----------------------
  round_idx <- which(smc_df$SMC == 1)
  for (i in round_idx) {
    m <- lubridate::month(smc_df$dates[i])
    cov_val  <- coverage[m]
    window_end <- min(i + 120, nrow(smc_df))
    smc_df$cov[i:window_end] <- cov_val
  }

  ## -------- 6. Compute decay ----------------------------------------------
  smc_df$decay <- calc_decay_arr(smc_df$SMC, const = -0.1806)

  return(smc_df)
}


#' Generate SMC Coverage and Decay Schedule
#'
#' Constructs a daily schedule of SMC (Seasonal Malaria Chemoprevention) coverage and decay values
#' over a specified time period. Missing values are filled forward using the last observed non-zero
#' coverage, and decay is applied according to a specified exponential decay constant.
#'
#' @param smc_cov A data frame with at least two columns:
#'   \itemize{
#'     \item \code{date_start}: A \code{Date} or character vector indicating the start of coverage.
#'     \item \code{coverage}: Numeric values indicating the coverage on each date.
#'   }
#' @param months_30_days Logical. If \code{TRUE}, generates a synthetic 360-day calendar (12 months of 30 days).
#'        Defaults to \code{FALSE}.
#' @param years A numeric vector of length 2 specifying the range of years to include (e.g., \code{c(2014, 2022)}).
#' @param const A numeric scalar for the exponential decay constant (default: \code{-0.1806}).
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item \code{dates}: A daily sequence of dates across the specified time range.
#'     \item \code{SMC}: Binary indicator (1 if SMC applied on that day, 0 otherwise).
#'     \item \code{cov}: Propagated SMC coverage values after filling and masking non-SMC months.
#'     \item \code{decay}: Daily decay values calculated using the specified constant.
#'   }
#'
#' @examples
#' \dontrun{
#' smc_data <- data.frame(
#'   date_start = as.Date(c("2014-06-01", "2014-07-01", "2015-06-01")),
#'   coverage = c(0.8, 0.85, 0.9)
#' )
#' schedule <- smc_schedule_from_data(smc_cov = smc_data,
#' months_30_days = FALSE, years = c(2014, 2015))
#' head(schedule)
#' }
#'
#' @export
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

#' Load and clean raw SMC coverage data
#'
#' Cleans and formats SMC coverage data. Allows for flexibility in input:
#' the user can supply either a data frame directly or a path to a file
#' (Excel or RDS). If both `data` and paths are `NULL`, the function returns an error.
#'
#' @param data A data frame containing raw SMC coverage data (default = NULL).
#' @param path_to_excel Optional. Path to an Excel file containing SMC data (default = NULL).
#' @param path_to_rds Optional. Path to an RDS file containing SMC data (default = NULL).
#'
#' @return A cleaned data frame with columns: `date_start` and `coverage`.
#'         Only one row per Year-Month is retained (the first entry chronologically).
#' @export
load_clean_smc_data <- function(data = NULL, path_to_excel = NULL, path_to_rds = NULL) {
  # Load data depending on inputs
  if (is.null(data)) {
    if (!is.null(path_to_rds)) {
      data <- readRDS(path_to_rds)
    } else if (!is.null(path_to_excel)) {
      data <- readxl::read_excel(path_to_excel)
    } else {
      stop("Please provide either a data frame or a path to an Excel or RDS file.")
    }
  }

  # Clean and format the data
  data %>%
    dplyr::select(date_start, smc_couv_tot) %>%
    dplyr::rename(coverage = smc_couv_tot) %>%
    dplyr::mutate(YearMonth = format(date_start, "%Y-%m")) %>%
    dplyr::group_by(YearMonth) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-YearMonth)
}
