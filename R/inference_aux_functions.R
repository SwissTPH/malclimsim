

#' Create a Parameter Transformation Function
#'
#' This function returns a transformation function that appends parameter values to a list. It is used for transforming model parameters into a format compatible with other functions.
#'
#' @param temp A numeric value representing the temperature parameter.
#' @param c_R_D A numeric value representing the decay rate.
#' @param SMC A numeric value representing the SMC (Scale of Mass Concentration).
#' @param decay A numeric value representing the decay parameter.
#' @param cov_SMC A numeric value representing the covariance for SMC.
#'
#' @return A function that takes a vector of parameters (`theta`) and combines them with the other provided values into a list.
#' @export
#' @examples
#' transform_function <- make_transform(temp = 25, c_R_D = 0.5, SMC = 0.1, decay = 0.03, cov_SMC = 0.05)
#' transformed_params <- transform_function(c(0.1, 0.2, 0.3))
make_transform <- function(temp, c_R_D, SMC, decay, cov_SMC) {
  function(theta) {
    c(list(temp = temp, c_R_D = c_R_D, SMC = SMC, decay = decay,
           cov_SMC = cov_SMC), as.list(theta))
  }
}

#' Create Index Mapping for Model Variables
#'
#' This function creates a list of indices for different model variables, organizing them into 'run' and 'state' categories based on the provided info.
#'
#' @param info A list containing index information, typically provided from a model setup.
#'
#' @return A list containing two elements, `run` and `state`, each with named indices for different model variables.
#'
#' @export
#' @examples
#' model_indices <- index(info)
index <- function(info) {
  list(
    run = c(
      month_inc_C = info$index$month_inc_C,
      month_inc_A = info$index$month_inc_A,
      month_inc_total = info$index$month_inc_total,
      wk_inc_C = info$index$wk_inc_C,
      wk_inc_A = info$index$wk_inc_A,
      wk_inc_total = info$index$wk_inc_total,
      prev_C_2 = info$index$prev_C_2,
      prev_A_2 = info$index$prev_A_2,
      prev_C_1 = info$index$prev_C_1,
      prev_A_1 = info$index$prev_A_1

    ),
    state = c(
      info$index$month_inc_total,
      info$index$month_inc_A,
      info$index$month_inc_C,
      info$index$wk_inc_total,
      info$index$wk_inc_C,
      info$index$wk_inc_A,
      info$index$prev_C_2,
      info$index$prev_A_2,
      info$index$prev_C_1,
      info$index$prev_A_1

    )
  )
}


#' Relate Time Steps to Observed Data for mcstate
#'
#' This function formats observed data into a form usable by the `mcstate` particle filter, adjusting for time steps (weeks or months) and ensuring the data is in the correct format for particle filtering.
#'
#' @param incidence_observed A data frame containing observed incidence data (cumulative incidence or other types of observations).
#' @param month A logical value indicating whether the data is provided in monthly intervals (`TRUE`) or weekly intervals (`FALSE`).
#' @param initial_time_obs An integer specifying the initial time step for monthly data (defaults to 0 for weekly data).
#' @param rate An integer specifying the rate (e.g., 7 for weekly data, 30 for monthly data) for resampling the time steps.
#'
#' @return A data frame with filtered data, including time intervals and the appropriate columns for particle filtering.
#' @export
#' @examples
#' # Assuming 'incidence_observed' is a data frame with your observed data
#' formatted_data <- relate_time_steps_to_observed_data(incidence_observed, month = FALSE, initial_time_obs = 0, rate = 7)
#'
filter_data <- function(incidence_observed, month = FALSE, initial_time_obs = 0, rate = 7) {
  # Check if the 'mcstate' package is loaded
  if (!requireNamespace("mcstate", quietly = TRUE)) {
    stop("The 'mcstate' package is required for this function.")
  }

  # Process for monthly data
  if (month) {
    filt_data <- mcstate::particle_filter_data(data = incidence_observed,
                                               time = "month_no",
                                               initial_time = initial_time_obs, rate = 30)

    filt_data$month_no_start <- (initial_time_obs):(nrow(incidence_observed) + initial_time_obs - 1)
    filt_data$month_no_end <- (initial_time_obs + 1) : (nrow(incidence_observed) + initial_time_obs)
    filt_data <- filt_data[-1,]

    #filt_data <- filt_data[!rowSums(apply(filt_data, 2, is.na)) > 0,]

  } else {
    # Process for weekly data
    filt_data <- mcstate::particle_filter_data(data = incidence_observed,
                                               time = "week_no",
                                               initial_time = 0, rate = 7)
    filt_data$week_no_start <- 0:(nrow(incidence_observed) - 1)
    filt_data$week_no_end <- 1:nrow(incidence_observed)
    filt_data <- filt_data[-1,]
  }

  # Remove the last row if it exceeds time steps
  filt_data <- filt_data[-nrow(filt_data),]

  return(filt_data)
}


#' #' Relate Time Steps to Observed Data for mcstate
#' #'
#' #' This function formats observed data into a form usable by the `mcstate` particle filter,
#' #' adjusting for time steps (weekly or monthly) and ensuring consistency for particle filtering.
#' #'
#' #' @param incidence_observed A data frame containing observed incidence data.
#' #' @param month Logical; TRUE for monthly data, FALSE for weekly data.
#' #' @param initial_time_obs Integer; initial time step (e.g., 0), used to align simulation and data.
#' #' @param rate Integer; time step rate (7 for weekly, 30 for monthly, etc.).
#' #'
#' #' @return A data frame formatted for particle filtering, with aligned start/end time steps.
#' #' @export
#' filter_data <- function(incidence_observed, month = FALSE, initial_time_obs = 0, rate = 7) {
#'   if (!requireNamespace("mcstate", quietly = TRUE)) {
#'     stop("The 'mcstate' package is required for this function.")
#'   }
#'
#'   time_col <- if (month) "month_no" else "week_no"
#'   prefix <- if (month) "month" else "week"
#'
#'   # Particle filter data
#'   filt_data <- mcstate::particle_filter_data(
#'     data = incidence_observed,
#'     time = time_col,
#'     initial_time = initial_time_obs,
#'     rate = rate
#'   )
#'
#'   n_obs <- nrow(incidence_observed)
#'
#'   # Add start and end columns with consistent naming
#'   filt_data[[paste0(prefix, "_no_start")]] <- initial_time_obs:(initial_time_obs + n_obs - 1)
#'   filt_data[[paste0(prefix, "_no_end")]] <- (initial_time_obs + 1):(initial_time_obs + n_obs)
#'
#'   # Drop the first row to align with intervals
#'   filt_data <- filt_data[-1, ]
#'
#'   # Drop last row if it exceeds time range (optional safety check)
#'   if (nrow(filt_data) > 0 && anyNA(filt_data[nrow(filt_data), ])) {
#'     filt_data <- filt_data[-nrow(filt_data), ]
#'   }
#'
#'   return(filt_data)
#' }

