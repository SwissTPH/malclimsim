

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
#'
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
      wk_inc_total = info$index$wk_inc_total
    ),
    state = c(
      info$index$month_inc_total,
      info$index$month_inc_A,
      info$index$month_inc_C,
      info$index$wk_inc_total,
      info$index$wk_inc_A,
      info$index$wk_inc_C
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
#'
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

