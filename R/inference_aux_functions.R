

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
      prev_A_1 = info$index$prev_A_1,
      r_C = info$index$r_C,
      r_A = info$index$r_A

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
      info$index$prev_A_1,
      info$index$r_C,
      info$index$r_A

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

#' Update proposal matrix, start values, and MCMC settings for next stage
#'
#' @param results_obj    Result object returned by inf_run (with MCMC traces)
#' @param proposal_matrix Current proposal variance matrix
#' @param param_names    Character vector of parameter names (rownames of proposal_matrix)
#' @param S_prev         Integer: number of samples to use for variance-covariance extraction
#' @param draw_n         Integer: number of random draws for new start values per chain
#' @param shrink         Numeric: shrinkage factor for bounds (e.g. 0.8 retains 80% of interval)
#' @param stage          Character: which stage to pull from create_mcmc_params()
#' @param n_steps        Integer or NULL: if not NULL, overrides control_params$n_steps
#' @return A list with elements:
#'   * proposal_matrix : updated proposal variance matrix
#'   * start_values    : matrix of new start values (chains x parameters)
#'   * adaptive_params : list of adaptive settings for this stage
#'   * control_params  : list of control settings for this stage (with n_steps possibly overridden)
#' @export
update_inf_stage <- function(results_obj,
                             proposal_matrix,
                             param_names,
                             S_prev    = 3000,
                             draw_n     = 4,
                             shrink     = 0.8,
                             stage,
                             n_steps   = NULL) {
  # 1. Extract the posterior draws and new proposal‐matrix
  vcv_res <- extract_vcv(
    results_obj,
    S_prev      = S_prev,
    save        = FALSE,
    param_names = param_names
  )
  param_samples    <- vcv_res[[1]]   # matrix: n_params × S_prev
  proposal_matrix2 <- vcv_res[[2]]   # updated proposal variances

  # 2. Compute shrunken [min, max] for each parameter
  start_param_names <- rownames(param_samples)
  n_par             <- nrow(param_samples)
  param_min <- numeric(n_par)
  param_max <- numeric(n_par)

  for (j in seq_len(n_par)) {
    vals <- param_samples[j, ]
    mn   <- min(vals); mx <- max(vals)
    # shrink the _range_ around zero, preserving sign
    param_min[j] <- if (mn < 0) mn * shrink else mn * (1/shrink)
    param_max[j] <- if (mx > 0) mx * shrink else mx * (1/shrink)
    # ensure ordering
    if (param_min[j] > param_max[j]) {
      temp         <- param_min[j]
      param_min[j] <- param_max[j]
      param_max[j] <- temp
    }
  }
  names(param_min) <- names(param_max) <- start_param_names

  # 3. Draw new starting values
  start_values <- matrix(
    nrow   = draw_n,
    ncol   = n_par,
    dimnames = list(NULL, start_param_names)
  )
  for (j in seq_len(n_par)) {
    start_values[, j] <- runif(
      n   = draw_n,
      min = param_min[j],
      max = param_max[j]
    )
  }

  start_values <- t(start_values)

  # 4. Pull in adaptive/control settings for this stage
  mcmc_params     <- create_mcmc_params(stage = stage)
  adaptive_params <- mcmc_params$adaptive_params
  control_params  <- mcmc_params$control_params

  # 5. Optionally override n_steps
  if (!is.null(n_steps)) {
    control_params$n_steps <- n_steps
  }

  # 6. Return everything needed for the next run_stage() call
  list(
    proposal_matrix = proposal_matrix2,
    start_values    = start_values,
    adaptive_params = adaptive_params,
    control_params  = control_params
  )
}



