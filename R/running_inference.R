#' Generate Synthetic Incidence Data (Optional)
#'
#' Optionally generates synthetic incidence data from a model, with the ability to add Poisson noise.
#'
#' @param model A simulation model object used to generate synthetic data.
#' @param param_inputs A list of model parameters to use for simulation.
#' @param dates A vector of two dates (`start_date`, `end_date`) specifying the simulation period.
#' @param month Logical; whether to aggregate simulated data by calendar months.
#' @param month_unequal_days Logical; whether to weight months by unequal days during simulation.
#' @param noise Logical; if `TRUE`, adds Poisson noise to the simulated outputs.
#' @param seed Integer; random seed used for reproducible noise.
#' @param synthetic Logical; if `TRUE`, synthetic data will be generated. If `FALSE`, returns the input `incidence_df`.
#' @param incidence_df A data frame of observed incidence; returned unchanged if `synthetic = FALSE`.
#'
#' @return A data frame of synthetic or observed incidence data with columns such as `inc`, `inc_C`, `inc_A`.
#'
#' @export
generate_synthetic_data <- function(model, param_inputs, dates, month, month_unequal_days, noise, seed, synthetic, incidence_df) {
  if (synthetic) {
    if (month) {
      incidence_df <- data_sim(model, param_inputs, start_date = dates[1], end_date = dates[2],
                               month = TRUE, save = FALSE, month_unequal_days = month_unequal_days)
    } else {
      incidence_df <- data_sim(model, param_inputs, start_date = dates[1], end_date = dates[2],
                               month = FALSE, save = FALSE, month_unequal_days = month_unequal_days)
    }

    if (noise) {
      set.seed(seed)
      incidence_df$inc <- rpois(nrow(incidence_df), lambda = incidence_df$inc)
      incidence_df$inc_C <- rpois(nrow(incidence_df), lambda = incidence_df$inc_C)
      incidence_df$inc_A <- rpois(nrow(incidence_df), lambda = incidence_df$inc_A)
    }
  }
  return(incidence_df)
}

# Define comparison function
# generate_comparison_function <- function(month, age_for_inf, incidence_observed, include_prev) {
#   comparison_fn <- generate_incidence_comparison(month, age_for_inf, incidence_observed, include_prev)
#   return(comparison_fn)
# }


#' Generate Model-Observation Comparison Function
#'
#' Wraps a call to `generate_incidence_comparison()` to produce a likelihood or loss function for comparing model output to observed data.
#'
#' @param month Logical; whether the model and data are aggregated by calendar month.
#' @param age_for_inf Character; specifies the age group used for inference (e.g., `"inc_C"` or `"inc_A"`).
#' @param incidence_observed A data frame of observed incidence data.
#' @param include_prev Logical; whether to include prevalence data in the comparison function.
#'
#' @return A comparison function object that can be passed to `mcstate::pmcmc()` or other inference routines.
#'
#' @export
generate_comparison_function <- function(month, age_for_inf, incidence_observed, include_prev) {
  comparison_fn <- generate_incidence_comparison(month, age_for_inf, incidence_observed, include_prev)
  return(comparison_fn)
}



#' Initialize Observation Time Index for Incidence Data
#'
#' Aligns simulated and observed data by setting an initial time index (`month_no`) for observations.
#'
#' @param simulated_result A data frame from simulation output, with a `date_ymd` column.
#' @param incidence_df A data frame of observed incidence data, with a `date_ymd` column.
#'
#' @return The integer `0`, and updates the `month_no` column of `incidence_df` in-place (used for alignment).
#'
#' @export
initialize_observation_time <- function(simulated_result, incidence_df) {
  initial_time_obs <- incidence_df$month_no[simulated_result$date_ymd == incidence_df$date_ymd[1]]
  if (length(initial_time_obs) == 0) initial_time_obs <- 0
  incidence_df$month_no <- seq(initial_time_obs, (nrow(incidence_df) + initial_time_obs - 1))
  return(0)
}


#' Set Up Filtered Data for Likelihood Evaluation
#'
#' Prepares the observed incidence data for comparison to simulated data.
#'
#' @param incidence_observed A data frame of observed incidence values.
#' @param month Logical; whether the data is aggregated by month.
#' @param initial_time_obs Integer; starting time index for aligning simulated and observed data.
#'
#' @return A filtered data object suitable for use in likelihood comparison.
#'
#' @export
filter_data_setup <- function(incidence_observed, month, initial_time_obs) {
  filt_data <- filter_data(incidence_observed, month = month, initial_time_obs = initial_time_obs)
  return(filt_data)
}


#' Define Transformation Function for MCMC Simulation
#'
#' Wraps the construction of a transformation function for model covariates (e.g., temperature, SMC coverage).
#'
#' @param temp Numeric vector or matrix of temperature values.
#' @param c_R_D Numeric vector representing case reduction from treatment or diagnostics.
#' @param SMC Numeric vector of SMC coverage.
#' @param decay Numeric vector for SMC efficacy decay over time.
#' @param cov_SMC Numeric vector of SMC target coverage.
#'
#' @return A transformation function compatible with `mcstate` model input.
#'
#' @export
define_transformations <- function(temp, c_R_D, SMC, decay, cov_SMC) {
  transform_fn <- make_transform(temp = temp, c_R_D = c_R_D, SMC = SMC, decay = decay, cov_SMC = cov_SMC)
  return(transform_fn)
}


#' Define Priors and Proposals for MCMC Parameters
#'
#' Sets up priors, proposals, and parameter fixing for MCMC estimation, including validation and transformation of parameters.
#'
#' @param param_inputs A named list of initial parameter values (may include fixed and to-be-estimated parameters).
#' @param proposal_matrix A named matrix of proposal standard deviations (row and column names should match parameters).
#' @param params_to_estimate A character vector of parameter names to estimate using MCMC.
#' @param transform_fn A function that transforms the parameter list into model input format (typically created with `define_transformations()`).
#' @param param_priors Optional list of prior definitions. If `NULL`, defaults will be generated using `initialize_priors()`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{mcmc_pars}{An MCMC parameter object created using `mcstate::pmcmc_parameters$new()`.}
#'   \item{paramFix}{A named vector of parameters fixed during inference (not estimated).}
#'   \item{param_priors}{The full list of prior definitions, either provided or generated.}
#' }
#'
#' @details
#' This function performs several validation steps:
#' \itemize{
#'   \item Ensures all `params_to_estimate` exist in `param_inputs` and `proposal_matrix`.
#'   \item Verifies each prior includes a `min`, `max`, and a valid density function `prior`.
#'   \item Ensures transformed parameters match those to be estimated.
#'   \item Identifies and fixes parameters in `param_inputs` not selected for estimation.
#' }
#'
#' @examples
#' \dontrun{
#' transform_fn <- define_transformations(...)
#' result <- define_priors_and_proposals(
#'   param_inputs = init_params,
#'   proposal_matrix = proposal_mat,
#'   params_to_estimate = c("beta", "sigma"),
#'   transform_fn = transform_fn
#' )
#' }
#'
#' @export
define_priors_and_proposals <- function(param_inputs, proposal_matrix, params_to_estimate, transform_fn, param_priors = NULL) {
  if(is.null(param_priors)){
    param_priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)
  }

  # Ensure all parameters are valid
  missing_from_inputs <- setdiff(params_to_estimate, names(param_inputs))
  if (length(missing_from_inputs) > 0) {
    stop(paste("The following parameters are missing from param_inputs:", paste(missing_from_inputs, collapse = ", ")))
  }

  missing_from_proposals <- setdiff(params_to_estimate, rownames(proposal_matrix))
  if (length(missing_from_proposals) > 0) {
    stop(paste("The following parameters are missing from proposal_matrix:", paste(missing_from_proposals, collapse = ", ")))
  }

  # Validate priors
  for (param in params_to_estimate) {
    prior <- param_priors[[param]]
    if (is.null(prior$min) || is.null(prior$max)) {
      stop(paste("Missing min or max for parameter:", param))
    }
    if (!is.function(prior$prior)) {
      stop(paste("Prior function is missing or invalid for parameter:", param))
    }
  }

  # Apply transformations (if any)
  transformed_params <- transform_fn(param_inputs)

  # Validate transformed parameters
  for (param in params_to_estimate) {
    if (!param %in% names(transformed_params)) {
      stop(paste("Transformed parameter missing for:", param))
    }
  }

  # Define priors and proposals (remaining logic here)
  if(is.null(param_priors)){
    param_priors <- initialize_priors(param_inputs, proposal_matrix, params_to_estimate)
  }

  valid_params <- names(param_inputs)[names(param_inputs) %in% names(param_priors)]
  paramFix <- param_inputs[valid_params]
  paramFix <- paramFix[setdiff(names(paramFix), params_to_estimate)]
  paramFix <- paramFix[sapply(paramFix, function(x) length(x) == 1)]
  paramFix <- unlist(paramFix)

  mcmc_pars <- mcstate::pmcmc_parameters$new(param_priors, proposal_matrix, transform_fn)
  mcmc_pars <- mcmc_pars$fix(paramFix)

  return(list(mcmc_pars = mcmc_pars, paramFix = paramFix, param_priors = param_priors))
}




#' Define MCMC Control Settings
#'
#' Creates two sets of MCMC control objects for use with `mcstate::pmcmc`, one with a short run for adaptation,
#' and one for the main run with full settings.
#'
#' @param control_params A list with MCMC control settings including `n_steps`, `n_burnin`, `n_chains`, `n_workers`, and `n_threads_total`.
#' @param adaptive_params A list or object specifying the adaptive proposal settings (e.g., from `mcstate::adaptive_proposal()`).
#' @param save_trajectories Logical; whether to save state trajectories from the MCMC run.
#' @param rerun_n Integer; frequency of rerunning adaptive proposal fitting.
#' @param rerun_random Logical; whether rerun intervals are randomized.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{control1}{An initial short `pmcmc_control` object (e.g., for proposal adaptation).}
#'   \item{control2}{The full `pmcmc_control` object used for the main MCMC run.}
#' }
#'
#' @export
define_mcmc_control <- function(control_params, adaptive_params, save_trajectories, rerun_n, rerun_random) {
  control1 <- mcstate::pmcmc_control(n_steps = 10, n_burnin = 0, progress = TRUE, n_chains = control_params$n_chains)
  control2 <- mcstate::pmcmc_control(n_steps = control_params$n_steps, n_burnin = control_params$n_burnin,
                                     progress = TRUE, n_chains = control_params$n_chains,
                                     n_workers = control_params$n_workers, n_threads_total = control_params$n_threads_total,
                                     adaptive_proposal = adaptive_params, save_trajectories = save_trajectories,
                                     rerun_every = rerun_n, rerun_random = rerun_random)
  return(list(control1 = control1, control2 = control2))
}


#' Run MCMC Simulation
#'
#' Executes the MCMC simulation using `mcstate::pmcmc()` with two control stages. The first run is typically used to adapt the proposal distribution,
#' and the second run performs the full MCMC sampling.
#'
#' @param mcmc_pars MCMC parameter object (e.g., generated by `mcstate::pmcmc_parameters()`).
#' @param filter A particle filter object, typically created using `mcstate::particle_filter()`.
#' @param start_values Initial values for the MCMC chains.
#' @param control1 A `pmcmc_control` object for the initial adaptation stage.
#' @param control2 A `pmcmc_control` object for the full MCMC run.
#'
#' @return An MCMC object (of class `mcstate_pmcmc`) containing posterior samples and diagnostics.
#'
#' @export
run_mcmc_simulation <- function(mcmc_pars, filter, start_values, control1, control2) {
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control1)
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control2)
  return(mcmc_run)
}


#' Filter Incidence Data by Date Range
#'
#' Filters a time series of incidence data to only include rows within a specified date range.
#'
#' @param incidence_df A data frame with a `date_ymd` column (must be coercible to Date).
#' @param dates A vector of two dates (start and end) as `Date` objects or strings in `"YYYY-MM-DD"` format.
#'
#' @return A filtered data frame containing only rows with `date_ymd` within the specified range.
#'
#' @examples
#' df <- data.frame(date_ymd = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 100),
#'                  incidence = rpois(100, 5))
#' filtered <- filter_incidence_by_dates(df, c("2020-01-10", "2020-01-20"))
#'
#' @export
filter_incidence_by_dates <- function(incidence_df, dates) {
  # Convert 'date_ymd' column to Date if not already
  incidence_df$date_ymd <- as.Date(incidence_df$date_ymd)

  # Parse the start and end dates
  start_date <- as.Date(dates[1])
  end_date <- as.Date(dates[2])

  # Filter the incidence_df based on the date range
  filtered_df <- incidence_df[incidence_df$date_ymd >= start_date & incidence_df$date_ymd <= end_date, ]

  return(filtered_df)
}


#' Run MCMC Inference Simulation
#'
#' Generates synthetic incidence data if required, sets up MCMC parameters and priors,
#' defines model comparison functions, and runs the MCMC process.
#'
#' @param model The model object used for simulation and inference.
#' @param param_inputs List of input parameters for the model, including temperature, decay, and coverage details.
#' @param control_params Control parameters for the MCMC process including chain, burn-in, and step counts.
#' @param params_to_estimate Character vector of parameters to estimate during inference.
#' @param proposal_matrix Covariance matrix for the MCMC proposal distribution.
#' @param adaptive_params Parameters for adaptive proposal distribution if adaptive MCMC is used.
#' @param start_values Initial parameter values for starting the MCMC chains.
#' @param noise Logical; if TRUE, adds noise to the synthetic incidence data.
#' @param seed Numeric; sets the random seed for reproducibility.
#' @param month_unequal_days Logical; indicates if monthly data has unequal days.
#' @param dates Vector of start and end dates for the simulation period.
#' @param synthetic Logical; if TRUE, generates synthetic incidence data.
#' @param incidence_df Data frame of observed incidence data if not synthetic.
#' @param save_trajectories Logical; if TRUE, saves MCMC trajectories.
#' @param rerun_n Numeric; frequency for re-running MCMC proposals.
#' @param rerun_random Logical; if TRUE, re-runs MCMC randomly.
#' @param param_priors Optional prior list for MCMC parameters.
#' @param n_years_warmup Number of years to run prior to start of observation period.
#' @param obs_config A named list specifying observation model configuration.
#'        Use `make_obs_config()` to construct this.
#'
#' @return A list containing MCMC results, including posterior samples, fixed parameters, priors, and incidence data.
#' @export
inf_run <- function(model, param_inputs, control_params, params_to_estimate, proposal_matrix,
                    adaptive_params, start_values, noise = FALSE, seed = 24,
                    month_unequal_days = FALSE, dates, synthetic = TRUE, incidence_df = NULL,
                    save_trajectories = TRUE, rerun_n = Inf, rerun_random = FALSE,
                    param_priors = NULL, n_years_warmup = 3, obs_config) {

  # --- Extend inputs to include warm-up period ---
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
                                                 years_to_extend = n_years_warmup)
  extend_dates <- dates
  extend_dates[1] <- paste0(year(as.Date(dates[1])) - n_years_warmup, "-", format(as.Date(dates[1]), "%m-%d"))

  # --- Filter incidence data to extended date range ---
  incidence_df <- filter_incidence_by_dates(incidence_df, extend_dates)

  # --- Optionally generate synthetic incidence data ---
  if (synthetic) {
    incidence_df <- generate_synthetic_data(
      model, param_inputs_ext, dates,
      month = (obs_config$time == "month"),
      month_unequal_days = month_unequal_days,
      noise = noise, seed = seed,
      synthetic = TRUE,
      incidence_df = incidence_df
    )
  }

  # --- Define parameter transforms and priors ---
  transform_fn <- define_transformations(
    temp = param_inputs_ext$temp,
    c_R_D = param_inputs_ext$c_R_D,
    SMC = param_inputs_ext$SMC,
    decay = param_inputs_ext$decay,
    cov_SMC = param_inputs_ext$cov_SMC
  )

  priors_and_proposals <- define_priors_and_proposals(
    param_inputs, proposal_matrix, params_to_estimate, transform_fn
  )

  mcmc_pars <- priors_and_proposals$mcmc_pars
  paramFix <- priors_and_proposals$paramFix


  # --- Simulate data using the model ---
  simulated_result <- data_sim_for_inference(
    model, param_inputs = param_inputs_ext, dates = extend_dates,
    noise = FALSE, month = (obs_config$time == "month")
  )

  if(obs_config$time == "month"){
    simulated_result$month_no <- 0:(nrow(simulated_result) - 1)
  }else{simulated_result$week_no <- 0:(nrow(simulated_result) - 1)}

  simulated_result <- filter_incidence_by_dates(simulated_result, dates)

  all_dates <- as.Date(intersect(simulated_result$date_ymd, incidence_df$date_ymd))

  incidence_df <- incidence_df %>% dplyr::filter(date_ymd %in% all_dates)

  incidence_observed <- filter_incidence_by_dates(incidence_df, dates)[-1]

  if(obs_config$time == "month"){
    incidence_observed$month_no <- simulated_result$month_no
  }else{incidence_observed$week_no <- simulated_result$week_no}


  #if(month){
  #  incidence_observed$month_no <- 0:(nrow(incidence_observed) - 1)
  #}else{incidence_observed$week_no <- 0:(nrow(incidence_observed) - 1)}

  #if(month){
  #  incidence_observed$month_no <- simulated_result$month_no
  # }else{incidence_observed$week_no <- simulated_result$week_no}

  initial_time_obs <- initialize_observation_time(simulated_result, incidence_df)
  filt_data <- filter_data_setup(incidence_observed, (obs_config$time == "month"), initial_time_obs)

  # --- Set up observation comparison function ---
  comparison_fn <- generate_incidence_comparison(
    month = obs_config$time == "month",
    age_for_inf = obs_config$age_group,
    incidence_df = incidence_observed,
    include_prev = obs_config$include_prev,
    use_SMC_as_covariate = obs_config$use_SMC_as_covariate
  )

  # --- Initialize particle filter ---
  filter <- mcstate::particle_deterministic$new(
    data = filt_data,
    model = model,
    index = index,
    compare = comparison_fn
  )

  filter$run(c(param_inputs))

  # --- MCMC control ---
  control_settings <- define_mcmc_control(
    control_params, adaptive_params,
    save_trajectories, rerun_n, rerun_random
  )

  control1 <- control_settings$control1
  control2 <- control_settings$control2

  # --- Run MCMC simulation ---
  start_values <- reorder_start_values(start_values, priors_and_proposals$param_priors)
  mcmc_run <- run_mcmc_simulation(mcmc_pars, filter, start_values, control1, control2)
  coda_pars <- coda::as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))

  # --- Return results ---
  results <- list(
    mcmc_run = mcmc_run,
    coda_pars = coda_pars,
    paramFix = paramFix,
    param_inputs = param_inputs,
    incidence_df = incidence_df,
    model = model,
    param_priors = priors_and_proposals$param_priors,
    n_chains = control_params$n_chains
  )

  return(results)
}
# inf_run <- function(model, param_inputs, control_params, params_to_estimate, proposal_matrix,
#                     adaptive_params, start_values, noise = FALSE, seed = 24,
#                     month_unequal_days = FALSE, dates, synthetic = TRUE, incidence_df = NULL,
#                     save_trajectories = TRUE, rerun_n = Inf, rerun_random = FALSE,
#                     param_priors = NULL, n_years_warmup = 3, obs_config) {
#
#   # --- Extend inputs to include warm-up period ---
#   param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
#                                                  years_to_extend = n_years_warmup)
#   extend_dates <- dates
#   extend_dates[1] <- paste0(year(as.Date(dates[1])) - n_years_warmup, "-", format(as.Date(dates[1]), "%m-%d"))
#
#   # --- Filter incidence data to extended date range ---
#   incidence_df <- filter_incidence_by_dates(incidence_df, extend_dates)
#
#   # --- Optionally generate synthetic incidence data ---
#   if (synthetic) {
#     incidence_df <- generate_synthetic_data(
#       model, param_inputs_ext, dates,
#       month = (obs_config$time == "month"),
#       month_unequal_days = month_unequal_days,
#       noise = noise, seed = seed,
#       synthetic = TRUE,
#       incidence_df = incidence_df
#     )
#   }
#
#   # --- Define parameter transforms and priors ---
#   transform_fn <- define_transformations(
#     temp = param_inputs_ext$temp,
#     c_R_D = param_inputs_ext$c_R_D,
#     SMC = param_inputs_ext$SMC,
#     decay = param_inputs_ext$decay,
#     cov_SMC = param_inputs_ext$cov_SMC
#   )
#
#   priors_and_proposals <- define_priors_and_proposals(
#     param_inputs, proposal_matrix, params_to_estimate, transform_fn
#   )
#
#   mcmc_pars <- priors_and_proposals$mcmc_pars
#   paramFix <- priors_and_proposals$paramFix
#
#   # --- Set up observation comparison function ---
#   is_month_time <- obs_config$time == "month"
#   incidence_observed <- incidence_df[-1]  # remove date column if needed
#   comparison_fn <- generate_incidence_comparison(
#     month = is_month_time,
#     age_for_inf = obs_config$age_group,
#     incidence_df = incidence_observed,
#     include_prev = obs_config$include_prev,
#     use_SMC_as_covariate = obs_config$use_SMC_as_covariate,
#     log_link = obs_config$log_link
#   )
#
#   # --- Simulate data using the model ---
#   simulated_result <- data_sim_for_inference(
#     model, param_inputs = param_inputs_ext, dates = extend_dates,
#     noise = FALSE, month = (obs_config$time == "month")
#   )
#
#   simulated_result <- filter_incidence_by_dates(simulated_result, dates)
#   simulated_result$month_no <- 0:(nrow(simulated_result) - 1)
#
#   incidence_observed <- filter_incidence_by_dates(incidence_df, dates)[-1]
#   incidence_observed$month_no <- 0:(nrow(incidence_observed) - 1)
#
#   initial_time_obs <- initialize_observation_time(simulated_result, incidence_df)
#   filt_data <- filter_data_setup(incidence_observed, (obs_config$time == "month"), initial_time_obs)
#
#   # --- Initialize particle filter ---
#   filter <- mcstate::particle_deterministic$new(
#     data = filt_data,
#     model = model,
#     index = index,
#     compare = comparison_fn
#   )
#
#   filter$run(c(param_inputs))
#
#   # --- MCMC control ---
#   control_settings <- define_mcmc_control(
#     control_params, adaptive_params,
#     save_trajectories, rerun_n, rerun_random
#   )
#
#   control1 <- control_settings$control1
#   control2 <- control_settings$control2
#
#   # --- Run MCMC simulation ---
#   start_values <- reorder_start_values(start_values, priors_and_proposals$param_priors)
#   mcmc_run <- run_mcmc_simulation(mcmc_pars, filter, start_values, control1, control2)
#   coda_pars <- coda::as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))
#
#   # --- Return results ---
#   results <- list(
#     mcmc_run = mcmc_run,
#     coda_pars = coda_pars,
#     paramFix = paramFix,
#     param_inputs = param_inputs,
#     incidence_df = incidence_df,
#     model = model,
#     param_priors = priors_and_proposals$param_priors,
#     n_chains = control_params$n_chains
#   )
#
#   return(results)
# }



#' Extend Time-Varying Inputs Backwards in Time
#'
#' This function extends selected time-varying inputs by repeating the first year's data
#' backwards for a specified number of years. This is useful when simulating models that
#' require initialization over an extended period.
#'
#' @param param_inputs A named list of model input vectors, where some elements may be time-varying (e.g., daily values).
#' @param days_per_year Integer specifying the number of time steps (e.g., days) per year. Defaults to 360.
#' @param years_to_extend Number of years to prepend to each time-varying input by repeating the first year's data. Defaults to 2.
#'
#' @return A list with the same structure as `param_inputs`, where selected time-varying inputs have been extended backward in time.
#'
#' @details
#' The following keys are treated as time-varying and extended:
#' `"cov_SMC"`, `"SMC"`, `"decay"`, `"c_R_D"`, `"temp"`.
#' Other keys are left unchanged.
#'
#' @examples
#' inputs <- list(
#'   cov_SMC = rep(0.5, 360 * 3),
#'   temp = rep(25, 360 * 3),
#'   constant = 42
#' )
#' extended <- extend_time_varying_inputs(inputs, days_per_year = 360, years_to_extend = 1)
#' length(extended$cov_SMC)  # should be 360 * 4
#'
#' @export
extend_time_varying_inputs <- function(param_inputs, days_per_year = 360, years_to_extend = 2) {
  # Identify the time-varying inputs
  time_varying_keys <- c("cov_SMC", "SMC", "decay", "c_R_D", "temp")

  # Initialize an empty list to store the extended inputs
  extended_inputs <- list()

  # Loop through each parameter in the input list
  for (key in names(param_inputs)) {
    if (key %in% time_varying_keys) {
      # Repeat the first year's data backwards for the specified number of years
      first_year_data <- head(param_inputs[[key]], days_per_year)
      extended_data <- do.call(c, replicate(years_to_extend, first_year_data, simplify = FALSE))
      # Combine the extended data with the original data
      extended_inputs[[key]] <- c(extended_data, param_inputs[[key]])
    } else {
      # Leave constants unchanged
      extended_inputs[[key]] <- param_inputs[[key]]
    }
  }

  return(extended_inputs)
}

#' Run Model Simulation with Pre-Warm Period
#'
#' Extends time-varying covariate inputs and runs the model over a pre-warm period to allow system stabilization before the target simulation window.
#'
#' @param model A model object used for simulation (compatible with `data_sim()`).
#' @param param_inputs A list of model parameters, including time-varying covariates.
#' @param dates A vector of two dates: start and end of the target simulation window (as `Date` or string in `"YYYY-MM-DD"` format).
#' @param prewarm_years Integer; number of years to pre-warm the model before the actual simulation window. Defaults to 2.
#' @param days_per_year Integer; number of time steps per year (typically 360 or 365). Defaults to 360.
#'
#' @return A data frame of simulated results filtered to include only dates from the original `start_date` onward.
#'
#' @details
#' Time-varying parameters are extended backward in time by repeating the first year's values for `prewarm_years`.
#' The simulation is then run over this extended window, and the pre-warm results are discarded to avoid initialization artifacts.
#'
#' @examples
#' \dontrun{
#' sim_results <- run_simulation_with_prewarm(
#'   model = malaria_model,
#'   param_inputs = parameters,
#'   dates = c("2020-01-01", "2022-12-31"),
#'   prewarm_years = 2
#' )
#' }
#'
#' @export
run_simulation_with_prewarm <- function(model, param_inputs, dates, prewarm_years = 2, days_per_year = 360) {
  # Extend time-varying inputs
  extended_param_inputs <- extend_time_varying_inputs(param_inputs, days_per_year, prewarm_years)

  # Adjust the start date to account for the pre-warm period
  extended_start_date <- as.Date(dates[1]) - (prewarm_years * days_per_year)

  # Run the model simulation with extended inputs
  extended_results <- data_sim(model, param_inputs = extended_param_inputs,
                               start_date = extended_start_date,
                               end_date = dates[2],
                               month = FALSE)

  # Filter the results to only include dates starting from the original start date
  filtered_results <- extended_results[extended_results$date_ymd >= dates[1], ]

  return(filtered_results)
}


