# Generate synthetic data if necessary
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
generate_comparison_function <- function(month, age_for_inf, incidence_observed) {
  comparison_fn <- generate_incidence_comparison(month, age_for_inf, incidence_observed)
  return(comparison_fn)
}

# Initialize observation time for alignment of observed and simulated data
initialize_observation_time <- function(simulated_result, incidence_df) {
  initial_time_obs <- simulated_result$month_no[simulated_result$date_ymd == incidence_df$month[1]]
  if (length(initial_time_obs) == 0) initial_time_obs <- 0
  incidence_df$month_no <- seq(initial_time_obs, (nrow(incidence_df) + initial_time_obs - 1))
  return(initial_time_obs)
}

# Setup filtered data
filter_data_setup <- function(incidence_observed, month, initial_time_obs) {
  filt_data <- filter_data(incidence_observed, month = month, initial_time_obs = initial_time_obs)
  return(filt_data)
}

# Define transformation function for MCMC
define_transformations <- function(temp, c_R_D, SMC, decay, cov_SMC) {
  transform_fn <- make_transform(temp = temp, c_R_D = c_R_D, SMC = SMC, decay = decay, cov_SMC = cov_SMC)
  return(transform_fn)
}

# Define priors and proposal parameters
define_priors_and_proposals <- function(param_inputs, proposal_matrix, params_to_estimate, transform_fn, param_priors = NULL) {
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

# Define MCMC control settings
define_mcmc_control <- function(control_params, adaptive_params, save_trajectories, rerun_n, rerun_random) {
  control1 <- mcstate::pmcmc_control(n_steps = 10, n_burnin = 0, progress = TRUE, n_chains = control_params$n_chains)
  control2 <- mcstate::pmcmc_control(n_steps = control_params$n_steps, n_burnin = control_params$n_burnin,
                                     progress = TRUE, n_chains = control_params$n_chains,
                                     n_workers = control_params$n_workers, n_threads_total = control_params$n_threads_total,
                                     adaptive_proposal = adaptive_params, save_trajectories = save_trajectories,
                                     rerun_every = rerun_n, rerun_random = rerun_random)
  return(list(control1 = control1, control2 = control2))
}

# Run the MCMC simulation
run_mcmc_simulation <- function(mcmc_pars, filter, start_values, control1, control2) {
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control1)
  mcmc_run <- mcstate::pmcmc(mcmc_pars, filter, initial = start_values, control = control2)
  return(mcmc_run)
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
#' @param month Logical; if TRUE, generates monthly data.
#' @param month_unequal_days Logical; indicates if monthly data has unequal days.
#' @param dates Vector of start and end dates for the simulation period.
#' @param age_for_inf Character; indicates type of age-specific incidence data ('total', 'sep_ages', 'total + sep_ages').
#' @param synthetic Logical; if TRUE, generates synthetic incidence data.
#' @param incidence_df Data frame of observed incidence data if not synthetic.
#' @param save_trajectories Logical; if TRUE, saves MCMC trajectories.
#' @param rerun_n Numeric; frequency for re-running MCMC proposals.
#' @param rerun_random Logical; if TRUE, re-runs MCMC randomly.
#'
#' @return A list containing MCMC results, including posterior samples, fixed parameters, priors, and incidence data.
#' @export
#'
#' @examples
#' inf_run(model, param_inputs, control_params, params_to_estimate, proposal_matrix, adaptive_params, start_values)
inf_run <- function(model, param_inputs, control_params, params_to_estimate, proposal_matrix,
                    adaptive_params, start_values, noise = FALSE, seed = 24, month = FALSE,
                    month_unequal_days = FALSE, dates, age_for_inf, synthetic = TRUE, incidence_df = NULL,
                    save_trajectories = TRUE, rerun_n = Inf, rerun_random = FALSE, param_priors = NULL) {

  # Generate synthetic data if necessary
  incidence_df <- generate_synthetic_data(model, param_inputs, dates, month, month_unequal_days, noise, seed, synthetic, incidence_df)

  # Define parameters and initialize transformation function
  transform_fn <- define_transformations(temp = param_inputs$temp, c_R_D = param_inputs$c_R_D, SMC = param_inputs$SMC,
                                         decay = param_inputs$decay, cov_SMC = param_inputs$cov_SMC)
  priors_and_proposals <- define_priors_and_proposals(param_inputs, proposal_matrix, params_to_estimate, transform_fn)
  mcmc_pars <- priors_and_proposals$mcmc_pars
  paramFix <- priors_and_proposals$paramFix

  # Set up comparison function and initial time
  incidence_observed <- incidence_df[-1]
  comparison_fn <- generate_comparison_function(month, age_for_inf, incidence_observed)
  simulated_result <- data_sim_for_inference(model, param_inputs = param_inputs, dates = dates, noise = FALSE, month = month)
  initial_time_obs <- initialize_observation_time(simulated_result, incidence_df)
  filt_data <- filter_data_setup(incidence_observed, month, initial_time_obs)

  # Define and set up filter for MCMC
  filter <- mcstate::particle_deterministic$new(data = filt_data, model = model, index = index, compare = comparison_fn)
  filter$run(c(param_inputs))

  # Define MCMC control settings
  control_settings <- define_mcmc_control(control_params, adaptive_params, save_trajectories, rerun_n, rerun_random)
  control1 <- control_settings$control1
  control2 <- control_settings$control2

  # Run MCMC
  start_values = reorder_start_values(start_values, priors_and_proposals$param_priors)
  mcmc_run <- run_mcmc_simulation(mcmc_pars, filter, start_values, control1, control2)
  coda_pars <- as.mcmc(cbind(mcmc_run$probabilities, mcmc_run$pars))

  # Collect results
  results <- list(mcmc_run = mcmc_run, coda_pars = coda_pars, paramFix = paramFix, param_inputs = param_inputs,
                  incidence_df = incidence_df, model = model, param_priors = priors_and_proposals$param_priors,
                  n_chains = control_params$n_chains)

  return(results)
}
