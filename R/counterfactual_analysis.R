#' Set eff_SMC Parameter to Zero
#'
#' This function takes a list of model parameters and sets the `eff_SMC` parameter to zero to simulate
#' a scenario without Seasonal Malaria Chemoprevention (SMC).
#'
#' @param param_inputs A list of parameter values for the model.
#' @return A list with the `eff_SMC` parameter set to zero.
#' @examples
#' no_smc_params <- set_eff_smc_to_zero(param_inputs)
set_eff_smc_to_zero <- function(param_inputs) {
  # Set the SMC effectiveness parameter to zero
  param_inputs$eff_SMC <- 0
  return(param_inputs)
}

#' Simulate Model with SMC
#'
#' This function simulates the model using the parameters that maximize the log posterior, which include SMC.
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @return A data frame containing the simulation output with SMC.
#' @export
simulate_with_smc <- function(results, start_date, end_date, model) {
  # Use the simulate_with_max_posterior_params function to simulate with SMC
  sim_data <- simulate_with_max_posterior_params(results, start_date, end_date, model)
  return(sim_data)
}

#' Simulate Model Without SMC
#'
#' This function simulates the model with `eff_SMC` set to zero, representing a scenario without SMC intervention.
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @return A data frame containing the simulation output without SMC.
#' @export
simulate_without_smc <- function(results, start_date, end_date, model) {
  # Extract the best-fit parameters
  max_posterior_params <- extract_max_posterior_params(results)

  # Update the parameters to set eff_SMC to zero
  param_inputs <- results$param_inputs
  updated_params <- update_param_list(param_inputs, max_posterior_params)
  no_smc_params <- set_eff_smc_to_zero(updated_params)

  # Run the model simulation without SMC
  simulation_output <- data_sim(model, no_smc_params, start_date, end_date,
                                month = TRUE, round = FALSE, save = FALSE,
                                month_unequal_days = FALSE)

  return(simulation_output)
}

#' Calculate Cases Averted by SMC
#'
#' This function calculates the number and percent of malaria cases averted due to SMC intervention by comparing
#' the simulated scenarios with and without SMC.
#'
#' @param with_smc_df Data frame containing the incidence data from the model simulated with SMC.
#' @param without_smc_df Data frame containing the incidence data from the model simulated without SMC.
#' @return A list containing the number of cases averted and the percentage of cases averted.
#' @export
#' @examples
#' cases_averted <- calculate_cases_averted(sim_with_smc, sim_without_smc)
calculate_cases_averted <- function(with_smc_df, without_smc_df) {

  # Ensure both data frames have the same dates for proper comparison
  if (!all(with_smc_df$date_ymd == without_smc_df$date_ymd)) {
    stop("Dates in with_smc_df and without_smc_df do not match.")
  }

  # Calculate the total number of cases for each scenario
  total_cases_with_smc <- sum(with_smc_df$inc, na.rm = TRUE)
  total_cases_without_smc <- sum(without_smc_df$inc, na.rm = TRUE)

  # Calculate the number and percentage of cases averted
  cases_averted <- total_cases_without_smc - total_cases_with_smc
  percent_averted <- (cases_averted / total_cases_without_smc) * 100

  return(list(
    cases_averted = cases_averted,
    percent_averted = percent_averted
  ))
}

#' Set SMC Coverage Level
#'
#' This function takes a list of model parameters and sets the `cov_SMC` parameter to a user-defined level.
#'
#' @param param_inputs A list of parameter values for the model.
#' @param coverage_level A numeric value (between 0 and 1) representing the desired coverage level. Default is 1 (100% coverage).
#' @return A list with the `cov_SMC` parameter updated to the specified coverage level.
#' @examples
#' full_coverage_params <- set_smc_coverage(param_inputs, coverage_level = 1)
set_smc_coverage <- function(param_inputs, coverage_level = 1) {
  # Update all values of cov_SMC to the desired level
  param_inputs$cov_SMC <- rep(coverage_level, length(param_inputs$cov_SMC))
  return(param_inputs)
}


#' Simulate Model with Full SMC Coverage
#'
#' This function simulates the model with `cov_SMC` set to a user-defined level (default is 100% coverage).
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @param coverage_level A numeric value between 0 and 1 for the desired SMC coverage. Default is 1 (100% coverage).
#' @return A data frame containing the simulation output with full or specified SMC coverage.
#' @export
#' @examples
#' sim_full_coverage <- simulate_with_full_coverage(results, start_date, end_date, model, coverage_level = 1)
simulate_with_full_coverage <- function(results, start_date, end_date, model, coverage_level = 1) {
  # Extract the best-fit parameters
  max_posterior_params <- extract_max_posterior_params(results)

  # Update the parameters with the extracted values
  param_inputs <- results$param_inputs
  updated_params <- update_param_list(param_inputs, max_posterior_params)

  # Set the SMC coverage level to full (or specified level)
  updated_params <- set_smc_coverage(updated_params, coverage_level)

  # Run the model simulation with the updated coverage level
  simulation_output <- data_sim(model, updated_params, start_date, end_date,
                                month = TRUE, round = FALSE, save = FALSE,
                                month_unequal_days = FALSE)

  return(simulation_output)
}


#' Calculate Confidence Intervals for eff_SMC
#'
#' This function calculates confidence intervals for eff_SMC by sampling
#' MCMC results, extracting parameter sets, updating them, and simulating outcomes.
#'
#' @param results The MCMC results object containing parameter samples.
#' @param param_inputs A list of parameter values for the model.
#' @param model The model function to simulate from.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param n_samples The number of parameter sets to sample from the MCMC results.
#' @return A data frame with confidence intervals for eff_SMC effectiveness.
#' @export
#' @examples
#' ci <- calculate_eff_smc_confidence_intervals(results, param_inputs, model, "2021-01-01", "2021-12-31", 100)
calculate_eff_smc_confidence_intervals <- function(results, param_inputs, model, start_date, end_date, n_samples = 100) {
  # Sample parameters from MCMC results
  sampled_params <- sample_params(results, n_samples)

  # Simulate outcomes with and without SMC for each parameter set
  simulations_with_smc <- simulate_models(model, param_inputs, sampled_params, start_date, end_date)
  simulations_without_smc <- simulate_models(model, param_inputs,
                                             lapply(sampled_params, set_eff_smc_to_zero),
                                             start_date, end_date)

  # Calculate cases averted for each simulation
  cases_averted_list <- mapply(calculate_cases_averted,
                               simulations_with_smc,
                               simulations_without_smc,
                               SIMPLIFY = FALSE)

  # Extract effectiveness estimates
  eff_smc_estimates <- sapply(cases_averted_list, function(x) x$percent_averted)

  # Calculate confidence intervals
  ci <- quantile(eff_smc_estimates, probs = c(0.025, 0.975))

  return(data.frame(
    lower_ci = ci[1],
    upper_ci = ci[2],
    median = median(eff_smc_estimates)
  ))
}

