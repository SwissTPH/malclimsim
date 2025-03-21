

#' Loading an Odin Model
#'
#' @param name name of the model written in Odin DSL to be loaded. This is a path to an R file,
#'
#' @return the loaded model to be used for simulation, inference, etc
#' @export
#'
#' @examples
#' load_model("model_det_1")
load_model <- function(name){
  file <- system.file("models", paste0(name, ".R"), package = "malclimsim")
  return(odin.dust::odin_dust(file))
}

# Helper function to model the decay of SMC (Seasonal Malaria Chemoprevention) efficacy over time.
# The decay follows a sigmoid curve where efficacy reduces as time (x) increases.
# Args:
#   - x: A numeric value representing the number of days since the last SMC treatment.
#   - const: The decay constant controlling the steepness of the sigmoid curve (default = -0.1806).
# Returns:
#   - A numeric value representing the remaining efficacy of SMC treatment at day 'x'.
decay_SMC <- function(x, const = -0.1806){
  return(0.90 / (1 + exp(const * (32.97 - x))))  # Sigmoid function for efficacy decay
}


calc_decay_arr <- function(SMC, decay_func = decay_SMC, const = -0.1806) {
  # Initialize decay array with zeros, length same as the input SMC
  decay_arr <- numeric(length(SMC))

  # Iterate over each point where SMC is administered
  i <- 1
  while (i <= length(SMC)) {
    if (SMC[i] > 0) {  # If SMC was administered
      # Find the index of the next SMC administration
      next_SMC <- which(SMC[(i + 1):length(SMC)] > 0)

      if (length(next_SMC) == 0) {
        # No subsequent SMC - apply decay until the end of the array
        end <- length(SMC)
      } else {
        # There is a subsequent SMC - apply decay until just before next SMC
        end <- i + next_SMC[1]
      }

      # Calculate decay values from index `i` to `end` (inclusive)
      decay_values <- decay_func(seq(0, end - i), const = const)

      # Ensure we do not exceed the boundaries of `decay_arr`
      decay_arr[i:end] <- decay_values[1:(end - i + 1)]
    }

    # Move to the next index
    i <- i + 1
  }

  # Return a decay array that has exactly the same length as `SMC`
  decay_arr <- decay_arr[1:length(SMC)]

  return(decay_arr)
}


# Function to run the malaria model simulation
# Arguments:
#   odin_mod: the model to simulate
#   pars: parameters for the model
#   time_start: initial time step
#   n_particles: number of particles (repetitions for stochastic simulation)
#   sim_time: total simulation time (days)
sim_mod <- function(odin_mod, pars, time_start, n_particles, sim_time){
  # Initialize the model with the provided parameters and particles
  model <- odin_mod$new(pars = pars, time = time_start, n_particles = n_particles)

  # Initialize an array to store simulation output for each time step
  x <- array(NA, dim = c(model$info()$len, n_particles, sim_time + 1))

  # Loop through each time step and run the model
  for (t in 0:sim_time) {
    x[ , , (t + 1)] <- model$run(t)  # Run the model at time t and store the result
  }

  # Return the simulation result array and the model object
  return(list(x, model))
}

#' Simulate Data from Model
#'
#' This function runs the simulation using the specified model and parameters,
#' and returns either monthly or weekly aggregated incidence data for adults (`inc_A`)
#' and children (`inc_C`), as well as total incidence (`inc`).
#'
#' @param model The model to be simulated.
#' @param param_inputs A list of parameters required for the simulation.
#' @param start_date The start date for the simulation, as a `Date` object or character string.
#' @param end_date The end date for the simulation, as a `Date` object or character string.
#' @param prewarm_years Integer, the number of years to simulate before `start_date` as a pre-warm period
#' (default is 2 years).
#' @param month Logical; if `TRUE`, aggregate data monthly; if `FALSE`, aggregate weekly.
#' @param round Logical; if `TRUE`, round the results to whole numbers.
#' @param save Logical; if `TRUE`, save the results to a file.
#' @param file Character; file name for saving the results, if `save` is `TRUE`.
#' @param month_unequal_days Logical; if `TRUE`, adjust aggregation for unequal days in months.
#' @return A data frame containing the incidence data aggregated by month or week.
#' @export
#' @examples
#' # Running the simulation and getting monthly data
#' inc_data <- data_sim(model, param_inputs, start_date = "2021-01-01", end_date = "2021-12-31", month = TRUE)
data_sim <- function(model, param_inputs, start_date, end_date,
                     prewarm_years = 2, month = FALSE, round = TRUE, save = TRUE, file = "",
                     month_unequal_days = FALSE, return_EIR = FALSE, return_compartments = FALSE){

  param_inputs <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
                                                 years_to_extend = prewarm_years)

  # Extend the start date backward by the pre-warm period
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Calculate the number of days in the full simulation period (prewarm + main simulation)
  n_days <- calculate_360_day_difference(prewarm_start_date, end_date) + 1


  # Run the simulation using the model and parameters
  results <- sim_mod(model, pars = param_inputs, time_start = 0,
                     n_particles = 1, sim_time = n_days)

  # Extract the simulation output and the model
  x <- results[[1]]
  mod <- results[[2]]

  # If monthly aggregation is selected
  if(month){
    if(month){
      # Create monthly indices for aggregation
      month_ind <- seq(1, n_days, by = 30)

      # Adjust for months with unequal days, if specified
      if(month_unequal_days){
        month_ind <- which(param_inputs$day_count == 0)
      }
    }

    inc_A <- x[mod$info()$index$month_inc_A,,][month_ind]
    inc_C <- x[mod$info()$index$month_inc_C,,][month_ind]

    month <- date_to_months(start_date = as.Date(prewarm_start_date), end_date = as.Date(end_date))

    n_months <- min(length(month), length(inc_A), length(inc_C))
    month <- month[1:n_months]
    inc_A <- inc_A[1:n_months]
    inc_C <- inc_C[1:n_months]

    month_no <- 0:(n_months - 1)

    if(return_EIR){
      EIR_monthly <- x[mod$info()$index$EIR_monthly,,][month_ind][main_sim_indices]
      EIR_monthly <- EIR_monthly[1:n_months]
      inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C, EIR_monthly = EIR_monthly)
    } else {
      inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C)
    }
  } else {
    # If weekly aggregation is selected
    wk_ind <- seq(1, n_days, by = 7)

    inc_A <- x[mod$info()$index$wk_inc_A,,][wk_ind]
    inc_C <- x[mod$info()$index$wk_inc_C,,][wk_ind]

    week <- date_to_weeks(start_date = as.Date(prewarm_start_date), end_date = as.Date(end_date))

    n_weeks <- min(length(week), length(inc_A), length(inc_C))
    week <- week[1:n_weeks]
    inc_A <- inc_A[1:n_weeks]
    inc_C <- inc_C[1:n_weeks]

    week_no <- 0:(n_weeks - 1)

    if(return_EIR){
      EIR_monthly <- x[mod$info()$index$EIR_monthly,,][wk_ind]
      EIR_monthly <- EIR_monthly[1:n_weeks]
      inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C, EIR_monthly = EIR_monthly)
    } else {
      inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C)
    }
  }


  # Round the results if specified
  if(round){
    inc_df[3:5] <- round(inc_df[3:5])
  }

  start_month <- which(inc_df$date_ymd == as.Date(format(as.Date(start_date), "%Y-%m-01")))
  end_month <- which(inc_df$date_ymd == as.Date(format(as.Date(end_date), "%Y-%m-01")))
  inc_df <- inc_df[start_month:end_month,]

  # Save the dataframe to a file if specified
  if(save){
    saveRDS(inc_df, paste0(dir, file))
  }

  return(inc_df)
}

data_sim <- function(model, param_inputs, start_date, end_date,
                     prewarm_years = 2, month = FALSE, round = TRUE, save = TRUE, file = "",
                     month_unequal_days = FALSE, return_EIR = FALSE, return_compartments = FALSE,
                     mu_transform_A = NULL, mu_transform_C = NULL, covariate_matrix = NULL) {

  param_inputs <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
                                             years_to_extend = prewarm_years)

  # Extend the start date backward by the pre-warm period
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Calculate the number of days in the full simulation period (prewarm + main simulation)
  n_days <- calculate_360_day_difference(prewarm_start_date, end_date) + 1

  # Run the simulation using the model and parameters
  results <- sim_mod(model, pars = param_inputs, time_start = 0,
                     n_particles = 1, sim_time = n_days)

  # Extract the simulation output and the model
  x <- results[[1]]
  mod <- results[[2]]

  if (month) {
    month_ind <- seq(1, n_days, by = 30)
    if (month_unequal_days) {
      month_ind <- which(param_inputs$day_count == 0)
    }

    inc_A <- x[mod$info()$index$month_inc_A,,][month_ind]
    inc_C <- x[mod$info()$index$month_inc_C,,][month_ind]

    month <- date_to_months(start_date = as.Date(prewarm_start_date), end_date = as.Date(end_date))
    n_months <- min(length(month), length(inc_A), length(inc_C))
    month <- month[1:n_months]
    inc_A <- inc_A[1:n_months]
    inc_C <- inc_C[1:n_months]

    month_no <- 0:(n_months - 1)

    if (return_EIR) {
      EIR_monthly <- x[mod$info()$index$EIR_monthly,,][month_ind][1:n_months]
      inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C, EIR_monthly = EIR_monthly)
    } else {
      inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C)
    }
  } else {
    wk_ind <- seq(1, n_days, by = 7)

    inc_A <- x[mod$info()$index$wk_inc_A,,][wk_ind]
    inc_C <- x[mod$info()$index$wk_inc_C,,][wk_ind]

    week <- date_to_weeks(start_date = as.Date(prewarm_start_date), end_date = as.Date(end_date))
    n_weeks <- min(length(week), length(inc_A), length(inc_C))
    week <- week[1:n_weeks]
    inc_A <- inc_A[1:n_weeks]
    inc_C <- inc_C[1:n_weeks]

    week_no <- 0:(n_weeks - 1)

    if (return_EIR) {
      EIR_monthly <- x[mod$info()$index$EIR_monthly,,][wk_ind][1:n_weeks]
      inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C, EIR_monthly = EIR_monthly)
    } else {
      inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C)
    }
  }

  # Merge with covariate matrix if provided
  if (!is.null(covariate_matrix)) {
    inc_df <- inc_df %>% left_join(covariate_matrix, by = c("date_ymd" = colnames(covariate_matrix)[1]))
  }

  # Apply transformations (if provided)
  if (!is.null(mu_transform_A)) {
    inc_df$inc_A <- mu_transform_A(inc_df, param_inputs)
  }
  if (!is.null(mu_transform_C)) {
    inc_df$inc_C <- mu_transform_C(inc_df, param_inputs)
  }

  # Update total incidence after transformations
  inc_df$inc <- inc_df$inc_A + inc_df$inc_C

  # Round the results if specified
  if (round) {
    inc_df[3:5] <- round(inc_df[3:5])
  }

  start_month <- which(inc_df$date_ymd == as.Date(format(as.Date(start_date), "%Y-%m-01")))
  end_month <- which(inc_df$date_ymd == as.Date(format(as.Date(end_date), "%Y-%m-01")))
  inc_df <- inc_df[start_month:end_month,]

  # Save the dataframe to a file if specified
  if (save) {
    saveRDS(inc_df, paste0(dir, file))
  }

  return(inc_df)
}





#' Simulate Data for Inference
#'
#' This function generates a simulated dataset for inference, based on the specified model and parameters.
#' It can include optional features such as monthly aggregation, unequal days in months, and adding noise to the simulated data.
#'
#' @param model The model to use for simulation. This should be a valid model object that `data_sim` can process.
#' @param param_inputs A list or vector of parameters required by the model for simulation.
#' @param dates A vector of two dates (`start_date` and `end_date`) specifying the time range for the simulation.
#' @param noise Logical. If `TRUE`, random noise is added to the simulated incidence data using a negative binomial distribution.
#' @param month Logical. If `TRUE`, incidence data is aggregated by month. Defaults to `FALSE`.
#' @param month_unequal_days Logical. If `TRUE`, the function accounts for months with unequal days in the simulation. Defaults to `FALSE`.
#'
#' @return A data frame containing simulated incidence data, formatted according to the output of `data_sim`.
#'         If `noise` is enabled, the `inc` column will include added noise.
#'
#' @details
#' - When `month = TRUE`, the simulation aggregates incidence by month, and the `month_unequal_days` parameter can control
#'   whether or not to adjust for months with different numbers of days.
#' - If `noise = TRUE`, a negative binomial distribution is used to add random noise to the incidence values.
#'   The size parameter for the negative binomial distribution is fixed at 100.
#'
#' @examples
#' # Example usage
#' model <- some_model_object
#' param_inputs <- list(beta = 0.3, gamma = 0.1)
#' dates <- c("2022-01-01", "2022-12-31")
#'
#' # Simulate data with no noise, not aggregated by month
#' sim_data <- data_sim_for_inference(model, param_inputs, dates, noise = FALSE, month = FALSE)
#'
#' # Simulate data with noise, aggregated by month
#' sim_data_with_noise <- data_sim_for_inference(model, param_inputs, dates, noise = TRUE, month = TRUE)
#'
#' @seealso
#' - `data_sim`: The underlying function used for simulation.
#'
#' @export
data_sim_for_inference <- function(model, param_inputs,
                                   dates, noise = FALSE, month = FALSE,
                                   month_unequal_days = FALSE){

  # incidence dataframe in specific format from data_sim
  if(month){
    incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                             end_date = dates[2], month = TRUE, save = FALSE,
                             month_unequal_days = month_unequal_days)
  }else{incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                                 end_date = dates[2], month = FALSE, save = FALSE,
                                 month_unequal_days = month_unequal_days)}

  if(noise){
    set.seed(seed)
    #incidence_df$inc <- rpois(nrow(incidence_df), lambda = incidence_df$inc)
    incidence_df$inc <- rnbinom(nrow(incidence_df), mu = incidence_df$inc, size = 100)
  }
  return(incidence_df)
}

#' Extract Parameters with Maximum Log Posterior
#'
#' This function finds and extracts the parameter set corresponding to the maximum log posterior value
#' from the MCMC results.
#'
#' @param results The MCMC results object containing parameter samples.
#' @return A named vector of parameter values corresponding to the maximum log posterior.
#' @export
#' @examples
#' # Assuming `results` contains the MCMC output
#' max_posterior_params <- extract_max_posterior_params(results)
extract_max_posterior_params <- function(results) {
  coda_pars <- results$coda_pars

  # Find the index of the row with the maximum log posterior
  max_posterior_index <- which.max(coda_pars[ , "log_posterior"])

  # Extract the row with the maximum log posterior
  max_posterior_row <- coda_pars[max_posterior_index, ]

  # Remove log_prior, log_likelihood, and log_posterior columns
  param_values <- max_posterior_row[!names(max_posterior_row) %in%
                                      c("log_prior", "log_likelihood", "log_posterior")]

  return(param_values)
}

#' Update Parameters in List with Values from Extracted Vector
#'
#' This function updates values in a list of parameters (e.g., `param_inputs`)
#' with values from a named vector of parameters obtained from an MCMC run.
#'
#' @param param_inputs A list of parameters where each element is a named parameter.
#' @param param_values A named vector of parameters, typically from the maximum
#' log posterior, to update values in `param_inputs`.
#' @return An updated list of parameters with replaced values where matches exist.
#' @export
#' @examples
#' # Assuming `param_inputs` is a list of parameters and `params_at_max_posterior`
#' # is a named vector with values to update
#' updated_params <- update_param_list(param_inputs, params_at_max_posterior)
update_param_list <- function(param_inputs, param_values) {
  # Loop over each element in the named vector 'param_values'
  for (param_name in names(param_values)) {
    # Check if the parameter exists in the list
    if (param_name %in% names(param_inputs)) {
      # Update the list's parameter with the value from 'param_values'
      param_inputs[[param_name]] <- param_values[[param_name]]
    }
  }
  return(param_inputs)
}

#' Simulate Model Using Parameters with Maximum Log Posterior
#'
#' This function takes the MCMC results, extracts the parameter values corresponding to the maximum
#' log posterior, updates the parameter list, and runs a simulation from the model. Additionally,
#' the function includes a pre-warm period, where the simulation starts earlier than the specified
#' `start_date` to allow the model to stabilize before the output period.
#'
#' @param results The MCMC results object containing parameter samples, including `log_posterior`.
#' @param start_date The start date for the desired simulation output, as a `Date` object or character string.
#' @param end_date The end date for the simulation, as a `Date` object or character string.
#' @param model The model function to simulate from. This function should accept parameters, start/end dates,
#'              and additional arguments for the simulation.
#' @param prewarm_years Integer, the number of years to simulate before `start_date` as a pre-warm period
#'                      (default is 2 years).
#' @param days_per_year Integer, the number of days in a year for the simulation (default is 360 days).
#'
#' @return A data frame containing the simulation results, filtered to include only the period
#'         from `start_date` to `end_date`.
#' @export
#'
#' @examples
#' # Assuming `results` contains the MCMC output, and `data_sim` is the simulation function
#' simulation_output <- simulate_with_max_posterior_params(
#'   results = results,
#'   start_date = "2021-01-01",
#'   end_date = "2021-12-31",
#'   model = data_sim,
#'   prewarm_years = 3,
#'   days_per_year = 360
#' )
simulate_with_max_posterior_params <- function(results, start_date, end_date, model, prewarm_years = 2, days_per_year = 360) {

  param_inputs <- results$param_inputs

  # Setup for allowing model to run some prior to inference (comparing to observations)
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
                                                 years_to_extend = prewarm_years)

  # Extract parameters with maximum log posterior
  max_posterior_params <- extract_max_posterior_params(results)

  # Update parameter list with the extracted parameters
  updated_params <- update_param_list(param_inputs_ext, max_posterior_params)

  # Extend the start date backward by the pre-warm period
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Run the model simulation for the pre-warm period + desired period
  extended_simulation_output <- data_sim(
    model = model,
    param_inputs = updated_params,
    start_date = prewarm_start_date,
    end_date = end_date,
    month = TRUE,
    round = FALSE,
    save = FALSE,
    month_unequal_days = FALSE
  )

  # Filter the simulation results to include only the desired period
  simulation_output <- extended_simulation_output[extended_simulation_output$date_ymd >= as.Date(start_date), ]

  return(simulation_output)
}



#' Sample Parameters from MCMC Results
#'
#' This function randomly selects `n` rows from the MCMC parameter output and extracts the
#' parameter values from each row.
#'
#' @param results The MCMC results object containing parameter samples.
#' @param n The number of samples to extract.
#' @return A list of named vectors, each containing one set of sampled parameters.
#' @examples
#' # Assuming `results` contains the MCMC output
#' sampled_parameters <- sample_params(results, 100)
sample_params <- function(results, n) {
  # Randomly sample `n` rows
  sampled_indices <- sample(1:nrow(results$coda_pars), size = n, replace = TRUE)
  sampled_params <- lapply(sampled_indices, function(idx) {
    # Extract row and remove columns for log_prior, log_likelihood, and log_posterior
    params <- results$coda_pars[idx, ]
    params <- params[!names(params) %in% c("log_prior", "log_likelihood", "log_posterior")]
    return(params)
  })
  return(sampled_params)
}

#' Simulate Models Using Sampled Parameters
#'
#' This function runs `n` model simulations, each using a different set of parameters sampled
#' from the MCMC results.
#'
#' @param model The model function to simulate from.
#' @param sampled_params A list of parameter sets to use for each simulation.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param param_inputs The initial list of parameters to update for each simulation.
#' @return A list of data frames, each containing the simulation output for a different parameter set.
#' @examples
#' # Assuming `sampled_params` is a list of sampled parameter sets
#' simulations <- simulate_models(model = data_sim, param_inputs = param_inputs,
#'                                sampled_params = sampled_params,
#'                                start_date = "2021-01-01", end_date = "2021-12-31")
# simulate_models <- function(model, param_inputs, sampled_params, start_date, end_date) {
#   simulations <- lapply(sampled_params, function(params) {
#     # Make a copy of param_inputs for each simulation to prevent overwriting
#     current_params <- param_inputs
#
#     # Update the copied parameters with the sampled parameters
#     updated_params <- update_param_list(current_params, params)
#
#     # Run the simulation with the updated parameters
#     sim_data <- data_sim(model, updated_params, start_date = start_date, end_date = end_date,
#                          month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
#     return(sim_data)
#   })
#   return(simulations)
# }
simulate_models <- function(model, param_inputs, sampled_params, start_date, end_date, prewarm_years = 2, days_per_year = 360) {
  # Extend parameter inputs to accommodate the prewarm period
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = days_per_year, years_to_extend = prewarm_years)

  # Calculate the prewarm start date
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Run simulations with prewarming
  simulations <- lapply(sampled_params, function(params) {
    # Make a copy of param_inputs for each simulation to prevent overwriting
    current_params <- param_inputs_ext

    # Update the copied parameters with the sampled parameters
    updated_params <- update_param_list(current_params, params)

    # Run the simulation for the prewarm period and the desired period
    extended_sim_data <- data_sim(
      model = model,
      param_inputs = updated_params,
      start_date = prewarm_start_date,
      end_date = end_date,
      month = TRUE,
      round = FALSE,
      save = FALSE,
      month_unequal_days = FALSE
    )

    # Filter the simulation results to include only the desired period
    sim_data <- extended_sim_data[extended_sim_data$date_ymd >= as.Date(start_date), ]

    return(sim_data)
  })

  return(simulations)
}

#' Calculate Quantiles for Incidence from Multiple Simulations
#'
#' This function calculates the 1st and 99th quantiles of the simulated incidence data across multiple simulations,
#' for each group (`inc_A`, `inc_C`, and `inc`).
#'
#' @param simulations A list of data frames, each containing the simulation output for different parameter sets.
#' @return A data frame containing the quantiles for each group (`inc_A`, `inc_C`, and `inc`) for each time point.
#' @examples
#' # Assuming `simulations` is a list of data frames with simulated incidence data
#' incidence_quantiles <- calculate_incidence_quantiles(simulations)
calculate_incidence_quantiles <- function(simulations) {
  # Extract incidence data for each group (inc_A, inc_C, inc) across all simulations
  inc_A_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc_A))
  inc_C_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc_C))
  inc_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc))

  # Calculate quantiles along rows (for each time point) for each group
  quantiles_inc_A <- apply(inc_A_data, 1, function(x) {
    c(quantile(x, probs = 0.005), quantile(x, probs = 0.995))
  })
  quantiles_inc_C <- apply(inc_C_data, 1, function(x) {
    c(quantile(x, probs = 0.005), quantile(x, probs = 0.995))
  })
  quantiles_inc <- apply(inc_data, 1, function(x) {
    c(quantile(x, probs = 0.005), quantile(x, probs = 0.995))
  })

  # Format the output into a data frame
  quantiles_df <- data.frame(
    date_ymd = simulations[[1]]$date_ymd,
    inc_A_q005 = quantiles_inc_A[1, ],
    inc_A_q995 = quantiles_inc_A[2, ],
    inc_C_q005 = quantiles_inc_C[1, ],
    inc_C_q995 = quantiles_inc_C[2, ],
    inc_q005 = quantiles_inc[1, ],
    inc_q995 = quantiles_inc[2, ]
  )

  return(quantiles_df)
}



#' Simulate Compartments Over Time with Prewarm Period
#'
#' This function simulates the compartments of an epidemiological model over a specified time period,
#' including an optional prewarm period to allow the model to reach equilibrium before the main simulation period.
#'
#' @param model An epidemiological model object used for simulation.
#' @param param_inputs A named list or vector of model parameters to be used in the simulation.
#' @param start_date A character string or Date object specifying the start date of the main simulation (e.g., "2014-01-01").
#' @param end_date A character string or Date object specifying the end date of the simulation (e.g., "2022-12-31").
#' @param prewarm_years Integer specifying the number of years to prewarm the model before the main simulation (default: 2).
#' @param days_per_year Integer specifying the number of days in a model year (default: 360).
#'
#' @return A data frame containing the simulated values of each compartment (SC, EC, IC, etc.)
#' over the simulation period with a weekly time resolution.
#'
#' @details The function simulates the dynamics of susceptible, exposed, infected, treated, and recovered
#' compartments for both children and adults. It first runs a prewarm period (if specified) to allow the model to stabilize,
#' before simulating the main study period. The results are aggregated at weekly intervals.
#'
#' @export
compartments_sim <- function(model, param_inputs, start_date, end_date, prewarm_years = 2, days_per_year = 360) {

  # Extend parameter inputs to accommodate the prewarm period
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = days_per_year, years_to_extend = prewarm_years)

  # Calculate the prewarm start date
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Calculate the total number of days including prewarm period
  total_days <- calculate_360_day_difference(prewarm_start_date, end_date) - 1

  # Run the simulation using the model and extended parameters
  results <- sim_mod(model, pars = c(param_inputs_ext), time_start = 0,
                     n_particles = 1, sim_time = total_days)

  # Extract the simulation output and the model
  x <- results[[1]]
  mod <- results[[2]]

  # Select weekly indices from the simulation results (excluding prewarm)
  full_wk_ind <- seq(1, total_days, by = 7)
  start_index <- calculate_360_day_difference(prewarm_start_date, start_date)
  wk_ind <- full_wk_ind[full_wk_ind >= start_index]

  # Ensure that extracted compartments match the length of wk_ind
  extract_compartment <- function(index_name) {
    values <- x[mod$info()$index[[index_name]],,]
    if (length(values) >= max(wk_ind)) {
      return(values[wk_ind])  # Ensure selection does not exceed available data
    } else {
      stop(paste("Mismatch in compartment size for:", index_name))
    }
  }

  # Extract compartment data using helper function
  SC <- extract_compartment("SC")
  EC <- extract_compartment("EC")
  IC <- extract_compartment("IC")
  TrC <- extract_compartment("TrC")
  RC <- extract_compartment("RC")
  SA <- extract_compartment("SA")
  EA <- extract_compartment("EA")
  IA <- extract_compartment("IA")
  TrA <- extract_compartment("TrA")
  RA <- extract_compartment("RA")
  PC <- extract_compartment("P_C")
  PA <- extract_compartment("P_A")
  EIR <- extract_compartment("EIR2")
  eff_SMC_cov <- extract_compartment("SMC_effect_2")
  prev_total_with_R <- extract_compartment("prev_total_1")
  prev_C_with_R <- extract_compartment("prev_C_1")
  prev_A_with_R <- extract_compartment("prev_A_1")
  prev_total_no_R <- extract_compartment("prev_total_2")
  prev_C_no_R <- extract_compartment("prev_C_2")
  prev_A_no_R <- extract_compartment("prev_A_2")
  SMC_effect <- extract_compartment("SMC_effect_2")
  mu_SE_C <- extract_compartment("mu_SE_C_2")
  mu_SE_A <- extract_compartment("mu_SE_A_2")
  X <- extract_compartment("X2")
  X_I <- extract_compartment("X_I")
  X_AP <- extract_compartment("X_AP")
  X_ASP <- extract_compartment("X_ASP")

  # Total population size
  P <- PC + PA

  # Generate correct weekly dates for the compartments
  compart_dates <- seq(as.Date(start_date), as.Date(end_date), by = "7 days")

  # Ensure the length of compart_dates matches the number of extracted weekly values
  min_length <- min(length(compart_dates), length(SC))

  # Trim all vectors to match the shortest length
  compart_dates <- compart_dates[1:min_length]
  SC <- SC[1:min_length]
  EC <- EC[1:min_length]
  IC <- IC[1:min_length]
  TrC <- TrC[1:min_length]
  RC <- RC[1:min_length]
  SA <- SA[1:min_length]
  EA <- EA[1:min_length]
  IA <- IA[1:min_length]
  TrA <- TrA[1:min_length]
  RA <- RA[1:min_length]
  PC <- PC[1:min_length]
  PA <- PA[1:min_length]
  P <- P[1:min_length]
  EIR <- EIR[1:min_length]
  eff_SMC_cov <- eff_SMC_cov[1:min_length]
  prev_total_with_R <- prev_total_with_R[1:min_length]
  prev_C_with_R <- prev_C_with_R[1:min_length]
  prev_A_with_R <- prev_A_with_R[1:min_length]
  prev_total_no_R <- prev_total_no_R[1:min_length]
  prev_C_no_R <- prev_C_no_R[1:min_length]
  prev_A_no_R <- prev_A_no_R[1:min_length]
  SMC_effect <- SMC_effect[1:min_length]
  mu_SE_C <- mu_SE_C[1:min_length]
  mu_SE_A <- mu_SE_A[1:min_length]
  X <- X[1:min_length]
  X_I <- X_I[1:min_length]
  X_AP <- X_AP[1:min_length]
  X_ASP <- X_ASP[1:min_length]

  # Create a data frame with compartment values
  compart_df <- data.frame(
    date = compart_dates,
    SC = SC,
    EC = EC,
    IC = IC,
    TrC = TrC,
    RC = RC,
    SA = SA,
    EA = EA,
    IA = IA,
    TrA = TrA,
    RA = RA,
    PC = PC,
    PA = PA,
    P = P,
    EIR = EIR,
    eff_SMC_cov = eff_SMC_cov,
    prev_total_with_R = prev_total_with_R,
    prev_C_with_R = prev_C_with_R,
    prev_A_with_R = prev_A_with_R,
    prev_total_no_R = prev_total_no_R,
    prev_C_no_R = prev_C_no_R,
    prev_A_no_R = prev_A_no_R,
    SMC_effect = SMC_effect,
    mu_SE_C = mu_SE_C,
    mu_SE_A = mu_SE_A,
    X = X,
    X_I = X_I,
    X_AP = X_AP,
    X_ASP = X_ASP
  )

  return(compart_df)
}

#' Simulate and Extract Yearly Prevalence Estimates for a District
#'
#' This function runs a compartmental model simulation for a specified district,
#' extracts yearly prevalence estimates, and returns a summarized data frame.
#'
#' @param model An epidemiological model object used for simulation.
#' @param results A list containing MCMC or model fitting results for the district,
#'        including parameter estimates and prior information.
#' @param start_date A character string or Date object specifying the start date of the simulation (e.g., "2014-01-01").
#' @param end_date A character string or Date object specifying the end date of the simulation (e.g., "2022-12-31").
#' @param district_name A character string specifying the name of the district.
#' @param param_inputs A named vector of parameters and corresponding values. If NULL, uses parameters that maximize likelihood in results.
#'
#' @return A data frame containing yearly average prevalence estimates for the given district.
#'         The data frame includes columns for `year`, various prevalence metrics (e.g., `prev_total_with_R`, `prev_C_with_R`),
#'         and the corresponding `District` name.
#'
#' @details
#' - The function first extracts the **maximum likelihood parameters** from the model fitting results.
#' - It updates the **parameter inputs** with these best estimates.
#' - The model is then **simulated** for the specified district and time period.
#' - The output is aggregated into **yearly averages** for all prevalence metrics.
#' - The district name is appended to the final data frame to facilitate later merging.
#'
#' @examples
#' # Example usage for the Koumra district
#' start_date <- "2014-01-01"
#' end_date <- "2022-12-31"
#' koumra_prev <- simulate_prevalence_for_district(model_simp, koumra_results, start_date, end_date, "Koumra")
#'
#' # View output
#' head(koumra_prev)
#'
#' @seealso
#' - `compartments_sim()`: Function used internally to simulate the model.
#' - `update_param_list()`, `extract_max_posterior_params()`: Functions to extract and update model parameters.
#'
#' @export
simulate_prevalence_for_district <- function(model, results, start_date, end_date, district_name, param_inputs = NULL) {

  if(is.null(param_inputs)){
    # Extract max likelihood parameters
    max_ll_params <- extract_max_posterior_params(results$results)

    # Update parameter inputs with max likelihood estimates
    param_inputs <- update_param_list(results$results$param_inputs, max_ll_params)
  }

  # Run the compartment simulation for the district
  compart_df <- compartments_sim(model = model, param_inputs = param_inputs,
                                 start_date = start_date, end_date = end_date)

  # Ensure date column is in Date format and extract year
  compart_df <- compart_df %>%
    mutate(date = as.Date(date), year = year(date))

  # Compute yearly averages for prevalence metrics
  compart_df_year_avg <- compart_df %>%
    group_by(year) %>%
    summarise(across(starts_with("prev"), mean, na.rm = TRUE), .groups = "drop")

  # Add district name for identification
  compart_df_year_avg <- compart_df_year_avg %>%
    mutate(District = district_name)

  return(compart_df_year_avg)
}

#' Sample Parameter Sets from MCMC Results
#'
#' This function randomly samples a specified number of parameter sets from the MCMC posterior distribution.
#'
#' @param mcmc_results A matrix or data frame containing MCMC posterior samples of parameters.
#' @param num_samples An integer specifying the number of samples to draw.
#' @return A matrix of sampled parameter sets.
#' @examples
#' sampled_params <- sample_mcmc_steps(MCMC_results_2_3$mcmc_run$pars, 100)
#' @export
sample_mcmc_steps <- function(mcmc_results, num_samples) {
  if (num_samples > nrow(mcmc_results)) {
    stop("Number of samples requested exceeds available MCMC steps.")
  }

  sampled_indices <- sample(1:nrow(mcmc_results), num_samples, replace = FALSE)
  return(mcmc_results[sampled_indices, , drop = FALSE])
}

#' Run Model Simulations for Each Parameter Set
#'
#' This function runs the epidemiological model using different parameter sets sampled from the MCMC posterior.
#'
#' @param model A function representing the epidemiological model.
#' @param param_samples A matrix of parameter sets sampled from the MCMC posterior.
#' @param start_date The simulation start date.
#' @param end_date The simulation end date.
#' @param prewarm_years The number of years to prewarm the model.
#' @param days_per_year The number of days in a simulation year.
#' @return A list of data frames containing simulation results for each sampled parameter set.
#' @examples
#' simulations <- run_mcmc_simulations(model, sampled_params, "2014-01-01", "2022-12-31", 2, 360)
#' @export
run_mcmc_simulations <- function(model, param_inputs, param_samples, start_date, end_date, prewarm_years = 2, days_per_year = 360) {
  simulation_results <- lapply(1:nrow(param_samples), function(i) {
    param_sample <- as.list(param_samples[i, ])
    param_inputs <- update_param_list(param_inputs, param_sample)
    compartments_sim(model, param_inputs, start_date, end_date, prewarm_years, days_per_year)
  })

  return(simulation_results)
}

#' Summarize Simulation Results with Confidence Intervals and Selective Variables
#'
#' This function calculates the median and confidence intervals for selected compartments
#' from the simulation results.
#'
#' @param simulation_results A list of data frames containing model outputs for different parameter sets.
#' @param ci_level The confidence level (e.g., 0.95 for 95% confidence intervals).
#' @param variables A character vector specifying which compartments to summarize (e.g., c("SC", "EC", "prev_total_with_R")).
#'                  If NULL, all compartments will be summarized.
#'
#' @return A data frame with median values and confidence intervals for the selected compartments at each time step.
#' @examples
#' # Summarize SC and EC compartments
#' summary_df <- summarize_simulations(simulations, ci_level = 0.95, variables = c("SC", "EC"))
#'
#' # Summarize all compartments
#' summary_df <- summarize_simulations(simulations, ci_level = 0.95)
#' @export
summarize_simulations <- function(simulation_results, ci_level = 0.95, variables = NULL) {
  # Combine all simulation runs into a single data frame
  all_results <- bind_rows(simulation_results, .id = "simulation_id")

  # If variables are NULL, use all columns except 'date' and 'simulation_id'
  if (is.null(variables)) {
    variables <- setdiff(colnames(all_results), c("date", "simulation_id"))
  } else {
    # Ensure selected variables exist in the data
    missing_vars <- setdiff(variables, colnames(all_results))
    if (length(missing_vars) > 0) {
      stop(paste("The following variables are not found in the simulation results:", paste(missing_vars, collapse = ", ")))
    }
  }

  # Group by date and calculate median and CI for selected variables
  summary_stats <- all_results %>%
    group_by(date) %>%
    summarize(across(all_of(variables), list(
      median = ~ median(.x, na.rm = TRUE),
      lower = ~ quantile(.x, probs = (1 - ci_level) / 2, na.rm = TRUE),
      upper = ~ quantile(.x, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")) %>%
    ungroup()

  return(summary_stats)
}


